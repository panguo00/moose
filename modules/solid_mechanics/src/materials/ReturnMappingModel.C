/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ReturnMappingModel.h"

#include "SymmIsotropicElasticityTensor.h"
#include "Conversion.h"

template <>
InputParameters
validParams<ReturnMappingModel>()
{
  InputParameters params = validParams<ConstitutiveModel>();

  // Return mapping iteration control parameters
  params.addParam<unsigned int>("max_its", 100, "Maximum number of return mapping iterations");
  params.addParam<bool>(
      "output_iteration_info", false, "Set true to output return mapping iteration information");
  params.addParam<bool>(
      "output_iteration_info_on_error",
      false,
      "Set true to output return mapping iteration information when a step fails");
  params.addParam<Real>(
      "relative_tolerance", 1e-12, "Relative convergence tolerance for return mapping iteration");
  params.addParam<Real>(
      "absolute_tolerance", 1e-15, "Absolute convergence tolerance for return mapping iteration");

  return params;
}

ReturnMappingModel::ReturnMappingModel(const InputParameters & parameters,
                                       const std::string inelastic_strain_name)
  : ConstitutiveModel(parameters),
    _max_its(parameters.get<unsigned int>("max_its")),
    _output_iteration_info(getParam<bool>("output_iteration_info")),
    _output_iteration_info_on_error(getParam<bool>("output_iteration_info_on_error")),
    _relative_tolerance(parameters.get<Real>("relative_tolerance")),
    _absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
    _effective_strain_increment(0),
    _effective_inelastic_strain(
        declareProperty<Real>("effective_" + inelastic_strain_name + "_strain")),
    _effective_inelastic_strain_old(
        declarePropertyOld<Real>("effective_" + inelastic_strain_name + "_strain"))
{
  if (_relative_tolerance > 1e-12)
    mooseWarning("relative_tolerance was set to: ",
                 _relative_tolerance,
                 " for model: ",
                 _name,
                 " Using values greater than the default tolerance (1e-12) is not recommended.");
  if (_absolute_tolerance > 1e-15)
    mooseWarning("absolute_tolerance was set to: ",
                 _absolute_tolerance,
                 " for model: ",
                 _name,
                 " Using values greater than the default tolerance (1e-18) is not recommended.");
  if (_max_its < 100)
    mooseWarning(
        "max_its was set to: " , _max_its , " for model: " , _name,
                               " Using values less than the default (100) is not recommended.");
}

void
ReturnMappingModel::initStatefulProperties(unsigned n_points)
{
  for (unsigned qp(0); qp < n_points; ++qp)
  {
    _effective_inelastic_strain[qp] = 0;
  }
}
void
ReturnMappingModel::computeStress(const Elem & current_elem,
                                  unsigned qp,
                                  const SymmElasticityTensor & elasticityTensor,
                                  const SymmTensor & stress_old,
                                  SymmTensor & strain_increment,
                                  SymmTensor & stress_new)
{
  // Given the stretching, compute the stress increment and add it to the old stress. Also update
  // the creep strain
  // stress = stressOld + stressIncrement
  if (_t_step == 0 && !_app.isRestarting())
    return;

  stress_new = elasticityTensor * strain_increment;
  stress_new += stress_old;

  SymmTensor inelastic_strain_increment;
  computeStress(current_elem,
                qp,
                elasticityTensor,
                stress_old,
                strain_increment,
                stress_new,
                inelastic_strain_increment);
}

void
ReturnMappingModel::computeStress(const Elem & /*current_elem*/,
                                  unsigned qp,
                                  const SymmElasticityTensor & elasticityTensor,
                                  const SymmTensor & stress_old,
                                  SymmTensor & strain_increment,
                                  SymmTensor & stress_new,
                                  SymmTensor & inelastic_strain_increment)
{
  // compute deviatoric trial stress
  SymmTensor dev_trial_stress(stress_new);
  dev_trial_stress.addDiag(-dev_trial_stress.trace() / 3.0);

  // compute effective trial stress
  Real dts_squared = dev_trial_stress.doubleContraction(dev_trial_stress);
  Real effective_trial_stress = std::sqrt(1.5 * dts_squared);

  // compute effective strain increment
  SymmTensor dev_strain_increment(strain_increment);
  dev_strain_increment.addDiag(-strain_increment.trace() / 3.0);
  _effective_strain_increment = dev_strain_increment.doubleContraction(dev_strain_increment);
  _effective_strain_increment = std::sqrt(2.0 / 3.0 * _effective_strain_increment);

  computeStressInitialize(qp, effective_trial_stress, elasticityTensor);

  // Use Newton iterations to determine inelastic strain increment

  Real scalar = 0.0;
  Real scalar_old = 0.0;
  Real scalar_increment = 0.0;
  Real scalar_upper_bound = std::numeric_limits<Real>::max();
  Real scalar_lower_bound = 0.0;
  Real resid_upper_bound = 0.0;
  Real resid_lower_bound = 0.0;
  unsigned int it = 0;
  Real residual_old = 0.0;
  Real reference_residual = 0.0;

  std::stringstream iter_output;

  iterationInitialize(qp, scalar);
  Real residual = computeResidual(qp, effective_trial_stress, scalar, reference_residual);
  iterationFinalize(qp, scalar);
  Real init_resid_sign = (residual < 0.0 ? -1.0 : 1.0);

  output_iter_info(iter_output, it, effective_trial_stress, scalar, residual, reference_residual);

  if (effective_trial_stress != 0.0)
  {
    while (it < _max_its && !converged(residual, reference_residual))
    {
      ++it;
      residual_old = residual;
      scalar_old = scalar;

      scalar_increment = -residual / computeDerivative(qp, effective_trial_stress, scalar);
      scalar = scalar_old + scalar_increment;

      iterationInitialize(qp, scalar);
      residual = computeResidual(qp, effective_trial_stress, scalar, reference_residual);
      iterationFinalize(qp, scalar);

      update_bounds(scalar,
                    residual,
                    init_resid_sign,
                    scalar_upper_bound,
                    scalar_lower_bound,
                    resid_upper_bound,
                    resid_lower_bound,
                    iter_output,
                    false);

      if (converged(residual, reference_residual))
      {
        output_iter_info(
            iter_output, it, effective_trial_stress, scalar, residual, reference_residual);
        break;
      }
      else
      {
        // Line Search
        bool modified_increment = false;
        if (residual_old - residual != 0.0)
        {
          Real alpha = residual_old / (residual_old - residual);
          if (alpha > 1.0) // upper bound for alpha
            alpha = 1.0;
          else if (alpha < 1e-2) // lower bound for alpha
            alpha = 1e-2;
          if (alpha != 1.0)
          {
            modified_increment = true;
            scalar_increment *= alpha;
            if (_output_iteration_info == true || _output_iteration_info_on_error == true)
              iter_output << "  Line search alpha = " << alpha
                          << " increment = " << scalar_increment << std::endl;
          }
        }

        bool bounded = false;
        // Check to see whether trial scalar_increment is outside the bounds, and set it to a point
        // within the bounds if it is
        if (scalar_old + scalar_increment >= scalar_upper_bound ||
            scalar_old + scalar_increment <= scalar_lower_bound)
        {
          if (scalar_upper_bound != std::numeric_limits<Real>::max() && scalar_lower_bound != 0.0)
          {
            // Real frac = resid_upper_bound / (resid_upper_bound - resid_lower_bound);
            Real frac = 0.5;
            scalar_increment =
                (1.0 - frac) * scalar_lower_bound + frac * scalar_upper_bound - scalar_old;
            // scalar_increment = (scalar_upper_bound + scalar_lower_bound) / 2.0 - scalar_old;
            modified_increment = true;
            bounded = true;
            if (_output_iteration_info == true || _output_iteration_info_on_error == true)
              iter_output << "  Trial scalar_increment exceeded bounds.  Setting between "
                             "lower/upper bounds. frac: "
                          << frac << std::endl;
          }
        }

        // Update the trial scalar and recompute residual if the line search or bounds checking
        // modified the increment
        if (modified_increment)
        {
          scalar = scalar_old + scalar_increment;
          iterationInitialize(qp, scalar);
          residual = computeResidual(qp, effective_trial_stress, scalar, reference_residual);
          iterationFinalize(qp, scalar);

          update_bounds(scalar,
                        residual,
                        init_resid_sign,
                        scalar_upper_bound,
                        scalar_lower_bound,
                        resid_upper_bound,
                        resid_lower_bound,
                        iter_output,
                        bounded);
        }

        //        iter_output << "BWS eps factor1 = "
        //                    << (scalar_upper_bound - scalar_lower_bound) / (PETSC_MACHINE_EPSILON
        //                    * scalar_upper_bound)
        //                    << " eps factor2 = "
        //                    << std::abs(scalar_increment) / (scalar_old  * PETSC_MACHINE_EPSILON)
        //                    << std::endl;

        // If the difference between the upper and lower bounds is too close to machine epsilon
        // or the scalar increment is too close to machine epsilon, break
        //        if (converged(1e3 * residual, reference_residual) &&
        //            ((scalar_lower_bound > 0.0 &&
        //              (scalar_upper_bound - scalar_lower_bound) / (PETSC_MACHINE_EPSILON *
        //              scalar_lower_bound) <= 10.0) ||
        //             (residual != 0.0 &&
        //              scalar_old != 0.0 &&
        //              std::abs(scalar_increment) < 10.0 * PETSC_MACHINE_EPSILON * scalar_old)))
        //        {
        //          output_iter_info(iter_output, it, effective_trial_stress, scalar, residual,
        //          reference_residual);
        //          break;
        //        }
      }

      output_iter_info(
          iter_output, it, effective_trial_stress, scalar, residual, reference_residual);
    }

    if (_output_iteration_info)
      _console << iter_output.str();

    if (residual != residual)
    {
      if (_output_iteration_info_on_error)
        Moose::err << iter_output.str();
      mooseError("Encountered nan in material: ",
                 _name,
                 ".  Rerun with  'output_iteration_info_on_error = true' for more information.");
    }

    if (it == _max_its && !converged(residual, reference_residual))
    {
      if (_output_iteration_info_on_error)
        Moose::err << iter_output.str();
      mooseError("Exceeded maximum iterations in ReturnMappingModel solve for material: ",
                 _name,
                 ".  Rerun with  'output_iteration_info_on_error = true' for more information.");
    }

    // compute inelastic and elastic strain increments
    inelastic_strain_increment = dev_trial_stress * (1.5 * scalar / effective_trial_stress);
    strain_increment -= inelastic_strain_increment;
  }
  else
    inelastic_strain_increment.zero();

  // compute stress increment
  stress_new = elasticityTensor * strain_increment;

  // update stress
  stress_new += stress_old;

  computeStressFinalize(qp, inelastic_strain_increment);
}

bool
ReturnMappingModel::converged(const Real & residual, const Real & reference)
{
  return (std::abs(residual) <= _absolute_tolerance ||
          (std::abs(residual) / reference) <= _relative_tolerance);
}

void
ReturnMappingModel::output_iter_info(std::stringstream & iter_output,
                                     const unsigned int & it,
                                     const Real & effective_trial_stress,
                                     const Real & scalar,
                                     const Real & residual,
                                     const Real & reference_residual)
{
  if (_output_iteration_info == true || _output_iteration_info_on_error == true)
  {
    iter_output << " it=" << it << " trl_strs=" << effective_trial_stress << " scalar=" << scalar
                << " residual=" << residual
                << " ref_res=" << reference_residual
                << " rel_res=" << std::abs(residual) / reference_residual
                << " rel_tol=" << _relative_tolerance << " abs_res=" << std::abs(residual)
                << " abs_tol=" << _absolute_tolerance << std::endl;
  }
}

void
ReturnMappingModel::update_bounds(const Real & scalar,
                                  const Real & residual,
                                  const Real & init_resid_sign,
                                  Real & scalar_upper_bound,
                                  Real & scalar_lower_bound,
                                  Real & resid_upper_bound,
                                  Real & resid_lower_bound,
                                  std::stringstream & iter_output,
                                  bool bounded)
{
  // Update upper/lower bounds as applicable
  if (residual * init_resid_sign < 0.0 && scalar < scalar_upper_bound)
  {
    scalar_upper_bound = scalar;
    resid_upper_bound = residual;
    if (scalar_upper_bound < scalar_lower_bound)
    {
      scalar_upper_bound = scalar_lower_bound;
      resid_upper_bound = resid_lower_bound;
      scalar_lower_bound = 0.0;
      resid_lower_bound = 0.0;
      if (_output_iteration_info == true || _output_iteration_info_on_error == true)
        iter_output << "  Corrected for scalar_upper_bound < scalar_lower_bound" << std::endl;
    }
  }
  // Don't permit setting scalar_lower_bound > scalar_upper_bound (but do permit the reverse).
  // This ensures that if we encounter multiple roots, we pick the lowest one.
  else if (residual * init_resid_sign > 0.0 && scalar > scalar_lower_bound &&
           scalar < scalar_upper_bound)
  {
    scalar_lower_bound = scalar;
    resid_lower_bound = residual;
  }

  if ((scalar_lower_bound > 0 &&
       (scalar_upper_bound - scalar_lower_bound) / (PETSC_MACHINE_EPSILON * scalar_lower_bound) <=
           10.0))
  {
    // Reset lower/upper bounds.  They're too close together and solution hasn't converged, so they
    // must be bad.
    scalar_lower_bound = 0.0;
    resid_lower_bound = 0.0;
    scalar_upper_bound = std::numeric_limits<Real>::max();
    resid_upper_bound = 0.0;
    if (_output_iteration_info == true || _output_iteration_info_on_error == true)
      iter_output << "  Reset lower/upper bounds because they're too close together" << std::endl;
  }
//  else if (bounded && (std::abs(residual) > 0.9 * std::abs(resid_lower_bound) ||
//                       std::abs(residual) > 0.9 * std::abs(resid_upper_bound)))
//  {
//    // Reset lower/upper bounds. The residual didn't decrease enough.
//    scalar_lower_bound = 0.0;
//    resid_lower_bound = 0.0;
//    scalar_upper_bound = std::numeric_limits<Real>::max();
//    resid_upper_bound = 0.0;
//    iter_output << "  Reset lower/upper bounds because residual didn't decrease enough"
//                << std::endl;
//  }
}
