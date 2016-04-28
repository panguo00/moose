/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ReturnMappingModel.h"

#include "SymmIsotropicElasticityTensor.h"

template<>
InputParameters validParams<ReturnMappingModel>()
{
  InputParameters params = validParams<ConstitutiveModel>();

  // Return mapping iteration control parameters
  params.addParam<unsigned int>("max_its", 30, "Maximum number of return mapping iterations");
  params.addParam<bool>("output_iteration_info", false, "Set true to output return mapping iteration information");
   params.addParam<bool>("output_iteration_info_on_error", false, "Set true to output return mapping iteration information when a step fails");
  params.addParam<Real>("relative_tolerance", 1e-12, "Relative convergence tolerance for return mapping iteration");
  params.addParam<Real>("absolute_tolerance", 1e-15, "Absolute convergence tolerance for return mapping iteration");

  return params;
}


ReturnMappingModel::ReturnMappingModel(const InputParameters & parameters) :
    ConstitutiveModel(parameters),
    _max_its(parameters.get<unsigned int>("max_its")),
    _output_iteration_info(getParam<bool>("output_iteration_info")),
    _output_iteration_info_on_error(getParam<bool>("output_iteration_info_on_error")),
    _relative_tolerance(parameters.get<Real>("relative_tolerance")),
    _absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
    _effective_strain_increment(0)
{
  if (_relative_tolerance > 1e-12)
    mooseWarning("relative_tolerance was set to: " << _relative_tolerance << " for model: " << _name << " Using values greater than the default tolerance (1e-12) is not recommended.");
  if (_absolute_tolerance > 1e-15)
    mooseWarning("absolute_tolerance was set to: " << _absolute_tolerance << " for model: " << _name << " Using values greater than the default tolerance (1e-15) is not recommended.");
}


void
ReturnMappingModel::computeStress(const Elem & current_elem,
                                  unsigned qp, const SymmElasticityTensor & elasticityTensor,
                                  const SymmTensor & stress_old, SymmTensor & strain_increment,
                                  SymmTensor & stress_new)
{
  // Given the stretching, compute the stress increment and add it to the old stress. Also update the creep strain
  // stress = stressOld + stressIncrement
  if (_t_step == 0) return;

  stress_new = elasticityTensor * strain_increment;
  stress_new += stress_old;

  SymmTensor inelastic_strain_increment;
  computeStress(current_elem, qp, elasticityTensor, stress_old,
                 strain_increment, stress_new, inelastic_strain_increment);
}

void
ReturnMappingModel::computeStress(const Elem & /*current_elem*/, unsigned qp,
                                  const SymmElasticityTensor & elasticityTensor,
                                  const SymmTensor & stress_old,
                                  SymmTensor & strain_increment,
                                  SymmTensor & stress_new,
                                  SymmTensor & inelastic_strain_increment)
{
  // compute deviatoric trial stress
  SymmTensor dev_trial_stress(stress_new);
  dev_trial_stress.addDiag(-dev_trial_stress.trace()/3.0);

  // compute effective trial stress
  Real dts_squared = dev_trial_stress.doubleContraction(dev_trial_stress);
  Real effective_trial_stress = std::sqrt(1.5 * dts_squared);

  // compute effective strain increment
  SymmTensor dev_strain_increment(strain_increment);
  dev_strain_increment.addDiag(-strain_increment.trace()/3.0);
  _effective_strain_increment = dev_strain_increment.doubleContraction(dev_strain_increment);
  _effective_strain_increment = std::sqrt(2.0/3.0 * _effective_strain_increment);

  computeStressInitialize(qp, effective_trial_stress, elasticityTensor);

  // Use Newton iterations to determine inelastic strain increment

  Real scalar = 0.0;
  Real scalar_old = 0.0;
  Real scalar_increment = 0.0;
  Real scalar_upper_bound = std::numeric_limits<Real>::max();
  Real scalar_lower_bound = 0.0;
  unsigned int it = 0;
  Real residual = 0.0;
  Real residual_old = 0.0;
  Real norm_residual = 10;
  Real first_norm_residual = 10;
  Real init_resid_sign = 1.0;

  std::stringstream iter_output;

  while (it < _max_its &&
        norm_residual > _absolute_tolerance &&
        (norm_residual / first_norm_residual) > _relative_tolerance)
  {
    iterationInitialize(qp, scalar);
    residual = computeResidual(qp, effective_trial_stress, scalar);
    if (it == 0)
    {
      if (residual < 0.0)
        init_resid_sign = -1.0;
    }
    else
    {
      // Update upper/lower bounds as applicable
      if (residual*init_resid_sign < 0.0 &&
          scalar < scalar_upper_bound)
        scalar_upper_bound = scalar;
      else if (residual*init_resid_sign > 0.0 &&
               scalar > scalar_lower_bound)
        scalar_lower_bound = scalar;
      if (scalar_upper_bound < scalar_lower_bound)
        mooseError("scalar_upper_bound < scalar_lower_bound in return mapping iteration. This is likely an issue with the derived material model");

      // Line Search
      bool modified_increment = false;
      if (residual_old - residual != 0.0)
      {
        Real alpha = residual_old / (residual_old - residual);
        if (alpha > 1.0)
          alpha = 1.0;
        else if (alpha < 0.0)
          alpha = 0.0;
        modified_increment = true;
        scalar_increment *= alpha;
      }

      // Check to see whether trial scalar_increment is outside bounds, and set it to bisection point of bounds if it is
      if (scalar_old + scalar_increment > scalar_upper_bound ||
          scalar_old + scalar_increment < scalar_lower_bound)
      {
        scalar_increment = (scalar_upper_bound + scalar_lower_bound) / 2.0 - scalar_old;
        modified_increment = true;

        //If the difference between the upper and lower bounds is too close to machine epsilon, break
        if (scalar_lower_bound > 0 &&
            (scalar_upper_bound - scalar_lower_bound) / (PETSC_MACHINE_EPSILON * scalar_lower_bound) < 5.0)
          break;
      }

      // Update the trial scalar and recompute residual if the line search or bounds checking modified the increment
      if (modified_increment)
      {
        scalar = scalar_old + scalar_increment;
        iterationInitialize(qp, scalar);
        residual = computeResidual(qp, effective_trial_stress, scalar);
        if (residual * init_resid_sign < 0.0 &&
            scalar < scalar_upper_bound)
          scalar_upper_bound = scalar;
        else if (residual * init_resid_sign > 0.0 &&
                 scalar > scalar_lower_bound)
          scalar_lower_bound = scalar;
        if (scalar_upper_bound < scalar_lower_bound)
          mooseError("scalar_upper_bound < scalar_lower_bound in return mapping iteration. This is likely an issue with the derived material model");
      }
    }

    norm_residual = std::abs(residual);
    if (it == 0)
    {
      first_norm_residual = norm_residual;
      if (first_norm_residual == 0)
        first_norm_residual = 1;
    }

    scalar_increment = -residual / computeDerivative(qp, effective_trial_stress, scalar);

    if (_output_iteration_info == true ||
        _output_iteration_info_on_error == true)
    {
      iter_output
        << " it="       << it
        << " trl_strs=" << effective_trial_stress
        << " scalar="   << scalar
        << " residual="  << residual
        << " rel_res="  << norm_residual / first_norm_residual
        << " rel_tol="  << _relative_tolerance
        << " abs_res="  << norm_residual
        << " abs_tol="  << _absolute_tolerance
        << std::endl;
    }

    iterationFinalize(qp, scalar);
    residual_old = residual;
    scalar_old = scalar;
    scalar = scalar_old + scalar_increment;

    ++it;
  }

  if (_output_iteration_info)
    _console << iter_output.str();


  if (it == _max_its &&
     norm_residual > _absolute_tolerance &&
     (norm_residual/first_norm_residual) > _relative_tolerance)
  {
    if (_output_iteration_info_on_error)
    {
      Moose::err << iter_output.str();
    }
    mooseError("Exceeded maximum iterations in ReturnMappingModel solve for material: " << _name << ".  Rerun with  'output_iteration_info_on_error = true' for more information.");
  }

  // compute inelastic and elastic strain increments (avoid potential divide by zero - how should this be done)?
  if (effective_trial_stress < 0.01)
  {
    effective_trial_stress = 0.01;
  }

  inelastic_strain_increment = dev_trial_stress;
  inelastic_strain_increment *= (1.5*scalar/effective_trial_stress);

  strain_increment -= inelastic_strain_increment;

  // compute stress increment
  stress_new = elasticityTensor * strain_increment;

  // update stress
  stress_new += stress_old;

  computeStressFinalize(qp, inelastic_strain_increment);
}
