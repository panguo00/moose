/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SimplePredictor.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<SimplePredictor>()
{
  InputParameters params = validParams<Predictor>();
  params.addRequiredParam<Real>("scale", "The scale factor for the predictor (can range from 0 to 1)");
  params.addParam<std::vector<Real> >("skip_times", "Time steps to skip");

  return params;
}

SimplePredictor::SimplePredictor(const InputParameters & parameters) :
    Predictor(parameters),
    _scale(getParam<Real>("scale")),
    _skip_times(getParam<std::vector<Real> >("skip_times"))
{
}

SimplePredictor::~SimplePredictor()
{
}

void
SimplePredictor::apply(NumericVector<Number> & sln)
{
  if (shouldApplyThisStep())
  {
    // Save the original stream flags
    std::ios_base::fmtflags out_flags = Moose::out.flags();

    _console << "  Applying predictor with scale factor = " << std::fixed << std::setprecision(2) << _scale << std::endl;

    // Restore the flags
    Moose::out.flags(out_flags);

    Real dt_adjusted_scale_factor = _scale * _dt / _dt_old;
    if (dt_adjusted_scale_factor != 0.0)
    {
      sln *= (1.0 + dt_adjusted_scale_factor);
      sln.add(-dt_adjusted_scale_factor, _solution_older);
    }
  }
  else
    _console << "  Skipping predictor this step" << std::endl;
}

bool
SimplePredictor::shouldApplyThisStep()
{
  bool should_apply = true;

  if (_dt_old <= 0)
    should_apply = false;

  if (should_apply)
  {
    const Real & current_time =  _fe_problem.time();
    for (unsigned int i=0; i<_skip_times.size(); ++i)
    {
      if (MooseUtils::absoluteFuzzyEqual(current_time, _skip_times[i]))
      {
        should_apply = false;
        break;
      }
    }
  }
  return should_apply;
}
