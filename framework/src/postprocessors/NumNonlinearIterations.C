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

#include "NumNonlinearIterations.h"

#include "FEProblem.h"
#include "SubProblem.h"

template<>
InputParameters validParams<NumNonlinearIterations>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<bool>("accumulate_over_step", false, "When set to true, accumulates to count the total over all Picard iterations for each step");
  return params;
}

NumNonlinearIterations::NumNonlinearIterations(const std::string & name, InputParameters parameters) :
    GeneralPostprocessor(name, parameters),
    _feproblem(dynamic_cast<FEProblem &>(_subproblem)),
    _accumulate_over_step(getParam<bool>("accumulate_over_step")),
    _time(-std::numeric_limits<Real>::max()),
    _num_iters(0)
{}

Real
NumNonlinearIterations::getValue()
{
  if (_accumulate_over_step)
  {
    if (_feproblem.time() != _time)
    {
      _num_iters = 0;
      _time = _feproblem.time();
    }
    _num_iters += _subproblem.nNonlinearIterations();
  }
  else
    _num_iters = _subproblem.nNonlinearIterations();

  return _num_iters;
}
