/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "WeibullRandomIC.h"
#include "MooseRandom.h"

#include "libmesh/point.h"


template <>
InputParameters
validParams<WeibullRandomIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("reference_volume", "Reference volume (of a test specimen)");
  params.addRequiredParam<Real>("weibull_modulus", "Weibull modulus");
  params.addParam<Real>("median", "Median value of property measured in a specimen of volume equal to reference_volume");
  return params;
}

WeibullRandomIC::WeibullRandomIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _reference_volume(getParam<Real>("reference_volume")),
    _weibull_modulus(getParam<Real>("weibull_modulus")),
    _median(getParam<Real>("median"))
{
}

Real
WeibullRandomIC::value(const Point & /*p*/)
{
  // Uniform distributed random number between 0 and 1
  const Real rand_num = MooseRandom::rand();

  // Volume of the current element
  const Real & element_volume = _current_elem->volume();

  Real sampled_value = _median*std::pow((_reference_volume*std::log(rand_num) / element_volume * std::log(0.5)), 1.0 / _weibull_modulus);

  return sampled_value;
}
