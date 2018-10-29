//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VolumeWeightedWeibullIC.h"

#include "libmesh/point.h"
#include "Distribution.h"

registerMooseObject("TensorMechanicsApp", VolumeWeightedWeibullIC);

template <>
InputParameters
validParams<VolumeWeightedWeibullIC>()
{
  InputParameters params = validParams<RandomICBase>();
  params.addRequiredParam<Real>("reference_volume", "Reference volume (of a test specimen)");
  params.addRequiredParam<Real>("weibull_modulus", "Weibull modulus");
  params.addParam<Real>("median", "Median value of property measured in a specimen of volume equal to reference_volume");
  params.addClassDescription("Initialize a variable with randomly generated numbers following "
                             "a volume-weighted Weibull distribution");
  return params;
}

VolumeWeightedWeibullIC::VolumeWeightedWeibullIC(const InputParameters & parameters)
  : RandomICBase(parameters),
  _reference_volume(getParam<Real>("reference_volume")),
  _weibull_modulus(getParam<Real>("weibull_modulus")),
  _median(getParam<Real>("median"))
{
}

Real
VolumeWeightedWeibullIC::value(const Point & /*p*/)
{
  const Real & element_volume = _current_elem->volume();

  return _median*std::pow((_reference_volume*std::log(generateRandom()) / element_volume * std::log(0.5)), 1.0 / _weibull_modulus);
}
