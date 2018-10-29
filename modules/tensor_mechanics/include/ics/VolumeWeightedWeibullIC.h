//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef VOLUMEWEIGHTEDWEIBULLIC_H
#define VOLUMEWEIGHTEDWEIBULLIC_H

#include "RandomICBase.h"

// Forward Declarations
class InputParameters;
class VolumeWeightedWeibullIC;
class Distribution;
namespace libMesh
{
class Point;
}

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<VolumeWeightedWeibullIC>();

/**
 * VolumeWeightedWeibullIC just returns a Random value.
 */
class VolumeWeightedWeibullIC : public RandomICBase
{
public:
  /**
   * Constructor
   * @param parameters The parameters object holding data for the class to use.
   */
  VolumeWeightedWeibullIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  ///
  const Real _reference_volume;
  const Real _weibull_modulus;
  const Real _median;
};

#endif // VOLUMEWEIGHTEDWEIBULLIC_H
