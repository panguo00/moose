/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef WEIBULLRANDOMIC_H
#define WEIBULLRANDOMIC_H

#include "InitialCondition.h"

// System includes
#include <string>

// Forward Declarations
class InputParameters;
class WeibullRandomIC;
namespace libMesh
{
class Point;
}

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<WeibullRandomIC>();

/**
 * WeibullRandomIC just returns a Random value.
 */
class WeibullRandomIC : public InitialCondition
{
public:
  /**
   * Constructor
   *
   * @param parameters The parameters object holding data for the class to use.
   */
  WeibullRandomIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  Real _reference_volume;
  Real _weibull_modulus;
  Real _median;
};

#endif // RANDOMIC_H
