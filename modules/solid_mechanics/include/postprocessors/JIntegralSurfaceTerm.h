/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef JINTEGRALSURFACETERM_H
#define JINTEGRALSURFACETERM_H

#include "SideIntegralPostprocessor.h"
#include "CrackFrontDefinition.h"
#include "Function.h"

//Forward Declarations
class JIntegralSurfaceTerm;

template<>
InputParameters validParams<JIntegralSurfaceTerm>();

/**
 * This postprocessor computes the surface term of the J-Integral
 *
 */
class JIntegralSurfaceTerm: public SideIntegralPostprocessor
{
public:
  JIntegralSurfaceTerm(const std::string & name, InputParameters parameters);
  virtual Real getValue();

protected:
  virtual void initialSetup();
  virtual Real computeQpIntegral();
  VariableValue & _scalar_q;
  VariableGradient & _grad_disp_x;
  VariableGradient & _grad_disp_y;
  VariableGradient & _grad_disp_z;
  Function & _pressure_function;
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  bool _treat_as_2d;
  bool _convert_J_to_K;
  bool _has_symmetry_plane;
  Real _poissons_ratio;
  Real _youngs_modulus;

};

#endif //JINTEGRALSURFACETERM_H
