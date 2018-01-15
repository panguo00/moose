/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef JINTEGRAL_H
#define JINTEGRAL_H

#include "ElementVectorPostprocessor.h"
#include "CrackFrontDefinition.h"

// Forward Declarations
class JIntegral;
class RankTwoTensor;

template <>
InputParameters validParams<JIntegral>();

/**
 * This vectorpostprocessor computes the J-Integral
 *
 */
class JIntegral : public ElementVectorPostprocessor
{
public:
  JIntegral(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  Real computeQpIntegral(const unsigned int crack_front_point_index);
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  bool _treat_as_2d;
  const MaterialProperty<RankTwoTensor> & _Eshelby_tensor;
  const MaterialProperty<RealVectorValue> * _J_thermal_term_vec;
  bool _convert_J_to_K;
  bool _has_symmetry_plane;
  Real _poissons_ratio;
  Real _youngs_modulus;
  unsigned int _ring_index;
  MooseEnum _q_function_type;
  std::vector<Real> _q_curr_elem;
  const std::vector<std::vector<Real>> * _phi_curr_elem;
  const std::vector<std::vector<RealGradient>> * _dphi_curr_elem;

  VectorPostprocessorValue & _j_integral;
};

#endif // JINTEGRAL3D_H
