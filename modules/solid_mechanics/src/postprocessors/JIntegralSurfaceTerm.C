/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "JIntegralSurfaceTerm.h"

template<>
InputParameters validParams<JIntegralSurfaceTerm>()
{
  InputParameters params = validParams<SideIntegralPostprocessor>();
  params.addCoupledVar("q", "The q function, aux variable");
  params.addCoupledVar("disp_x", "The x displacement variable");
  params.addCoupledVar("disp_y", "The y displacement variable");
  params.addCoupledVar("disp_z", "The z displacement variable");
  params.addRequiredParam<UserObjectName>("crack_front_definition","The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index","The index of the point on the crack front corresponding to this q function");
  params.addParam<bool>("convert_J_to_K",false,"Convert J-integral to stress intensity factor K.");
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.addParam<Real>("poissons_ratio","Poisson's ratio");
  params.addParam<Real>("youngs_modulus","Young's modulus of the material.");
  params.addParam<FunctionName>("pressure_function","The function describing the pressure on the crack surface");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

JIntegralSurfaceTerm::JIntegralSurfaceTerm(const std::string & name, InputParameters parameters):
    SideIntegralPostprocessor(name, parameters),
    _scalar_q(coupledValue("q")),
    _grad_disp_x(coupledGradient("disp_x")),
    _grad_disp_y(coupledGradient("disp_y")),
    _grad_disp_z(coupledGradient("disp_z")),
    _pressure_function(getFunction("pressure_function")),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _treat_as_2d(false),
    _convert_J_to_K(getParam<bool>("convert_J_to_K")),
    _has_symmetry_plane(isParamValid("symmetry_plane")),
    _poissons_ratio(isParamValid("poissons_ratio") ? getParam<Real>("poissons_ratio") : 0),
    _youngs_modulus(isParamValid("youngs_modulus") ? getParam<Real>("youngs_modulus") : 0)
{
}

void
JIntegralSurfaceTerm::initialSetup()
{
  _treat_as_2d = _crack_front_definition->treatAs2D();

  if (_treat_as_2d)
  {
    if (_has_crack_front_point_index)
    {
      mooseWarning("crack_front_point_index ignored because CrackFrontDefinition is set to treat as 2D");
    }
  }
  else
  {
    if (!_has_crack_front_point_index)
    {
      mooseError("crack_front_point_index must be specified in JIntegralSurfaceTerm");
    }
  }

  if (_convert_J_to_K && (!isParamValid("youngs_modulus") || !isParamValid("poissons_ratio")))
    mooseError("youngs_modulus and poissons_ratio must be specified if convert_J_to_K = true");
}

Real
JIntegralSurfaceTerm::computeQpIntegral()
{
  const RealVectorValue& crack_direction = _crack_front_definition->getCrackDirection(_crack_front_point_index);

  RealVectorValue grad_disp_in_crack_dir;

  grad_disp_in_crack_dir(0) = _grad_disp_x[_qp](0) * crack_direction(0) +
                              _grad_disp_x[_qp](1) * crack_direction(1) +
                              _grad_disp_x[_qp](2) * crack_direction(2);
  grad_disp_in_crack_dir(1) = _grad_disp_y[_qp](0) * crack_direction(0) +
                              _grad_disp_y[_qp](1) * crack_direction(1) +
                              _grad_disp_y[_qp](2) * crack_direction(2);
  grad_disp_in_crack_dir(2) = _grad_disp_z[_qp](0) * crack_direction(0) +
                              _grad_disp_z[_qp](1) * crack_direction(1) +
                              _grad_disp_z[_qp](2) * crack_direction(2);

  Real q_avg_seg = 1.0;
  if (!_crack_front_definition->treatAs2D())
  {
    q_avg_seg = (_crack_front_definition->getCrackFrontForwardSegmentLength(_crack_front_point_index) +
                 _crack_front_definition->getCrackFrontBackwardSegmentLength(_crack_front_point_index)) / 2.0;
  }

  Real grad_disp_dot_normal = grad_disp_in_crack_dir * _normals[_qp];

  return _pressure_function.value(_t,_q_point[_qp]) * _scalar_q[_qp] * grad_disp_dot_normal / q_avg_seg;
}

Real
JIntegralSurfaceTerm::getValue()
{
  gatherSum(_integral_value);
  if (_has_symmetry_plane)
    _integral_value *= 2.0;

  Real sign = (_integral_value > 0.0) ? 1.0 : ((_integral_value < 0.0) ? -1.0: 0.0);
  if (_convert_J_to_K)
    _integral_value = sign * std::sqrt(std::abs(_integral_value) * _youngs_modulus / (1 - std::pow(_poissons_ratio,2)));

  return _integral_value;
}
