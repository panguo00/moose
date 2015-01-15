#include "TrussMaterial.h"

#include "Material.h"
#include "ColumnMajorMatrix.h"
#include "SymmIsotropicElasticityTensor.h"
#include "VolumetricModel.h"

template<>
InputParameters validParams<TrussMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("disp_x", "Variable containing the x displacement");
  params.addCoupledVar("disp_y", 0.0, "Variable containing the y displacement");
  params.addCoupledVar("disp_z", 0.0, "Variable containing the z displacement");
  params.addParam<Real>("youngs_modulus", "Young's Modulus");
  params.addCoupledVar("youngs_modulus_var","Variable containing Young's modulus");
  params.addParam<Real>("t_ref", 0.0, "The reference temperature at which this material has zero strain.");
  params.addParam<Real>("thermal_expansion", 0.0, "The thermal expansion coefficient.");
  params.addCoupledVar("temp", "The temperature if you want thermal expansion.");
  return params;
}

TrussMaterial::TrussMaterial(const std::string  & name,
                             InputParameters parameters)
  :Material(name, parameters),
   _disp_x_nodal_val(coupledValue("disp_x")),
   _disp_y_nodal_val(coupledValue("disp_y")),
   _disp_z_nodal_val(coupledValue("disp_z")),
   _axial_stress(declareProperty<Real>("axial_stress")),
   _e_over_l(declareProperty<Real>("e_over_l")),
   _youngs_modulus(isParamValid("youngs_modulus") ? getParam<Real>("youngs_modulus") : 0),
   _youngs_modulus_coupled(isCoupled("youngs_modulus_var")),
   _youngs_modulus_var(_youngs_modulus_coupled ? coupledValue("youngs_modulus_var"): _zero),
   _has_temp(isCoupled("temp")),
   _temp(_has_temp ? coupledValue("temp") : _zero),
   _t_ref(getParam<Real>("t_ref")),
   _alpha(getParam<Real>("thermal_expansion")),
   _dim(1)
{
// This doesn't work: if there are just line elements, this
// returns 1 even if the mesh has higher dimensionality
//  unsigned int dim = _subproblem.mesh().dimension();


  if (isCoupled("disp_y"))
  {
    _dim = 2;
    if (isCoupled("disp_z"))
      _dim = 3;
  }

  if (parameters.isParamValid("youngs_modulus"))
  {
    if (_youngs_modulus_coupled)
      mooseError("Cannot specify both youngs_modulus and youngs_modulus_var");
  }
  else
  {
    if (!_youngs_modulus_coupled)
      mooseError("Must specify either youngs_modulus or youngs_modulus_var");
  }
}

TrussMaterial::~TrussMaterial()
{
}

void
TrussMaterial::computeProperties()
{
  const Node* const node0=_current_elem->get_node(0);
  const Node* const node1=_current_elem->get_node(1);

  Real dx=(*node1)(0)-(*node0)(0);
  Real dy=0;
  Real dz=0;
  if (_dim > 1)
  {
    dy=(*node1)(1)-(*node0)(1);
    if (_dim > 2)
      dz=(*node1)(2)-(*node0)(2);
  }
  Real orig_length=std::sqrt( dx*dx + dy*dy + dz*dz );

  RealVectorValue disp_vec0(_disp_x_nodal_val[0], _disp_y_nodal_val[0], _disp_z_nodal_val[0]);
  RealVectorValue disp_vec1(_disp_x_nodal_val[1], _disp_y_nodal_val[1], _disp_z_nodal_val[1]);
  std::cout<<"BWS dim: "<<_dim<<std::endl;
  std::cout<<"BWS disp_vec0: "<<disp_vec0(0)<<" "<<disp_vec0(1)<<" "<<disp_vec0(2)<<std::endl;
  std::cout<<"BWS disp_vec1: "<<disp_vec1(0)<<" "<<disp_vec1(1)<<" "<<disp_vec1(2)<<std::endl;

  Real ddx=dx+disp_vec1(0)-disp_vec0(0);
  Real ddy=0;
  Real ddz=0;
  if (_dim > 1)
  {
    ddy=dy+disp_vec1(1)-disp_vec0(1);
    if (_dim > 2)
    {
      ddz=dz+disp_vec1(2)-disp_vec0(2);
    }
  }
  Real new_length = std::sqrt( ddx*ddx + ddy*ddy + ddz*ddz );
  Real strain = (new_length-orig_length)/orig_length;

  Real thermal_strain = 0.0;

  for (_qp=0; _qp < _qrule->n_points(); ++_qp)
  {
    Real youngs_modulus(_youngs_modulus_coupled ? _youngs_modulus_var[_qp] : _youngs_modulus);
    if (_has_temp)
    {
      thermal_strain = _alpha * (_t_ref - _temp[_qp]);
    }
    _axial_stress[_qp] = youngs_modulus*(strain+thermal_strain);
    _e_over_l[_qp] = youngs_modulus/orig_length;
  }
}
