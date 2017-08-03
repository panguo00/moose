/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMTMMaterialTensorMarkerUserObject.h"

#include "libmesh/quadrature.h"
#include "RankTwoTensor.h"
#include "RankTwoScalarTools.h"

template <>
InputParameters
validParams<XFEMTMMaterialTensorMarkerUserObject>()
{
  InputParameters params = validParams<XFEMMarkerUserObject>();
  params.addParam<MooseEnum>("scalar_type", RankTwoScalarTools::scalarOptions(), "Scalar quantity to be computed from tensor and used as a failure criterion");
  params.addRequiredParam<std::string>("tensor", "The material tensor name.");
  params.addRequiredParam<Real>("threshold", "The threshold for crack growth.");
  params.addRequiredParam<bool>(
      "average", "Should the tensor quantity be averaged over the quadrature points?");
  params.addParam<Real>(
      "random_range", 0.0, "Range of a uniform random distribution for the threshold");
  params.addParam<Point>("point1", Point(0, 0, 0), "Start point for axis used to calculate some cylindrical material tensor quantities");
  params.addParam<Point>("point2", Point(0, 1, 0), "End point for axis used to calculate some cylindrical material tensor quantities");
  return params;
}

XFEMTMMaterialTensorMarkerUserObject::XFEMTMMaterialTensorMarkerUserObject(
    const InputParameters & parameters)
  : XFEMMarkerUserObject(parameters),
    _tensor(getMaterialProperty<RankTwoTensor>(getParam<std::string>("tensor"))),
    _scalar_type(getParam<MooseEnum>("scalar_type")),
    _point1(parameters.get<Point>("point1")),
    _point2(parameters.get<Point>("point2")),
    _threshold(getParam<Real>("threshold")),
    _average(getParam<bool>("average")),
    _random_range(getParam<Real>("random_range"))
{
  setRandomResetFrequency(EXEC_INITIAL);
}

bool
XFEMTMMaterialTensorMarkerUserObject::doesElementCrack(RealVectorValue & direction)
{
  bool does_it_crack = false;
  unsigned int numqp = _qrule->n_points();

  Real rnd_mult = (1.0 - _random_range / 2.0) + _random_range * getRandomReal();

  if (_average)
  {
    RankTwoTensor average_tensor;
    Point average_point;
    for (unsigned int qp = 0; qp < numqp; ++qp)
    {
      average_tensor += _tensor[qp];
      average_point += _q_point[qp];
    }
    average_tensor *= 1.0 / (Real)numqp;
    average_point *= 1.0 / (Real)numqp;
    Point point_dir;
    Real tensor_quantity = RankTwoScalarTools::getQuantity(average_tensor, _scalar_type, _point1, _point2, average_point, point_dir);
    direction(0) = point_dir(0);
    direction(1) = point_dir(1);
    direction(2) = point_dir(2);
    if (tensor_quantity > _threshold * rnd_mult)
      does_it_crack = true;
  }
  else
  {
    unsigned int max_index = 999999;
    std::vector<Real> tensor_quantities;
    tensor_quantities.reserve(numqp);
    Real max_quantity = 0;
    std::vector<Point> directions;
    directions.resize(numqp);
    for (unsigned int qp = 0; qp < numqp; ++qp)
    {
      tensor_quantities[qp] = RankTwoScalarTools::getQuantity(_tensor[qp], _scalar_type, _point1, _point2, _q_point[qp], directions[qp]);
      if (directions[qp](0) == 0 && directions[qp](1) == 0 && directions[qp](2) == 0)
      {
        mooseError("Direction has zero length in XFEMTMMaterialTensorMarkerUserObject");
      }
      if (tensor_quantities[qp] > max_quantity)
      {
        max_quantity = tensor_quantities[qp];
        max_index = qp;
      }
    }
    if (max_quantity > _threshold * rnd_mult)
    {
      does_it_crack = true;
      direction(0) = directions[max_index](0);
      direction(1) = directions[max_index](1);
      direction(2) = directions[max_index](2);
//      direction = directions[max_index];
    }
  }

  return does_it_crack;
}
