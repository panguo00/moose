/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMTMMATERIALTENSORMARKERUSEROBJECT_H
#define XFEMTMMATERIALTENSORMARKERUSEROBJECT_H

#include "XFEMMarkerUserObject.h"

class XFEMTMMaterialTensorMarkerUserObject;
class RankTwoTensor;

template <>
InputParameters validParams<XFEMTMMaterialTensorMarkerUserObject>();

class XFEMTMMaterialTensorMarkerUserObject : public XFEMMarkerUserObject
{
public:
  XFEMTMMaterialTensorMarkerUserObject(const InputParameters & parameters);
  virtual ~XFEMTMMaterialTensorMarkerUserObject() {}

protected:
  const MaterialProperty<RankTwoTensor> & _tensor;
  MooseEnum _scalar_type;
  const Point _point1;
  const Point _point2;
  Real _threshold;
  bool _average;
  Real _random_range;

  virtual bool doesElementCrack(RealVectorValue & direction);
};

#endif // XFEMTMMATERIALTENSORMARKERUSEROBJECT_H
