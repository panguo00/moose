/***************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "XFEMCutElem.h"

#include "EFANode.h"
#include "EFAElement.h"

XFEMCutElem::XFEMCutElem(Elem* elem, unsigned int n_qpoints):
  _n_nodes(elem->n_nodes()),
  _n_qpoints(n_qpoints),
  _nodes(_n_nodes,NULL)
{
  for (unsigned int i = 0; i < _n_nodes; ++i)
    _nodes[i] = elem->get_node(i);
  _elem_volume = elem->volume();
}

XFEMCutElem::~XFEMCutElem()
{
}

Real
XFEMCutElem::getPhysicalVolumeFraction() const
{
  return _physical_volfrac;
}

Real
XFEMCutElem::getMomentFittingWeight(unsigned int i_qp) const
{
  return _new_weights[i_qp];
}

void XFEMCutElem::setQuadraturePointsAndWeights(const std::vector<Point> &qp_points, const std::vector<Real> &qp_weights)
{
  _qp_points = qp_points;
  _qp_weights = qp_weights;
}

