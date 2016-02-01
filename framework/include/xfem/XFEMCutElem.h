/****************************************************************/
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

#ifndef XFEMCUTELEM_H
#define XFEMCUTELEM_H

// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/location_maps.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/vector_value.h"
#include "libmesh/point_locator_base.h"

// C++ Includes   -----------------------------------
#include <vector>

using namespace libMesh;

class EFANode;
class EFAElement;

class XFEMCutElem
{
public:

  XFEMCutElem(Elem* elem, unsigned int n_qpoints);
  virtual ~XFEMCutElem();

protected:

  unsigned int _n_nodes;
  unsigned int _n_qpoints;
  std::vector<Node*> _nodes;
  std::vector<Point> _qp_points;
  std::vector<Real> _qp_weights;
  Real _elem_volume;
  Real _physical_volfrac;
  std::vector<Real> _new_weights; // quadrature weights from moment fitting
  virtual Point getNodeCoordinates(EFANode* node, MeshBase* displaced_mesh = NULL) const = 0;

public:

  void setQuadraturePointsAndWeights(const std::vector<Point> &qp_points, const std::vector<Real> &qp_weights);
  virtual void computePhysicalVolumeFraction() = 0;
  Real getPhysicalVolumeFraction() const;
  virtual void computeMomentFittingWeights() = 0;
  Real getMomentFittingWeight(unsigned int i_qp) const;
  virtual Point getCutPlaneOrigin(unsigned int plane_id, MeshBase* displaced_mesh=NULL) const = 0;
  virtual Point getCutPlaneNormal(unsigned int plane_id, MeshBase* displaced_mesh=NULL) const = 0;
  virtual void getCrackTipOriginAndDirection(unsigned tip_id, Point & origin, Point & direction) const = 0;
  virtual void getFragmentFaces(std::vector<std::vector<Point> > &frag_faces, MeshBase* displaced_mesh=NULL) const = 0;
  virtual const EFAElement * getEFAElement() const = 0;
  virtual unsigned int numCutPlanes() const = 0;
};
#endif
