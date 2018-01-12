/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef GEOMETRICCUTUSEROBJECT_H
#define GEOMETRICCUTUSEROBJECT_H

// MOOSE includes
#include "CrackFrontPointsProvider.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/elem.h"

using namespace libMesh;

class XFEM;

namespace Xfem
{
struct CutEdge
{
  unsigned int id1;
  unsigned int id2;
  Real distance;
  unsigned int host_side_id;
};

struct CutNode
{
  unsigned int id;
  unsigned int host_id;
};

struct CutFace
{
  unsigned int face_id;
  std::vector<unsigned int> face_edge;
  std::vector<Real> position;
};

struct GeomMarkedElemInfo2D
{
  std::vector<CutEdge> elem_cut_edges;
  std::vector<CutNode> elem_cut_nodes;
  std::vector<CutEdge> frag_cut_edges;
  std::vector<std::vector<Point>> frag_edges;
};
//TODO add 3d version

} // namespace Xfem

// Forward declarations
class GeometricCutUserObject;

template <>
InputParameters validParams<GeometricCutUserObject>();

class GeometricCutUserObject : public CrackFrontPointsProvider
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  GeometricCutUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes,
                                    Real time) const = 0;
  virtual bool
  cutElementByGeometry(const Elem * elem, std::vector<Xfem::CutFace> & cut_faces, Real time) const = 0;

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges,
                                     Real time) const = 0;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces,
                                     Real time) const = 0;

protected:
  MooseSharedPointer<XFEM> _xfem;
  std::map<const Elem *, Xfem::GeomMarkedElemInfo2D> _marked_elems_2d;
};

#endif // GEOMETRICCUTUSEROBJECT_H
