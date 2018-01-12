/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GeometricCutUserObject.h"

// MOOSE includes
#include "MooseError.h"
#include "XFEM.h"
#include "EFAElement2D.h"

template <>
InputParameters
validParams<GeometricCutUserObject>()
{
  InputParameters params = validParams<CrackFrontPointsProvider>();
  params.addClassDescription("Base UserObject class for XFEM Geometric Cuts");
  return params;
}

GeometricCutUserObject::GeometricCutUserObject(const InputParameters & parameters)
    : CrackFrontPointsProvider(parameters)
{
  FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblemBase in XFEMMaterialStateMarkerBase");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMMaterialStateMarkerBase");
}

void
GeometricCutUserObject::initialize()
{
  _marked_elems_2d.clear();
}

void
GeometricCutUserObject::execute()
{
  if (_current_elem->dim() == 2)
  {
    std::vector<Xfem::CutEdge> elem_cut_edges;
    std::vector<Xfem::CutNode> elem_cut_nodes;
    std::vector<Xfem::CutEdge> frag_cut_edges;
    std::vector<std::vector<Point>> frag_edges;

    EFAElement2D * EFAElem = _xfem->getEFAElem2D(_current_elem);

    // Don't cut again if elem has been already cut twice
    if (!EFAElem->isFinalCut())
    {
      // get fragment edges
      _xfem->getFragmentEdges(_current_elem, EFAElem, frag_edges);

      // mark cut edges for the element and its fragment
      bool cut = cutElementByGeometry(_current_elem, elem_cut_edges, elem_cut_nodes, _t);
      if (EFAElem->numFragments() > 0)
        cut |= cutFragmentByGeometry(frag_edges, frag_cut_edges, _t);

      if (cut)
      {
        if (_marked_elems_2d.find(_current_elem) != _marked_elems_2d.end())
          mooseError("Element already marked for crack growth");
        auto & me = _marked_elems_2d[_current_elem];
        me.elem_cut_edges = elem_cut_edges;
        me.elem_cut_nodes = elem_cut_nodes;
        me.frag_cut_edges = frag_cut_edges;
        me.frag_edges = frag_edges;
      }
    }
  }
  else if (_current_elem->dim() == 3)
  {
    std::vector<Xfem::CutFace> elem_cut_faces;
    std::vector<Xfem::CutFace> frag_cut_faces;
    std::vector<std::vector<Point>> frag_faces;

    EFAElement3D * EFAElem = _xfem->getEFAElem3D(_current_elem);

    // Don't cut again if elem has been already cut twice
    if (!EFAElem->isFinalCut())
    {
      // get fragment edges
      _xfem->getFragmentFaces(_current_elem, EFAElem, frag_faces);

      // mark cut faces for the element and its fragment
      bool cut = cutElementByGeometry(_current_elem, elem_cut_faces, _t);
      // TODO: This would be done for branching, which is not yet supported in 3D
      //if (EFAElem->numFragments() > 0)
      //  cut |= cutFragmentByGeometry(frag_faces, frag_cut_faces, _t);

      if (cut)
      {
        if (_marked_elems_3d.find(_current_elem) != _marked_elems_3d.end())
          mooseError("Element already marked for crack growth");
        auto & me = _marked_elems_3d[_current_elem];
        me.elem_cut_faces = elem_cut_faces;
        me.frag_cut_faces = frag_cut_faces;
        me.frag_faces = frag_faces;
      }
    }
  }
}

void
GeometricCutUserObject::threadJoin(const UserObject & y)
{
  const GeometricCutUserObject & gcuo = dynamic_cast<const GeometricCutUserObject &>(y);

  for (const auto & it : gcuo._marked_elems_2d)
  {
    mooseAssert(_marked_elems_2d.find(it.first) == _marked_elems_2d.end(),
                "Element already marked for crack growth");
    _marked_elems_2d[it.first] = it.second;
  }
  for (const auto & it : gcuo._marked_elems_3d)
  {
    mooseAssert(_marked_elems_3d.find(it.first) == _marked_elems_3d.end(),
                "Element already marked for crack growth");
    _marked_elems_3d[it.first] = it.second;
  }
}

void
GeometricCutUserObject::finalize()
{

  // Long way of doing this:
  // _communicator.set_union(_marked_elems_2d);
  std::map<const Elem *, Xfem::CutEdge> elem_cut_edges_map;
  std::map<const Elem *, Xfem::CutNode> elem_cut_nodes_map;
  std::map<const Elem *, Xfem::CutEdge> frag_cut_edges_map;
  std::map<std::vector<std::vector<Point>>> frag_edges_map;

  for (auto & it : _marked_elems_2d)
  {
    elem_cut_edges_map[it.first] = it.second.elem_cut_edges;
    elem_cut_nodes_map[it.first] = it.second.elem_cut_nodes;
    frag_cut_edges_map[it.first] = it.second.frag_cut_edges;
    frag_edges_map[it.first] = it.second.frag_edges;
  }

  _communicator.set_union(elem_cut_edges_map);
  _communicator.set_union(elem_cut_nodes_map);
  _communicator.set_union(frag_cut_edges_map);
  _communicator.set_union(frag_edges_map);

  for (auto & it : elem_cut_edges_map)
  {
    auto & me = _marked_elems_2d.find(it.first);
    if (me == _marked_elems_2d.end())
    {
      mooseAssert(elem_cut_nodes_map.find(it.first) != elem_cut_nodes_map.end(),
                  "Map should have entry for this elem");
      mooseAssert(frag_cut_edges_map.find(it.first) != frag_cut_edges_map.end(),
                  "Map should have entry for this elem");
      mooseAssert(frag_edges_map.find(it.first) != frag_edges_map.end(),
                  "Map should have entry for this elem");

      me = _marked_elems_2d[it.first];
      me.elem_cut_edges = elem_cut_edges_map[it.first].elem_cut_edges;
      me.elem_cut_nodes = elem_cut_nodes_map[it.first].elem_cut_nodes;
      me.frag_cut_edges = frag_cut_edges_map[it.first].frag_cut_edges;
      me.frag_edges = frag_edges_map[it.first].frag_edges;
    }
  }

  // Long way of doing this:
  // _communicator.set_union(_marked_elems_3d);
        me.elem_cut_faces = elem_cut_faces;
        me.frag_cut_faces = frag_cut_faces;
        me.frag_faces = frag_faces;

  std::map<const Elem *, Xfem::CutEdge> elem_cut_faces_map;
  std::map<const Elem *, Xfem::CutEdge> frag_cut_faces_map;
  std::map<std::vector<std::vector<Point>>> frag_faces_map;

  for (auto & it : _marked_elems_3d)
  {
    elem_cut_faces_map[it.first] = it.second.elem_cut_faces;
    frag_cut_faces_map[it.first] = it.second.frag_cut_faces;
    frag_faces_map[it.first] = it.second.frag_faces;
  }

  _communicator.set_union(elem_cut_faces_map);
  _communicator.set_union(frag_cut_faces_map);
  _communicator.set_union(frag_faces_map);

  for (auto & it : elem_cut_faces_map)
  {
    auto & me = _marked_elems_3d.find(it.first);
    if (me == _marked_elems_3d.end())
    {
      mooseAssert(frag_cut_faces_map.find(it.first) != frag_cut_faces_map.end(),
                  "Map should have entry for this elem");
      mooseAssert(frag_faces_map.find(it.first) != frag_faces_map.end(),
                  "Map should have entry for this elem");

      me = _marked_elems_3d[it.first];
      me.elem_cut_faces = elem_cut_faces_map[it.first].elem_cut_faces;
      me.frag_cut_faces = frag_cut_faces_map[it.first].frag_cut_faces;
      me.frag_faces = frag_faces_map[it.first].frag_faces;
    }
  }


  //BWS TODO this works only if there is a single geometric cut object
  // Actually -- I don't think we need this because XFEM does it.
  //_xfem->clearGeomMarkedElems();

  for (const auto & it : _marked_elems_2d)
    _xfem->addGeomMarkedElem2D(it.first, it.second);

  for (const auto & it : _marked_elems_3d)
    _xfem->addGeomMarkedElem3D(it.first, it.second);

  _marked_elems_2d.clear();
  _marked_elems_3d.clear();
}
