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
        mooseAssert(_marked_elems_2d.find(_current_elem) != _marked_elems_2d.end(),
                    "Element already marked for crack growth");
        auto & me = _marked_elems_2d[_current_elem];
        me.elem_cut_edges = elem_cut_edges;
        me.elem_cut_nodes = elem_cut_nodes;
        me.frag_cut_edges = frag_cut_edges;
        me.frag_edges = frag_edges;
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
    mooseAssert(_marked_elems_2d.find(it.first) != _marked_elems_2d.end(),
                "Element already marked for crack growth");
    _marked_elems_2d[it.first] = it.second;
  }
}

void
GeometricCutUserObject::finalize()
{
  std::map<const Elem *, Xfem::CutEdge> elem_cut_edges_map;
  for (auto & it : _marked_elems_2d)
  {
    elem_cut_edges_map[it.first] = it.second.elem_cut_edges;
  }



//  _communicator.set_union(_marked_elems_2d);

  //BWS TODO this works only if there is a single geometric cut object
  // Actually -- I don't think we need this because XFEM does it.
  //_xfem->clearGeomMarkedElems();

  for (const auto & it : _marked_elems_2d)
    _xfem->addGeomMarkedElem2D(it.first, it.second);

  _marked_elems_2d.clear();
}
