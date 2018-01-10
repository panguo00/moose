/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GeometricCutUserObject.h"

// MOOSE includes
#include "MooseError.h"

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
  _marked_elems.clear();
}

void
GeometricCutUserObject::execute()
{
  if (_current_elem->dim() == 2)
  {
    std::vector<CutEdge> elem_cut_edges;
    std::vector<CutNode> elem_cut_nodes;
    std::vector<CutEdge> frag_cut_edges;
    std::vector<std::vector<Point>> frag_edges;

    EFAElement2D * EFAElem = getEFAElem2D(elem);

    // Don't cut again if elem has been already cut twice
    if (EFAElem->isFinalCut())
      continue;

    // get fragment edges
    getFragmentEdges(elem, EFAElem, frag_edges);

    // mark cut edges for the element and its fragment
    for (unsigned int i = 0; i < active_geometric_cuts.size(); ++i)
    {
      active_geometric_cuts[i]->cutElementByGeometry(elem, elem_cut_edges, elem_cut_nodes, time);
      if (EFAElem->numFragments() > 0)
	active_geometric_cuts[i]->cutFragmentByGeometry(frag_edges, frag_cut_edges, time);
    }
  }
}

void
GeometricCutUserObject::threadJoin(const UserObject & y)
{
  const GeometricCutUserObject & gcuo = dynamic_cast<const GeometricCutUserObject &>(y);

  for (std::map<unsigned int, RealVectorValue>::const_iterator mit = gcuo._marked_elems.begin();
       mit != gcuo._marked_elems.end();
       ++mit)
  {
    _marked_elems[mit->first] = mit->second; // TODO do error checking for duplicates here too
  }

  for (std::set<unsigned int>::const_iterator mit = gcuo._marked_frags.begin();
       mit != gcuo._marked_frags.end();
       ++mit)
  {
    _marked_frags.insert(*mit); // TODO do error checking for duplicates here too
  }

  for (std::map<unsigned int, unsigned int>::const_iterator mit = gcuo._marked_elem_sides.begin();
       mit != gcuo._marked_elem_sides.end();
       ++mit)
  {
    _marked_elem_sides[mit->first] = mit->second; // TODO do error checking for duplicates here too
  }
}

void
GeometricCutUserObject::finalize()
{
  _communicator.set_union(_marked_elems);
  _communicator.set_union(_marked_frags);
  _communicator.set_union(_marked_elem_sides);

  _xfem->clearStateMarkedElems();
  std::map<unsigned int, RealVectorValue>::iterator mit;
  for (mit = _marked_elems.begin(); mit != _marked_elems.end(); ++mit)
  {
    if (_marked_elem_sides.find(mit->first) != _marked_elem_sides.end())
      _xfem->addStateMarkedElem(mit->first, mit->second, _marked_elem_sides[mit->first]);
    else if (_marked_frags.find(mit->first) != _marked_frags.end())
      _xfem->addStateMarkedFrag(mit->first, mit->second);
    else
      _xfem->addStateMarkedElem(mit->first, mit->second);
  }

  _marked_elems.clear();
  _marked_frags.clear();
  _marked_elem_sides.clear();
}
