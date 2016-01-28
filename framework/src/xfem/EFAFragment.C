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

#include "EFAFragment.h"

#include "EFAElement.h"
#include "EFAfuncs.h"

EFAFragment::EFAFragment()
{}

EFAFragment::~EFAFragment()
{}

std::vector<EFANode*>
EFAFragment::get_common_nodes(EFAFragment* other) const
{
  std::set<EFANode*> frag1_nodes = get_all_nodes();
  std::set<EFANode*> frag2_nodes = other->get_all_nodes();
  std::vector<EFANode*> common_nodes = get_common_elems(frag1_nodes, frag2_nodes);
  return common_nodes;
}
