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

#ifndef EFAFRAGMENT_H
#define EFAFRAGMENT_H

#include "EFAEdge.h"

class EFAFragment
{
public:

  EFAFragment();
  virtual ~EFAFragment();

  virtual void switchNode(EFANode *new_node, EFANode *old_node) = 0;
  virtual bool containsNode(EFANode *node) const = 0;
  virtual unsigned int get_num_cuts() const = 0;
  virtual std::set<EFANode*> get_all_nodes() const = 0;
  virtual bool isConnected(EFAFragment *other_fragment) const = 0;
  virtual void remove_invalid_embedded(std::map<unsigned int, EFANode*> &EmbeddedNodes) = 0;

  // common methods
  std::vector<EFANode*> get_common_nodes(EFAFragment* other) const;
};

#endif
