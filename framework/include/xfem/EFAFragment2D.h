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

#ifndef EFAFRAGMENT2D_H
#define EFAFRAGMENT2D_H

#include "EFAEdge.h"
#include "EFAFace.h"
#include "EFAFragment.h"

class EFAElement2D;
class EFAFace;

class EFAFragment2D : public EFAFragment
{
public:

  EFAFragment2D(EFAElement2D * host, bool create_boundary_edges,
                const EFAElement2D * from_host,
                unsigned int frag_id = std::numeric_limits<unsigned int>::max());
  EFAFragment2D(EFAElement2D* host, const EFAFace * from_face);
  ~EFAFragment2D();

private:

  EFAElement2D * _host_elem;
  std::vector<EFAEdge*> _boundary_edges;

public:
  // override pure virtual methods
  virtual void switchNode(EFANode *new_node, EFANode *old_node);
  virtual bool containsNode(EFANode *node) const;
  virtual unsigned int get_num_cuts() const;
  virtual std::set<EFANode*> get_all_nodes() const;
  virtual bool isConnected(EFAFragment *other_fragment) const;
  virtual void remove_invalid_embedded(std::map<unsigned int, EFANode*> &EmbeddedNodes);

  // EFAfragment2D specific methods
  void combine_tip_edges();
  bool is_edge_interior(unsigned int edge_id) const;
  std::vector<unsigned int> get_interior_edge_id() const;
  bool isSecondaryInteriorEdge(unsigned int edge_id) const;
  unsigned int num_edges() const;
  EFAEdge* get_edge(unsigned int edge_id) const;
  void add_edge(EFAEdge* new_edge);
  std::set<EFANode*> get_edge_nodes(unsigned int edge_id) const;
  EFAElement2D * get_host() const;
  std::vector<EFAFragment2D*> split();
};

#endif
