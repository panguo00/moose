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

#ifndef EFAEDGE_H
#define EFAEDGE_H

#include "FaceNode.h"

class EFAEdge
{
  public:

  EFAEdge(EFANode * node1, EFANode * node2);
  EFAEdge(const EFAEdge & other_edge);

  ~EFAEdge();

  private:

  EFANode * _edge_node1;
  EFANode * _edge_node2;
  std::vector<EFANode*> _embedded_nodes;
  std::vector<double> _intersection_x;

  public:

  bool equivalent(const EFAEdge & other) const;
  bool isPartialOverlap(const EFAEdge & other) const;
  bool containsEdge(const EFAEdge & other) const;
  bool getNodeMasters(EFANode* node, std::vector<EFANode*> &master_nodes,
                      std::vector<double> &master_weights) const;
//  bool operator < (const EFAEdge & other) const;

  void add_intersection(double position, EFANode * embedded_node_tmp, EFANode * from_node);
  void reset_intersection(double position, EFANode * embedded_node_tmp, EFANode * from_node);
  void copy_intersection(const EFAEdge & other, unsigned int from_node_id);
  EFANode * get_node(unsigned int index) const;
  void reverse_nodes();

  bool has_intersection() const;
  bool has_intersection_at_position(double position, EFANode * from_node) const;
  double get_intersection(unsigned int emb_id, EFANode * from_node) const;
  double distance_from_node1(EFANode * node) const;
  bool is_embedded_node(const EFANode * node) const;
  unsigned int get_embedded_index(EFANode * node) const;
  unsigned int get_embedded_index(double position, EFANode* from_node) const;

  EFANode * get_embedded_node(unsigned int index) const;
  unsigned int num_embedded_nodes() const;
  void consistency_check();
  void switchNode(EFANode *new_node, EFANode *old_node);
  bool containsNode(const EFANode *node) const;
  void remove_embedded_node();
  void remove_embedded_node(EFANode * node);
};

#endif
