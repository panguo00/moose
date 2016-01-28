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

#ifndef EFAFACE_H
#define EFAFACE_H

#include "EFAEdge.h"
#include "EFAFragment2D.h"
#include "FaceNode.h"

class EFAFragment2D;

class EFAFace
{
public:

  EFAFace(unsigned int n_nodes);
  EFAFace(const EFAFace & other_face);
  EFAFace(const EFAFragment2D * frag);

  ~EFAFace();

private:

  unsigned int _num_nodes;
  std::vector<EFANode*> _nodes;
  unsigned int _num_edges;
  std::vector<EFAEdge*> _edges;
  std::vector<FaceNode*> _interior_nodes;

public:

  unsigned int num_nodes() const;
  void set_node(unsigned int node_id, EFANode* node);
  EFANode* get_node(unsigned int node_id) const;
  void switchNode(EFANode *new_node, EFANode *old_node);
  bool getMasterInfo(EFANode* node, std::vector<EFANode*> &master_nodes,
                     std::vector<double> &master_weights) const;
  bool getEdgeNodeParaCoor(EFANode* node, std::vector<double> &xi_2d) const;
  bool getFaceNodeParaCoor(EFANode* node, std::vector<double> &xi_2d) const;
  unsigned int num_interior_nodes() const;
  void createNodes();

  unsigned int num_edges() const;
  EFAEdge* get_edge(unsigned int edge_id) const;
  void set_edge(unsigned int edge_id, EFAEdge* new_edge);
  void createEdges();
  void combine_two_edges(unsigned int edge_id1, unsigned int edge_id2);
  void sort_edges();
  void reverse_edges();
  bool is_trig_quad() const;

  bool equivalent(const EFAFace* other_face) const;
  bool containsNode(const EFANode* node) const;
  bool containsFace(const EFAFace* other_face) const;
  bool doesOwnEdge(const EFAEdge* other_edge) const;
  void remove_embedded_node(EFANode* emb_node);
  std::vector<EFAFace*> split() const;
  EFAFace* combine_with(const EFAFace* other_face) const;
  void reset_edge_intersection(const EFAFace* ref_face);

  unsigned int get_num_cuts() const;
  bool has_intersection() const;
  void copy_intersection(const EFAFace &from_face);
  bool isAdjacent(const EFAFace* other_face) const;
  unsigned int adjacentCommonEdge(const EFAFace* other_face) const;
  bool is_same_orientation(const EFAFace* other_face) const;
  FaceNode* get_interior_node(unsigned int index) const;

private:
  void mapParaCoorFrom1Dto2D(unsigned int edge_id, double xi_1d,
                             std::vector<double> &xi_2d) const;
};

#endif
