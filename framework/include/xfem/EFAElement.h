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

#ifndef EFAELEMENT_H
#define EFAELEMENT_H

#include "EFAEdge.h"
#include "EFAFragment.h"
#include "FaceNode.h"

class EFAElement
{
public:

  EFAElement(unsigned int eid, unsigned int n_nodes);

  virtual ~EFAElement();

protected:

  unsigned int _id;
  unsigned int _num_nodes;
  std::vector<EFANode*> _nodes;
  std::vector<EFANode*> _local_nodes;
  EFAElement* _parent;
  std::vector<EFAElement*> _children;
  bool _crack_tip_split_element;
  std::vector<unsigned int> _crack_tip_neighbors;
  std::vector<EFAElement*> _general_neighbors; // all elements sharing at least one node with curr elem

public:

  // common methods
  unsigned int id() const;
  unsigned int num_nodes() const;
  void set_node(unsigned int node_id, EFANode* node);
  EFANode* get_node(unsigned int node_id) const;
  bool containsNode(EFANode* node) const;
  void display_nodes() const;
  EFANode * create_local_node_from_global_node(const EFANode * global_node) const;
  EFANode * get_global_node_from_local_node(const EFANode * local_node) const;
  unsigned int getLocalNodeIndex(EFANode * node) const;
  std::vector<EFANode*> get_common_nodes(const EFAElement* other_elem) const;

  void set_crack_tip_split();
  bool is_crack_tip_split() const;
  unsigned int num_crack_tip_neighbors() const;
  unsigned int get_crack_tip_neighbor(unsigned int index) const;
  void add_crack_tip_neighbor(EFAElement * neighbor_elem);

  EFAElement* parent() const;
  EFAElement* get_child(unsigned int child_id) const;
  void set_parent(EFAElement* parent);
  unsigned int num_children() const;
  void add_child(EFAElement* child);
  void remove_parent_children();
  std::vector<EFAElement*> get_general_neighbors(std::map<EFANode*, std::set<EFAElement*> > &InverseConnectivity) const;
  EFAElement* get_general_neighbor(unsigned int index) const;
  unsigned int num_general_neighbors() const;

  // pure virtual methods
  virtual unsigned int num_frags() const = 0;
  virtual bool is_partial() const = 0;
  virtual void get_non_physical_nodes(std::set<EFANode*> &non_physical_nodes) const = 0;

  virtual void switchNode(EFANode *new_node, EFANode *old_node, bool descend_to_parent) = 0;
  virtual void switchEmbeddedNode(EFANode *new_node, EFANode *old_node) = 0;
  virtual void getMasterInfo(EFANode* node, std::vector<EFANode*> &master_nodes,
                             std::vector<double> &master_weights) const = 0;
  virtual unsigned int num_interior_nodes() const = 0;

  virtual bool overlays_elem(const EFAElement* other_elem) const = 0;
  virtual unsigned int get_neighbor_index(const EFAElement * neighbor_elem) const = 0;
  virtual void clear_neighbors() = 0;
  virtual void setup_neighbors(std::map<EFANode*, std::set<EFAElement*> > &InverseConnectivityMap) = 0;
  virtual void neighbor_sanity_check() const = 0;

  virtual void init_crack_tip(std::set< EFAElement*> &CrackTipElements) = 0;
  virtual bool should_duplicate_for_crack_tip(const std::set<EFAElement*> &CrackTipElements) = 0;
  virtual bool shouldDuplicateCrackTipSplitElem(const std::set<EFAElement*> &CrackTipElements) = 0;
  virtual bool shouldDuplicateForPhantomCorner() = 0;
  virtual bool will_crack_tip_extend(std::vector<unsigned int> &split_neighbors) const = 0;
  virtual bool is_crack_tip_elem() const = 0;

  virtual unsigned int get_num_cuts() const = 0;
  virtual bool is_final_cut() const = 0;
  virtual void update_fragments(const std::set<EFAElement*> &CrackTipElements,
                                std::map<unsigned int, EFANode*> &EmbeddedNodes) = 0;
  virtual void fragment_sanity_check(unsigned int n_old_frag_edges,
                                     unsigned int n_old_frag_cuts) const = 0;
  virtual void restore_fragment(const EFAElement* const from_elem) = 0;

  virtual void create_child(const std::set<EFAElement*> &CrackTipElements,
                            std::map<unsigned int, EFAElement*> &Elements,
                            std::map<unsigned int, EFAElement*> &newChildElements,
                            std::vector<EFAElement*> &ChildElements,
                            std::vector<EFAElement*> &ParentElements,
                            std::map<unsigned int, EFANode*> &TempNodes) = 0;
  virtual void remove_phantom_embedded_nodes() = 0;
  virtual void connect_neighbors(std::map<unsigned int, EFANode*> &PermanentNodes,
                                 std::map<unsigned int, EFANode*> &TempNodes,
                                 std::map<EFANode*, std::set<EFAElement*> > &InverseConnectivityMap,
                                 bool merge_phantom_edges) = 0;
  virtual void print_elem() = 0;

protected:

  // common methods
  void mergeNodes(EFANode* &childNode, EFANode* &childOfNeighborNode, EFAElement* childOfNeighborElem,
                  std::map<unsigned int, EFANode*> &PermanentNodes, std::map<unsigned int, EFANode*> &TempNodes);
};

#endif
