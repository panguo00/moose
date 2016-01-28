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

#include "EfaEdge.h"
#include "FaceNode.h"
#include "EFAfragment.h"

class EfaElement
{
public:

  EfaElement(unsigned int eid, unsigned int n_nodes);

  virtual ~EfaElement();

protected:

  unsigned int _id;
  unsigned int _num_nodes;
  std::vector<EFAnode*> _nodes;
  std::vector<EFAnode*> _local_nodes;
  EfaElement* _parent;
  std::vector<EfaElement*> _children;
  bool _crack_tip_split_element;
  std::vector<unsigned int> _crack_tip_neighbors;
  std::vector<EfaElement*> _general_neighbors; // all elements sharing at least one node with curr elem

public:

  // common methods
  unsigned int id() const;
  unsigned int num_nodes() const;
  void set_node(unsigned int node_id, EFAnode* node);
  EFAnode* get_node(unsigned int node_id) const;
  bool containsNode(EFAnode* node) const;
  void display_nodes() const;
  EFAnode * create_local_node_from_global_node(const EFAnode * global_node) const;
  EFAnode * get_global_node_from_local_node(const EFAnode * local_node) const;
  unsigned int getLocalNodeIndex(EFAnode * node) const;
  std::vector<EFAnode*> get_common_nodes(const EfaElement* other_elem) const;

  void set_crack_tip_split();
  bool is_crack_tip_split() const;
  unsigned int num_crack_tip_neighbors() const;
  unsigned int get_crack_tip_neighbor(unsigned int index) const;
  void add_crack_tip_neighbor(EfaElement * neighbor_elem);

  EfaElement* parent() const;
  EfaElement* get_child(unsigned int child_id) const;
  void set_parent(EfaElement* parent);
  unsigned int num_children() const;
  void add_child(EfaElement* child);
  void remove_parent_children();
  std::vector<EfaElement*> get_general_neighbors(std::map<EFAnode*, std::set<EfaElement*> > &InverseConnectivity) const;
  EfaElement* get_general_neighbor(unsigned int index) const;
  unsigned int num_general_neighbors() const;

  // pure virtual methods
  virtual unsigned int num_frags() const = 0;
  virtual bool is_partial() const = 0;
  virtual void get_non_physical_nodes(std::set<EFAnode*> &non_physical_nodes) const = 0;

  virtual void switchNode(EFAnode *new_node, EFAnode *old_node, bool descend_to_parent) = 0;
  virtual void switchEmbeddedNode(EFAnode *new_node, EFAnode *old_node) = 0;
  virtual void getMasterInfo(EFAnode* node, std::vector<EFAnode*> &master_nodes,
                             std::vector<double> &master_weights) const = 0;
  virtual unsigned int num_interior_nodes() const = 0;

  virtual bool overlays_elem(const EfaElement* other_elem) const = 0;
  virtual unsigned int get_neighbor_index(const EfaElement * neighbor_elem) const = 0;
  virtual void clear_neighbors() = 0;
  virtual void setup_neighbors(std::map<EFAnode*, std::set<EfaElement*> > &InverseConnectivityMap) = 0;
  virtual void neighbor_sanity_check() const = 0;

  virtual void init_crack_tip(std::set< EfaElement*> &CrackTipElements) = 0;
  virtual bool should_duplicate_for_crack_tip(const std::set<EfaElement*> &CrackTipElements) = 0;
  virtual bool shouldDuplicateCrackTipSplitElem(const std::set<EfaElement*> &CrackTipElements) = 0;
  virtual bool shouldDuplicateForPhantomCorner() = 0;
  virtual bool will_crack_tip_extend(std::vector<unsigned int> &split_neighbors) const = 0;
  virtual bool is_crack_tip_elem() const = 0;

  virtual unsigned int get_num_cuts() const = 0;
  virtual bool is_final_cut() const = 0;
  virtual void update_fragments(const std::set<EfaElement*> &CrackTipElements,
                                std::map<unsigned int, EFAnode*> &EmbeddedNodes) = 0;
  virtual void fragment_sanity_check(unsigned int n_old_frag_edges,
                                     unsigned int n_old_frag_cuts) const = 0;
  virtual void restore_fragment(const EfaElement* const from_elem) = 0;

  virtual void create_child(const std::set<EfaElement*> &CrackTipElements,
                            std::map<unsigned int, EfaElement*> &Elements,
                            std::map<unsigned int, EfaElement*> &newChildElements,
                            std::vector<EfaElement*> &ChildElements,
                            std::vector<EfaElement*> &ParentElements,
                            std::map<unsigned int, EFAnode*> &TempNodes) = 0;
  virtual void remove_phantom_embedded_nodes() = 0;
  virtual void connect_neighbors(std::map<unsigned int, EFAnode*> &PermanentNodes,
                                 std::map<unsigned int, EFAnode*> &TempNodes,
                                 std::map<EFAnode*, std::set<EfaElement*> > &InverseConnectivityMap,
                                 bool merge_phantom_edges) = 0;
  virtual void print_elem() = 0;

protected:

  // common methods
  void mergeNodes(EFAnode* &childNode, EFAnode* &childOfNeighborNode, EfaElement* childOfNeighborElem,
                  std::map<unsigned int, EFAnode*> &PermanentNodes, std::map<unsigned int, EFAnode*> &TempNodes);
};

#endif
