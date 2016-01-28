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

#ifndef EFAELEMENT2D_H
#define EFAELEMENT2D_H

#include "EFAElement.h"
#include "EFAFace.h"
#include "EFAFragment2D.h"

class EFAElement2D : public EFAElement
{
public:

  EFAElement2D(unsigned int eid, unsigned int n_nodes);
  EFAElement2D(const EFAElement2D * from_elem, bool convert_to_local);
  EFAElement2D(const EFAFace * from_face);

  ~EFAElement2D();

private:

  unsigned int _num_edges;
  std::vector<EFAEdge*> _edges;
  std::vector<FaceNode*> _interior_nodes;
  std::vector<std::vector<EFAElement2D*> >_edge_neighbors;
  std::vector<EFAFragment2D*> _fragments;

public:

  // override virtual methods in base class
  virtual unsigned int num_frags() const;
  virtual bool is_partial() const;
  virtual void get_non_physical_nodes(std::set<EFANode*> &non_physical_nodes) const;

  virtual void switchNode(EFANode *new_node, EFANode *old_node, bool descend_to_parent);
  virtual void switchEmbeddedNode(EFANode *new_node, EFANode *old_node);
  virtual void getMasterInfo(EFANode* node, std::vector<EFANode*> &master_nodes,
                             std::vector<double> &master_weights) const;
  virtual unsigned int num_interior_nodes() const;

  virtual bool overlays_elem(const EFAElement* other_elem) const;
  virtual unsigned int get_neighbor_index(const EFAElement* neighbor_elem) const;
  virtual void clear_neighbors();
  virtual void setup_neighbors(std::map<EFANode*, std::set<EFAElement*> > &InverseConnectivityMap);
  virtual void neighbor_sanity_check() const;

  virtual void init_crack_tip(std::set<EFAElement*> &CrackTipElements);
  virtual bool should_duplicate_for_crack_tip(const std::set<EFAElement*> &CrackTipElements);
  virtual bool shouldDuplicateCrackTipSplitElem(const std::set<EFAElement*> &CrackTipElements);
  virtual bool shouldDuplicateForPhantomCorner();
  virtual bool will_crack_tip_extend(std::vector<unsigned int> &split_neighbors) const;
  virtual bool is_crack_tip_elem() const;

  virtual unsigned int get_num_cuts() const;
  virtual bool is_final_cut() const;
  virtual void update_fragments(const std::set<EFAElement*> &CrackTipElements,
                                std::map<unsigned int, EFANode*> &EmbeddedNodes);
  virtual void fragment_sanity_check(unsigned int n_old_frag_edges, unsigned int n_old_frag_cuts) const;
  virtual void restore_fragment(const EFAElement* const from_elem);

  virtual void create_child(const std::set<EFAElement*> &CrackTipElements,
                            std::map<unsigned int, EFAElement*> &Elements,
                            std::map<unsigned int, EFAElement*> &newChildElements,
                            std::vector<EFAElement*> &ChildElements,
                            std::vector<EFAElement*> &ParentElements,
                            std::map<unsigned int, EFANode*> &TempNodes);
  virtual void remove_phantom_embedded_nodes();
  virtual void connect_neighbors(std::map<unsigned int, EFANode*> &PermanentNodes,
                                 std::map<unsigned int, EFANode*> &TempNodes,
                                 std::map<EFANode*, std::set<EFAElement*> > &InverseConnectivityMap,
                                 bool merge_phantom_edges);
  virtual void print_elem();

  // EFAelement2D specific methods
  EFAFragment2D* get_fragment(unsigned int frag_id) const;
  std::set<EFANode*> get_edge_nodes(unsigned int edge_id) const;
  bool getEdgeNodeParaCoor(EFANode* node, std::vector<double> &para_coor) const;
  FaceNode* get_interior_node(unsigned int interior_node_id) const;
  void delete_interior_nodes();

  unsigned int num_edges() const;
  void set_edge(unsigned int edge_id, EFAEdge* edge);
  void createEdges();
  EFAEdge* get_edge(unsigned int edge_id) const;

  EFAEdge* get_frag_edge(unsigned int frag_id, unsigned int edge_id) const;
  std::set<EFANode*> getPhantomNodeOnEdge(unsigned int edge_id) const;
  bool getFragmentEdgeID(unsigned int elem_edge_id, unsigned int &frag_edge_id) const;
  bool is_edge_phantom(unsigned int edge_id) const;

  unsigned int num_edge_neighbors(unsigned int edge_id) const;
  EFAElement2D* get_edge_neighbor(unsigned int edge_id, unsigned int neighbor_id) const;

  unsigned int get_crack_tip_split_element_id() const;

  bool frag_has_tip_edges() const;
  unsigned int get_tip_edge_id() const;
  EFANode* get_tip_embedded() const;
  bool edge_contains_tip(unsigned int edge_id) const;
  bool frag_edge_already_cut(unsigned int ElemEdgeID) const;

  void add_edge_cut(unsigned int edge_id, double position, EFANode* embedded_node,
                    std::map< unsigned int, EFANode*> &EmbeddedNodes, bool add_to_neighbor);
  void add_frag_edge_cut(unsigned int frag_edge_id, double position,
                         std::map< unsigned int, EFANode*> &EmbeddedNodes);
  std::vector<EFAFragment2D*> branching_split(std::map<unsigned int, EFANode*> &EmbeddedNodes);

private:

  // EFAelement2D unique methods
  void mapParaCoorFrom1Dto2D(unsigned int edge_id, double xi_1d,
                             std::vector<double> &para_coor) const;
};

#endif
