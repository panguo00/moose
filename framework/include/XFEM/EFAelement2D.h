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

#include "EFAfragment2D.h"
#include "EFAelement.h"
#include "EFAface.h"

class EFAelement2D : public EFAelement
{
public:

  EFAelement2D(unsigned int eid, unsigned int n_nodes);
  EFAelement2D(const EFAelement2D * from_elem, bool convert_to_local);
  EFAelement2D(const EFAface * from_face);

  ~EFAelement2D();

private:

  unsigned int _num_edges;
  std::vector<EFAedge*> _edges;
  std::vector<FaceNode*> _interior_nodes;
  std::vector<std::vector<EFAelement2D*> >_edge_neighbors;
  std::vector<EFAfragment2D*> _fragments;

public:

  // override virtual methods in base class
  unsigned int num_frags() const;
  bool is_partial() const;
  void get_non_physical_nodes(std::set<EFAnode*> &non_physical_nodes) const;

  void switchNode(EFAnode *new_node, EFAnode *old_node, bool descend_to_parent);
  void switchEmbeddedNode(EFAnode *new_node, EFAnode *old_node);
  void getMasterInfo(EFAnode* node, std::vector<EFAnode*> &master_nodes,
                     std::vector<double> &master_weights) const;
  unsigned int num_interior_nodes() const;

  bool overlays_elem(const EFAelement* other_elem) const;
  unsigned int get_neighbor_index(const EFAelement* neighbor_elem) const;
  void clear_neighbors();
  void setup_neighbors(std::map<EFAnode*, std::set<EFAelement*> > &InverseConnectivityMap);
  void neighbor_sanity_check() const;

  void init_crack_tip(std::set<EFAelement*> &CrackTipElements);
  bool should_duplicate_for_crack_tip(const std::set<EFAelement*> &CrackTipElements);
  bool shouldDuplicateCrackTipSplitElem(const std::set<EFAelement*> &CrackTipElements);
  bool shouldDuplicateForPhantomCorner();
  bool will_crack_tip_extend(std::vector<unsigned int> &split_neighbors) const;
  bool is_crack_tip_elem() const;

  unsigned int get_num_cuts() const;
  bool is_cut_twice() const;
  void update_fragments(const std::set<EFAelement*> &CrackTipElements,
                        std::map<unsigned int, EFAnode*> &EmbeddedNodes);
  void fragment_sanity_check(unsigned int n_old_frag_edges, unsigned int n_old_frag_cuts) const;
  void restore_fragment(const EFAelement* const from_elem);

  void create_child(const std::set<EFAelement*> &CrackTipElements,
                    std::map<unsigned int, EFAelement*> &Elements,
                    std::map<unsigned int, EFAelement*> &newChildElements,
                    std::vector<EFAelement*> &ChildElements,
                    std::vector<EFAelement*> &ParentElements,
                    std::map<unsigned int, EFAnode*> &TempNodes);
  void remove_phantom_embedded_nodes();
  void connect_neighbors(std::map<unsigned int, EFAnode*> &PermanentNodes,
                         std::map<unsigned int, EFAnode*> &EmbeddedNodes,
                         std::map<unsigned int, EFAnode*> &TempNodes,
                         std::map<EFAnode*, std::set<EFAelement*> > &InverseConnectivityMap,
                         bool merge_phantom_edges);
  void print_elem();

  // EFAelement2D specific methods
  EFAfragment2D* get_fragment(unsigned int frag_id) const;
  std::set<EFAnode*> get_edge_nodes(unsigned int edge_id) const;
  bool getEdgeNodeParaCoor(EFAnode* node, std::vector<double> &para_coor) const;
  FaceNode* get_interior_node(unsigned int interior_node_id) const;

  unsigned int num_edges() const;
  void set_edge(unsigned int edge_id, EFAedge* edge);
  void createEdges();
  EFAedge* get_edge(unsigned int edge_id) const;

  EFAedge* get_frag_edge(unsigned int frag_id, unsigned int edge_id) const;
  std::set<EFAnode*> getPhantomNodeOnEdge(unsigned int edge_id) const;
  bool getFragmentEdgeID(unsigned int elem_edge_id, unsigned int &frag_edge_id) const;
  bool is_edge_phantom(unsigned int edge_id) const;

  unsigned int num_edge_neighbors(unsigned int edge_id) const;
  EFAelement2D* get_edge_neighbor(unsigned int edge_id, unsigned int neighbor_id) const;

  bool frag_has_tip_edges() const;
  unsigned int get_tip_edge_id() const;
  EFAnode* get_tip_embedded() const;
  bool edge_contains_tip(unsigned int edge_id) const;

  void add_edge_cut(unsigned int edge_id, double position, EFAnode* embedded_node,
                    std::map< unsigned int, EFAnode*> &EmbeddedNodes, bool add_to_neighbor);
  void add_frag_edge_cut(unsigned int frag_edge_id, double position,
                         std::map< unsigned int, EFAnode*> &EmbeddedNodes);

private:

  // EFAelement2D unique methods
  void mapParaCoorFrom1Dto2D(unsigned int edge_id, double xi_1d,
                             std::vector<double> &para_coor) const;
};

#endif