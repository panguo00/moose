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

#ifndef EFAELEMENT3D_H
#define EFAELEMENT3D_H

#include "EFAElement.h"
#include "EFAFragment3D.h"
#include "VolumeNode.h"

class EFAElement3D : public EFAElement
{
public:

  EFAElement3D(unsigned int eid, unsigned int n_nodes, unsigned int n_faces);
  EFAElement3D(const EFAElement3D * from_elem, bool convert_to_local);

  ~EFAElement3D();

private:

  unsigned int _num_faces;
  std::vector<EFAFace*> _faces;
  std::vector<VolumeNode*> _interior_nodes;
  std::vector<std::vector<EFAElement3D*> > _face_neighbors;
  std::vector<EFAFragment3D*> _fragments;
  std::vector<std::vector<EFAFace*> > _adjacent_face_ix;

public:

  // override virtual methods in base class
  virtual unsigned int numFragments() const;
  virtual bool isPartial() const;
  virtual void getNonPhysicalNodes(std::set<EFANode*> &non_physical_nodes) const;

  virtual void switchNode(EFANode *new_node, EFANode *old_node, bool descend_to_parent);
  virtual void switchEmbeddedNode(EFANode *new_node, EFANode *old_node);
  virtual void getMasterInfo(EFANode* node, std::vector<EFANode*> &master_nodes,
                             std::vector<double> &master_weights) const;
  virtual unsigned int numInteriorNodes() const;

  virtual bool overlaysElement(const EFAElement* other_elem) const;
  virtual unsigned int getNeighborIndex(const EFAElement* neighbor_elem) const;
  virtual void clearNeighbors();
  virtual void setupNeighbors(std::map<EFANode*, std::set<EFAElement*> > &InverseConnectivityMap);
  virtual void neighborSanityCheck() const;

  virtual void initCrackTip(std::set<EFAElement*> &CrackTipElements);
  virtual bool shouldDuplicateForCrackTip(const std::set<EFAElement*> &CrackTipElements);
  virtual bool shouldDuplicateCrackTipSplitElement(const std::set<EFAElement*> &CrackTipElements);
  virtual bool shouldDuplicateForPhantomCorner();
  virtual bool willCrackTipExtend(std::vector<unsigned int> &split_neighbors) const;
  virtual bool isCrackTipElement() const;

  virtual unsigned int getNumCuts() const;
  virtual bool isFinalCut() const;
  virtual void updateFragments(const std::set<EFAElement*> &CrackTipElements,
                                std::map<unsigned int, EFANode*> &EmbeddedNodes);
  virtual void fragmentSanityCheck(unsigned int n_old_frag_faces, unsigned int n_old_frag_cuts) const;
  virtual void restoreFragment(const EFAElement* const from_elem);

  virtual void createChild(const std::set<EFAElement*> &CrackTipElements,
                            std::map<unsigned int, EFAElement*> &Elements,
                            std::map<unsigned int, EFAElement*> &newChildElements,
                            std::vector<EFAElement*> &ChildElements,
                            std::vector<EFAElement*> &ParentElements,
                            std::map<unsigned int, EFANode*> &TempNodes);
  virtual void removePhantomEmbeddedNode();
  virtual void connectNeighbors(std::map<unsigned int, EFANode*> &PermanentNodes,
                                 std::map<unsigned int, EFANode*> &TempNodes,
                                 std::map<EFANode*, std::set<EFAElement*> > &InverseConnectivityMap,
                                 bool merge_phantom_faces);
  virtual void printElement();

  // EFAelement3D specific methods
  EFAFragment3D* get_fragment(unsigned int frag_id) const;
  std::set<EFANode*> get_face_nodes(unsigned int face_id) const;
  bool getFaceNodeParaCoor(EFANode* node, std::vector<double> &xi_3d) const;
  VolumeNode* get_interior_node(unsigned int interior_node_id) const;
  void remove_embedded_node(EFANode* emb_node, bool remove_for_neighbor);

  unsigned int num_faces() const;
  void set_face(unsigned int face_id, EFAFace* face);
  void createFaces();
  EFAFace* get_face(unsigned int face_id) const;
  unsigned int get_face_id(EFAFace* face) const;
  std::vector<unsigned int> get_common_face_id(const EFAElement3D* other_elem) const;
  unsigned int getNeighborFaceNodeID(unsigned int face_id, unsigned int node_id,
                                     EFAElement3D* neighbor_elem) const;
  unsigned int getNeighborFaceEdgeID(unsigned int face_id, unsigned int edg_id,
                                     EFAElement3D* neighbor_elem) const;
  void create_adjacent_face_ix();
  EFAFace* get_adjacent_face(unsigned int face_id, unsigned int edge_id) const;

  EFAFace* get_frag_face(unsigned int frag_id, unsigned int face_id) const;
  std::set<EFANode*> getPhantomNodeOnFace(unsigned int face_id) const;
  bool getFragmentFaceID(unsigned int elem_face_id, unsigned int &frag_face_id) const;
  bool getFragmentFaceEdgeID(unsigned int ElemFaceID, unsigned int ElemFaceEdgeID,
                             unsigned int &FragFaceID, unsigned int &FragFaceEdgeID) const;
  bool is_real_edge_cut(unsigned int ElemFaceID, unsigned int ElemFaceEdgeID, double position) const;
  bool is_face_phantom(unsigned int face_id) const;
  unsigned int num_face_neighbors(unsigned int face_id) const;
  EFAElement3D* get_face_neighbor(unsigned int face_id, unsigned int neighbor_id) const;

  bool frag_has_tip_faces() const;
  std::vector<unsigned int> get_tip_face_id() const;
  std::set<EFANode*> get_tip_embedded_nodes() const;
  bool face_contains_tip(unsigned int face_id) const;
  bool frag_face_already_cut(unsigned int ElemFaceID) const;

  void addFaceEdgeCut(unsigned int face_id, unsigned int edge_id, double position,
                      EFANode* embedded_node, std::map<unsigned int, EFANode*> &EmbeddedNodes,
                      bool add_to_neighbor, bool add_to_adjacent);
  void addFragFaceEdgeCut(unsigned int frag_face_id, unsigned int frag_edge_id, double position,
                          std::map<unsigned int, EFANode*> &EmbeddedNodes, bool add_to_neighbor,
                          bool add_to_adjacent);

private:

  // EFAelement3D specific methods
  void checkNeighborFaceCut(unsigned int face_id, unsigned int edge_id, double position,
                            EFANode* from_node, EFANode* embedded_node, EFANode* &local_embedded);
  void mapParaCoorFrom2Dto3D(unsigned int face_id, std::vector<double> &xi_2d,
                             std::vector<double> &xi_3d) const;
};

#endif
