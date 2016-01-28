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

#ifndef EFAFRAGMENT3D_H
#define EFAFRAGMENT3D_H

#include "EFAEdge.h"
#include "EFAFace.h"
#include "EFAFragment.h"

class EFAElement3D;

class EFAFragment3D : public EFAFragment
{
public:

  EFAFragment3D(EFAElement3D* host, bool create_faces, const EFAElement3D * from_host,
                unsigned int frag_id = std::numeric_limits<unsigned int>::max());
  ~EFAFragment3D();

private:

  EFAElement3D * _host_elem;
  std::vector<EFAFace*> _faces;
  std::vector<std::vector<EFAFace*> > _adjacent_face_ix;

public:
  // override pure virtual methods
  virtual void switchNode(EFANode *new_node, EFANode *old_node);
  virtual bool containsNode(EFANode *node) const;
  virtual unsigned int get_num_cuts() const;
  virtual std::set<EFANode*> get_all_nodes() const;
  virtual bool isConnected(EFAFragment *other_fragment) const;
  virtual void remove_invalid_embedded(std::map<unsigned int, EFANode*> &EmbeddedNodes);

  // EFAfragment3D specific methods
  void combine_tip_faces();
  bool is_face_interior(unsigned int face_id) const;
  std::vector<unsigned int> get_interior_face_id() const;
  bool isThirdInteriorFace(unsigned int face_id) const;
  unsigned int num_faces() const;
  EFAFace* get_face(unsigned int face_id) const;
  unsigned int get_face_id(EFAFace* face) const;
  void add_face(EFAFace* new_face);
  std::set<EFANode*> get_face_nodes(unsigned int face_id) const;
  EFAElement3D * get_host() const;
  std::vector<EFAFragment3D*> split();
  void create_adjacent_face_ix();
  EFAFace* get_adjacent_face(unsigned int face_id, unsigned int edge_id) const;
  void remove_embedded_node(EFANode* emb_node);
  bool hasFaceWithOneCut() const;
  void get_node_info(std::vector<std::vector<unsigned int> > &face_node_ix,
                     std::vector<EFANode*> &nodes) const;

private:

  EFAFragment3D* connect_subfaces(EFAFace* start_face, unsigned int startOldFaceID,
                                  std::vector<std::vector<EFAFace*> > &subfaces);
  EFAEdge* lonelyEdgeOnFace(unsigned int face_id) const;
  void combine_two_faces(unsigned int face_id1, unsigned int face_id2, const EFAFace* elem_face);
};

#endif
