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

#ifndef ELEMENTFRAGMENTALGORITHM_H
#define ELEMENTFRAGMENTALGORITHM_H

#include <cstddef>
#include <iostream>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <limits>

#include "EFAnode.h"
#include "FaceNode.h"
#include "EfaEdge.h"
#include "EfaElement2D.h"
#include "EfaElement3D.h"
#include "EFAfragment.h"

class ElementFragmentAlgorithm
{
public:

  /**
   * Constructor
   **/
  ElementFragmentAlgorithm();

  ~ElementFragmentAlgorithm();

private:
  //unsigned int MaxElemId;
  std::map< unsigned int, EFAnode*> _permanent_nodes;
  std::map< unsigned int, EFAnode*> _embedded_nodes;
  std::map< unsigned int, EFAnode*> _temp_nodes;
  std::map< unsigned int, EfaElement*> _elements;
//  std::map< std::set< EFAnode* >, std::set< EFAelement* > > _merged_edge_map;
  std::set< EfaElement*> _crack_tip_elements;
  std::vector< EFAnode* > _new_nodes;
  std::vector< EfaElement* > _child_elements;
  std::vector< EfaElement* > _parent_elements;
  std::map< EFAnode*, std::set< EfaElement *> > _inverse_connectivity;

public:

  unsigned int add2DElements( std::vector< std::vector<unsigned int> > &quads );
  EfaElement* add2DElement( std::vector<unsigned int> quad, unsigned int id );
  EfaElement* add3DElement( std::vector<unsigned int> quad, unsigned int id );

  void updateEdgeNeighbors();
  void initCrackTipTopology();
  void addElemEdgeIntersection(unsigned int elemid, unsigned int edgeid, double position);
  void addFragEdgeIntersection(unsigned int elemid, unsigned int frag_edge_id, double position);
  void addElemFaceIntersection(unsigned int elemid, unsigned int faceid,
                               std::vector<unsigned int> edgeid, std::vector<double> position);
  void addFragFaceIntersection(unsigned int ElemID, unsigned int FragFaceID,
                               std::vector<unsigned int> FragFaceEdgeID, std::vector<double> position);

  void updatePhysicalLinksAndFragments();

  void updateTopology(bool mergeUncutVirtualEdges=true);
  void reset();
  void clearAncestry();
  void restoreFragmentInfo(EfaElement * const elem, const EfaElement * const from_elem);

  void createChildElements();
  void connectFragments(bool mergeUncutVirtualEdges);

  void sanityCheck();
  void updateCrackTipElements();
  void printMesh();
  void error(const std::string &error_string);

  const std::vector<EfaElement*> &getChildElements(){return _child_elements;};
  const std::vector<EfaElement*> &getParentElements(){return _parent_elements;};
  const std::vector<EFAnode*> &getNewNodes(){return _new_nodes;};
  const std::set<EfaElement*> &getCrackTipElements(){return _crack_tip_elements;};
  const std::map<unsigned int, EFAnode*> &getPermanentNodes(){return _permanent_nodes;}
  const std::map<unsigned int, EFAnode*> &getTempNodes(){return _temp_nodes;}
  const std::map<unsigned int, EFAnode*> &getEmbeddedNodes(){return _embedded_nodes;}
  EfaElement* getElemByID(unsigned int id);
  unsigned int getElemIdByNodes(unsigned int * node_id);
  void clearPotentialIsolatedNodes();
};

#endif // #ifndef ELEMENTFRAGMENTALGORITHM_H
