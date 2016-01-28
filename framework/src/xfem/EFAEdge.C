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

#include "EFAEdge.h"

EFAEdge::EFAEdge(EFANode * node1, EFANode * node2):
  _edge_node1(node1),
  _edge_node2(node2)
{
  _embedded_nodes.clear();
  _intersection_x.clear();
  consistency_check();
}

EFAEdge::EFAEdge(const EFAEdge & other_edge)
{
  _edge_node1 = other_edge._edge_node1;
  _edge_node2 = other_edge._edge_node2;
  _intersection_x = other_edge._intersection_x;
  _embedded_nodes = other_edge._embedded_nodes;
  consistency_check();
}

EFAEdge::~EFAEdge() // do not delete edge node - they will be deleted
{}                  // in EFAelement's destructor

bool
EFAEdge::equivalent(const EFAEdge & other) const
{
  bool isEqual = false;
  if (other._edge_node1 == _edge_node1 &&
      other._edge_node2 == _edge_node2)
    isEqual = true;
  else if (other._edge_node2 == _edge_node1 &&
           other._edge_node1 == _edge_node2)
    isEqual = true;
  return isEqual;
}

bool
EFAEdge::isPartialOverlap(const EFAEdge & other) const
{
  bool partially_overlap = false;
  if (containsEdge(other) || other.containsEdge(*this))
    partially_overlap = true;
  return partially_overlap;
}

bool
EFAEdge::containsEdge(const EFAEdge & other) const
{
  bool contains = false;
  if (containsNode(other._edge_node1) &&
      containsNode(other._edge_node2))
    contains = true;
  return contains;
}

bool
EFAEdge::getNodeMasters(EFANode* node, std::vector<EFANode*> &master_nodes,
                        std::vector<double> &master_weights) const
{
  master_nodes.clear();
  master_weights.clear();
  bool masters_found = false;
  if (_edge_node1 == node)
  {
    master_nodes.push_back(_edge_node1);
    master_weights.push_back(1.0);
    masters_found = true;
  }
  else if (_edge_node2 == node)
  {
    master_nodes.push_back(_edge_node2);
    master_weights.push_back(1.0);
    masters_found = true;
  }
  else
  {
    for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
    {
      if (_embedded_nodes[i] == node)
      {
        master_nodes.push_back(_edge_node1);
        master_nodes.push_back(_edge_node2);
        master_weights.push_back(1.0-_intersection_x[i]);
        master_weights.push_back(_intersection_x[i]);
        masters_found = true;
        break;
      }
    } // i
  }
  return masters_found;
}

// TODO: Saving because I don't want to throw it away, but it needs more work to be used.
//bool
//CutElemMesh::edge_t::operator < (const edge_t & other) const
//{
//  node_t * this_min_node;
//  node_t * this_max_node;
//  node_t * other_min_node;
//  node_t * other_max_node;
//
//  int this_node1_unique_index = ((int) edge_node1->category + 1) * edge_node1->id;
//  int this_node2_unique_index = ((int) edge_node2->category + 1) * edge_node2->id;
//  int other_node1_unique_index = ((int) other.edge_node1->category + 1) * edge_node1->id;
//  int other_node2_unique_index = ((int) other.edge_node2->category + 1) * edge_node2->id;
//  int this_min_index = std::min(this_node1_unique_index, this_node2_unique_index);
//  int other_min_index = std::min(other_node1_unique_index, other_node2_unique_index);
//
//  if (this_min_index < other_min_index)
//    return true;
//  else if (this_min_index == other_min_index)
//  {
//    int this_max_index = std::max(this_node1_unique_index, this_node2_unique_index);
//    int other_max_index = std::max(other_node1_unique_index, other_node2_unique_index);
//
//    if (this_max_index < other_max_index)
//      return true;
//  }
//  return false;
//}

void
EFAEdge::add_intersection(double position, EFANode * embedded_node_tmp, EFANode * from_node)
{
  _embedded_nodes.push_back(embedded_node_tmp);
  if (from_node == _edge_node1)
    _intersection_x.push_back(position);
  else if (from_node == _edge_node2)
    _intersection_x.push_back(1.0 - position);
  else
    mooseError("In add_intersection from_node does not exist on edge");
}

void
EFAEdge::reset_intersection(double position, EFANode * embedded_node_tmp, EFANode * from_node)
{
  for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
  {
    if (_embedded_nodes[i] == embedded_node_tmp)
    {
      if (from_node == _edge_node1)
        _intersection_x[i] = position;
      else if (from_node == _edge_node2)
        _intersection_x[i] = 1.0 - position;
      else
        mooseError("In set_intersection from_node does not exist on edge");
      break;
    }
  } // i
}

void
EFAEdge::copy_intersection(const EFAEdge & other, unsigned int from_node_id)
{
  _embedded_nodes.clear();
  _intersection_x.clear();
  _embedded_nodes = other._embedded_nodes;
  if (from_node_id == 0)
    _intersection_x = other._intersection_x;
  else if (from_node_id == 1)
  {
    for (unsigned int i = 0; i < other.num_embedded_nodes(); ++i)
      _intersection_x.push_back(1.0 - other._intersection_x[i]);
  }
  else
    mooseError("from_node_id out of bounds");
  if (_embedded_nodes.size() != _intersection_x.size())
    mooseError("in copy_intersection() num of emb_nodes must be = num of inters_x");
}

EFANode *
EFAEdge::get_node(unsigned int index) const
{
  if (index == 0)
    return _edge_node1;
  else if (index == 1)
    return _edge_node2;
  else
    mooseError("In get_node index out of bounds");
}

void
EFAEdge::reverse_nodes()
{
  EFANode* tmp = _edge_node1;
  _edge_node1 = _edge_node2;
  _edge_node2 = tmp;
  for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
    _intersection_x[i] = 1.0 - _intersection_x[i];
}

bool
EFAEdge::has_intersection() const
{
  return (_embedded_nodes.size() > 0);
}

bool
EFAEdge::has_intersection_at_position(double position, EFANode * from_node) const
{
  double tol = 1.e-4;
  bool has_int = false;
  if (has_intersection())
  {
    double tmp_intersection_x = -1.0;
    if (from_node == _edge_node1)
      tmp_intersection_x = position;
    else if (from_node == _edge_node2)
      tmp_intersection_x = 1.0 - position;
    else
      mooseError("In has_intersection from_node does not exist on edge");

    for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
    {
      if (std::abs(tmp_intersection_x - _intersection_x[i]) < tol)
      {
        has_int = true;
        break;
      }
    } // i
  }
  return has_int;
}

double
EFAEdge::get_intersection(unsigned int emb_id, EFANode * from_node) const
{
  if (from_node == _edge_node1)
    return _intersection_x[emb_id];
  else if (from_node == _edge_node2)
    return 1.0 - _intersection_x[emb_id];
  else
    mooseError("In get_intersection node not in edge");
}

double
EFAEdge::distance_from_node1(EFANode * node) const
{
  double xi = -100.0;
  if (_edge_node1 == node)
    xi = 0.0;
  else if (_edge_node2 == node)
    xi = 1.0;
  else if (is_embedded_node(node))
  {
    unsigned int embedded_node_id = get_embedded_index(node);
    xi = _intersection_x[embedded_node_id];
  }
  else
    mooseError("the given node is not found in the current edge");
  return xi;
}

bool
EFAEdge::is_embedded_node(const EFANode * node) const
{
  bool is_emb = false;
  for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
  {
    if (_embedded_nodes[i] == node)
    {
      is_emb = true;
      break;
    }
  }
  return is_emb;
}

unsigned int
EFAEdge::get_embedded_index(EFANode * node) const
{
  unsigned int index = 99999;
  for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
  {
    if (_embedded_nodes[i] == node)
    {
      index = i;
      break;
    }
  }
  return index;
}

unsigned int
EFAEdge::get_embedded_index(double position, EFANode* from_node) const
{
  unsigned int index = 99999;
  double tol = 1.e-4;
  if (has_intersection())
  {
    double tmp_intersection_x = -1.0; // dist from edge_node1
    if (from_node == _edge_node1)
      tmp_intersection_x = position;
    else if (from_node == _edge_node2)
      tmp_intersection_x = 1.0 - position;
    else
      mooseError("In get_embedded_index(), from_node does not exist on edge");

    for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
    {
      if (std::abs(tmp_intersection_x - _intersection_x[i]) < tol)
      {
        index = i;
        break;
      }
    } // i
  }
  return index;
}

EFANode *
EFAEdge::get_embedded_node(unsigned int index) const
{
  if (index < _embedded_nodes.size())
    return _embedded_nodes[index];
  else
    mooseError("in get_embedded_node() index out of bounds");
}

unsigned int
EFAEdge::num_embedded_nodes() const
{
  return _embedded_nodes.size();
}

void
EFAEdge::consistency_check()
{
  bool consistent = true;
  if ((_edge_node1->category() == N_CATEGORY_PERMANENT ||
       _edge_node1->category() == N_CATEGORY_TEMP) &&
       _edge_node2->category() == N_CATEGORY_LOCAL_INDEX)
    consistent = false;
  else if ((_edge_node2->category() == N_CATEGORY_PERMANENT ||
            _edge_node2->category() == N_CATEGORY_TEMP) &&
            _edge_node1->category() == N_CATEGORY_LOCAL_INDEX)
    consistent = false;
  if (!consistent)
    mooseError("In consistency_check nodes on edge are not consistent");
  if (_embedded_nodes.size() != _intersection_x.size())
    mooseError("In consistency_check num of emb_nodes must be = num of inters_x");
}

void
EFAEdge::switchNode(EFANode *new_node, EFANode *old_node)
{
  if (_edge_node1 == old_node)
    _edge_node1 = new_node;
  else if (_edge_node2 == old_node)
    _edge_node2 = new_node;
  else if (is_embedded_node(old_node))
  {
    unsigned int id = get_embedded_index(old_node);
    _embedded_nodes[id] = new_node;
  }
}

bool
EFAEdge::containsNode(const EFANode *node) const
{
  if (_edge_node1 == node || _edge_node2 == node ||
      is_embedded_node(node))
    return true;
  else
    return false;
}

void
EFAEdge::remove_embedded_node()
{
  _embedded_nodes.clear();
  _intersection_x.clear();
}

void
EFAEdge::remove_embedded_node(EFANode * node)
{
  unsigned int index = 0;
  bool node_found = false;
  for (unsigned int i = 0; i < _embedded_nodes.size(); ++i)
  {
    if (_embedded_nodes[i] == node)
    {
      index = i;
      node_found = true;
      break;
    }
  }// i
  if (node_found)
  {
    _embedded_nodes.erase(_embedded_nodes.begin() + index);
    _intersection_x.erase(_intersection_x.begin() + index);
  }
}
