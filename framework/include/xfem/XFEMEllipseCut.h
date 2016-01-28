
#ifndef XFEM_ELLIPSE_CUT_H
#define XFEM_ELLIPSE_CUT_H

#include "XFEMGeometricCut.h"

class XFEMEllipseCut : public XFEMGeometricCut
{
public:

  XFEMEllipseCut(std::vector<Real> square_nodes);
  ~XFEMEllipseCut();

  virtual bool cut_elem_by_geometry(const Elem* elem, std::vector<cutEdge> & cutEdges, Real time);
  virtual bool cut_elem_by_geometry(const Elem* elem, std::vector<cutFace> & cutFaces, Real time);

  virtual bool cut_frag_by_geometry(std::vector<std::vector<Point> > & frag_edges, std::vector<cutEdge> & cutEdges, Real time);
  virtual bool cut_frag_by_geometry(std::vector<std::vector<Point> > & frag_faces, std::vector<cutFace> & cutFaces, Real time);

private:

  std::vector<Point> _vertices;
  Point _center;
  Point _normal;
  Point _unit_vec1;
  Point _unit_vec2;
  Real _long_axis;
  Real _short_axis;
//  Real _angle;

private:

  bool intersect_with_edge(Point p1, Point p2, Point &pint);
  bool isInsideCutPlane(Point p);
};
  
#endif

