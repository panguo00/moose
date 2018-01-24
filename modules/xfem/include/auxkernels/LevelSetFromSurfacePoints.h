/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef LEVELSETFROMSURFACEPOINTS_H
#define LEVELSETFROMSURFACEPOINTS_H

#include "AuxKernel.h"

class LevelSetFromSurfacePoints : public AuxKernel
{
public:
  LevelSetFromSurfacePoints(const InputParameters & parameters);

  virtual ~LevelSetFromSurfacePoints() {}

  virtual void initialSetup() override;
  virtual Real computeValue() override;

protected:
  Real minimumDistanceForObject(const std::vector<Point> & point_data_this_obj,
                                const Point & center_this_obj);

  FileName _point_data_file_name;
  FileName _center_data_file_name;
  std::string _object_id_header;
  std::string _x_header;
  std::string _y_header;

  std::map<int, std::vector<Point>> _point_data;
  std::map<int, Point> _object_centers;

};

template <>
InputParameters validParams<LevelSetFromSurfacePoints>();

#endif // LEVELSETFROMSURFACEPOINTS_H
