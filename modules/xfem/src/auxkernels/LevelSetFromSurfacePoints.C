/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "LevelSetFromSurfacePoints.h"
#include "DelimitedFileReader.h"
#include "Conversion.h"

template <>
InputParameters
validParams<LevelSetFromSurfacePoints>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredParam<FileName>("point_data_file", "Name of the CSV file containing the point coordinate data");
  params.addRequiredParam<FileName>("center_data_file", "Name of the CSV file containing the object center coordinate data");
  params.addRequiredParam<std::string>("object_id_header", "Header name for the column in the CSV file with the object IDs");
  params.addRequiredParam<std::string>("x_header", "Header name for the column in the CSV file with the x coordinate");
  params.addRequiredParam<std::string>("y_header", "Header name for the column in the CSV file with the y coordinate");

  return params;
}

LevelSetFromSurfacePoints::LevelSetFromSurfacePoints(const InputParameters & parameters) : AuxKernel(parameters),
  _point_data_file_name(getParam<FileName>("point_data_file")),
  _center_data_file_name(getParam<FileName>("center_data_file")),
  _object_id_header(getParam<std::string>("object_id_header")),
  _x_header(getParam<std::string>("x_header")),
  _y_header(getParam<std::string>("y_header"))
{
  if (!isNodal())
    mooseError("LevelSetFromSurfacePoints can only be run on a nodal variable");
}

void
LevelSetFromSurfacePoints::initialSetup()
{
  // Read and store the point coordinate data for each object
  MooseUtils::DelimitedFileReader pd_csv_reader(_point_data_file_name);
  pd_csv_reader.setDelimiter(",");
  pd_csv_reader.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::ON);
  pd_csv_reader.read();
  std::vector<Real> object_id_real = pd_csv_reader.getData(_object_id_header);

  std::vector<int> object_id(object_id_real.size());
  for (unsigned int i = 0; i < object_id_real.size(); ++i)
    object_id[i] = static_cast<int>(object_id_real[i]);

  std::vector<Real> x  = pd_csv_reader.getData(_x_header);
  std::vector<Real> y  = pd_csv_reader.getData(_y_header);

  if ((object_id.size() != x.size()) ||
      (object_id.size() != y.size()))
    mooseError("Columns of data in point data csv file must all be the same size");

  for (unsigned int i = 0; i < object_id.size(); ++i)
    _point_data[object_id[i]].push_back(Point(x[i], y[i], 0.0));

  // Read and store the center coordinate data for each object
  MooseUtils::DelimitedFileReader oc_csv_reader(_center_data_file_name);
  oc_csv_reader.setDelimiter(",");
  oc_csv_reader.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::ON);
  oc_csv_reader.read();
  object_id_real = oc_csv_reader.getData(_object_id_header);

  object_id.resize(object_id_real.size());
  if (_point_data.size() != object_id.size())
    mooseError("Number of objects for point data and center data must match");

  for (unsigned int i = 0; i < object_id_real.size(); ++i)
  {
    object_id[i] = static_cast<int>(object_id_real[i]);
    if (_point_data.find(object_id[i]) == _point_data.end())
      mooseError("No point data for object " + Moose::stringify(object_id[i]) + " in center data file");
  }

  x  = oc_csv_reader.getData(_x_header);
  y  = oc_csv_reader.getData(_y_header);

  if ((object_id.size() != x.size()) ||
      (object_id.size() != y.size()))
    mooseError("Columns of data in object center csv file must all be the same size");

  for (unsigned int i = 0; i < object_id.size(); ++i)
    _object_centers[object_id[i]] = Point(x[i], y[i], 0.0);


//  Real real_max = std::numeric_limits<Real>::max();
//  for (auto & pd : _point_data)
//  {
//    Point min(real_max, real_max, 0.0);
//    Point max(-real_max, -real_max, 0.0);
//    for (auto & point : pd.second)
//    {
//      if (point(0) < min(0))
//        min(0) = point(0);
//      if (point(1) < min(1))
//        min(1) = point(1);
//      if (point(0) > max(0))
//        max(0) = point(0);
//      if (point(1) > max(1))
//        max(1) = point(1);
//    }
//    Point center = (min + max) / 2.0;
//    std::cout<<"BWS min: "<<min(0)<<" "<<min(1)<<std::endl;
//    std::cout<<"BWS max: "<<max(0)<<" "<<max(1)<<std::endl;
//    std::cout<<"BWS center: "<<center(0)<<" "<<center(1)<<std::endl;
//    _object_centers[pd.first] = center;
//  }
}

Real
LevelSetFromSurfacePoints::computeValue()
{
  Real min_dist = std::numeric_limits<Real>::max();

  for (auto & pd : _point_data)
  {
    Real this_obj_min_dist = minimumDistanceForObject(pd.second, _object_centers[pd.first]);
    min_dist = std::min(this_obj_min_dist, min_dist);
  }

  return min_dist;
}


Real LevelSetFromSurfacePoints::minimumDistanceForObject(const std::vector<Point> & point_data_this_obj, const Point & center_this_obj)
{
  Real min_dist = std::numeric_limits<Real>::max();

  for (const auto & point : point_data_this_obj)
  {
    Point point_to_node = *_current_node - point;
    Real signed_distance = point_to_node.norm();
    if (signed_distance < std::abs(min_dist))
    {
      Point center_to_point = point - center_this_obj;
      Real dot_prod = point_to_node * center_to_point;
      if (dot_prod < 0.0)
        signed_distance *= -1.0;
      min_dist = signed_distance;
    }
  }

  return min_dist;
}
