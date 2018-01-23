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

#include "VectorPostprocessorComponent.h"
#include "VectorPostprocessorInterface.h"

template <>
InputParameters
validParams<VectorPostprocessorComponent>()
{
  InputParameters params = validParams<GeneralPostprocessor>();

  params.addRequiredParam<VectorPostprocessorName>("vectorpostprocessor",
                                                   "The vectorpostprocessor from which a value is extracted");
  params.addRequiredParam<std::string>("vector_name", "Name of the vector for which to report a value");
  params.addRequiredParam<unsigned int>("index", "Index of the vectorpostprocessor for which to report a value");
  params.addClassDescription("Reports the value of the specified component of another VectorPostprocessor");

  return params;
}

VectorPostprocessorComponent::VectorPostprocessorComponent(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
  _vpp_name(getParam<VectorPostprocessorName>("vectorpostprocessor")),
  _vector_name(getParam<std::string>("vector_name")),
  _vpp_values(getVectorPostprocessorValue("vectorpostprocessor", _vector_name)),
  _vpp_index(getParam<unsigned int>("index"))
{
}

Real
VectorPostprocessorComponent::getValue()
{
  if (_vpp_index >= _vpp_values.size())
  {
    std::cout<<"BWS vpp: "<<_vpp_name<<" vname: "<<_vector_name<<std::endl;
    mooseError("In VectorPostprocessorComponent index greater than size of vector");
  }
  return _vpp_values[_vpp_index];
}
