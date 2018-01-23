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

#ifndef VECTORPOSTPROCESSORCOMPONENT_H
#define VECTORPOSTPROCESSORCOMPONENT_H

#include "GeneralPostprocessor.h"

// Forward Declarations
class VectorPostprocessorComponent;

template <>
InputParameters validParams<VectorPostprocessorComponent>();

class VectorPostprocessorComponent : public GeneralPostprocessor
{
public:
  VectorPostprocessorComponent(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}

  virtual Real getValue() override;

protected:
  const VectorPostprocessorName _vpp_name;
  const std::string _vector_name;
  const VectorPostprocessorValue & _vpp_values;
  const unsigned int _vpp_index;
};

#endif // VECTORPOSTPROCESSORCOMPONENT_H
