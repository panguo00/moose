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

#ifndef VECTOR_OF_POSTPROCESSORS_H
#define VECTOR_OF_POSTPROCESSORS_H

#include "GeneralVectorPostprocessor.h"

//Forward Declarations
class VectorOfPostprocessors;

template<>
InputParameters validParams<VectorOfPostprocessors>();

class VectorOfPostprocessors :
  public GeneralVectorPostprocessor
{
public:
  VectorOfPostprocessors(const std::string & name, InputParameters parameters);

  virtual ~VectorOfPostprocessors() {}

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

  virtual void threadJoin(const UserObject & y);

protected:
  VectorPostprocessorValue & _pp_vec;
  std::vector<const PostprocessorValue *> _postprocessor_values;
};

#endif
