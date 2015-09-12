/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef FRICTIONALCONTACTDAMPERPROBLEM_H
#define FRICTIONALCONTACTDAMPERPROBLEM_H

#include "ReferenceResidualProblem.h"

class FrictionalContactDamperProblem;

template<>
InputParameters validParams<FrictionalContactDamperProblem>();

/**
 * FEProblem derived class for frictional contact-specific callbacks
 */
class FrictionalContactDamperProblem : public ReferenceResidualProblem
{
public:
  FrictionalContactDamperProblem(const InputParameters & params);
  virtual ~FrictionalContactDamperProblem();

  virtual void initialSetup();
  virtual void timestepSetup();
  virtual bool shouldUpdateSolution();
  virtual bool updateSolution(NumericVector<Number>& vec_solution, NumericVector<Number>& ghosted_solution);
  bool limitSlip(NumericVector<Number>& vec_solution, NumericVector<Number>& ghosted_solution);
  unsigned int numLocalFrictionalConstraints();
  void updateContactPoints(NumericVector<Number>& ghosted_solution,
                           bool update_incremental_slip);

protected:
  std::set<std::pair<int,int> > _interactions;
  NonlinearVariableName _disp_x;
  NonlinearVariableName _disp_y;
  NonlinearVariableName _disp_z;

  int _num_nl_iterations;
  int _num_contact_nodes;
  int _num_sticking;
  int _num_slipping;
  int _num_slipping_friction;
  int _num_slip_reversed;
  int _num_modified;
  Real _inc_slip_norm;
  Real _it_slip_norm;
  Real _max_iterative_slip;

  /// Convenient typedef for frequently used iterator
  typedef std::map<std::pair<unsigned int, unsigned int>, PenetrationLocator *>::iterator pl_iterator;
};

#endif /* FRICTIONALCONTACTDAMPERPROBLEM_H */
