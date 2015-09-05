/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#include "FrictionalContactDamperProblem.h"

#include "NonlinearSystem.h"
#include "DisplacedProblem.h"
#include "PenetrationLocator.h"
#include "NearestNodeLocator.h"
#include "MooseApp.h"

#include <limits>

template<>
InputParameters validParams<FrictionalContactDamperProblem>()
{
  InputParameters params = validParams<ReferenceResidualProblem>();
  params.addRequiredParam<std::vector<int> >("master","IDs of the master surfaces for which the slip should be calculated");
  params.addRequiredParam<std::vector<int> >("slave","IDs of the slave surfaces for which the slip should be calculated");
  params.addRequiredParam<NonlinearVariableName>("disp_x","Variable containing the x displacement");
  params.addRequiredParam<NonlinearVariableName>("disp_y","Variable containing the y displacement");
  params.addParam<NonlinearVariableName>        ("disp_z","Variable containing the z displacement");
  return params;
}

FrictionalContactDamperProblem::FrictionalContactDamperProblem(const InputParameters & params) :
    ReferenceResidualProblem(params),
    _num_contact_nodes(0),
    _num_slip_reversed(0),
    _num_modified(0),
    _inc_slip_norm(0.0),
    _it_slip_norm(0.0)
{
  std::vector<int> master = params.get<std::vector<int> >("master");
  std::vector<int> slave = params.get<std::vector<int> >("slave");

  unsigned int dim = getNonlinearSystem().subproblem().mesh().dimension();

  _disp_x = params.get<NonlinearVariableName>("disp_x");
  _disp_y = params.get<NonlinearVariableName>("disp_y");

  if (dim == 3)
  {
    if (!params.isParamValid("disp_z"))
      mooseError("Missing disp_z in FrictionalContactDamperProblem");
    _disp_z = params.get<NonlinearVariableName>("disp_z");
  }

  unsigned int num_interactions = master.size();
  if (num_interactions != slave.size())
    mooseError("Sizes of master surface and slave surface lists must match in FrictionalContactDamperProblem");

  for (unsigned int i=0; i<master.size(); ++i)
  {
    std::pair<int,int> ms_pair(master[i],slave[i]);
    _interactions.insert(ms_pair);
  }
}

FrictionalContactDamperProblem::~FrictionalContactDamperProblem()
{}

void
FrictionalContactDamperProblem::initialSetup()
{
  ReferenceResidualProblem::initialSetup();
}

void
FrictionalContactDamperProblem::timestepSetup()
{
  _num_nl_iterations = 0;
  ReferenceResidualProblem::timestepSetup();
}

bool
FrictionalContactDamperProblem::shouldUpdateSolution()
{
  return true;
}

bool
FrictionalContactDamperProblem::updateSolution(NumericVector<Number>& vec_solution, NumericVector<Number>& ghosted_solution)
{
  bool solution_modified = false;

  unsigned int nfc = numLocalFrictionalConstraints();
  unsigned int dim = getNonlinearSystem().subproblem().mesh().dimension();
  updateContactPoints(ghosted_solution,false);

  _console << "Slip Update: " << _num_nl_iterations << std::endl;
  _console << "Iter  #Cont     #SlipRev     #Mod" << std::endl;

  solution_modified = limitSlip(vec_solution, ghosted_solution);

    _console << std::setw(10) << _num_nl_iterations
             << std::setw(10) << _num_contact_nodes
             << std::setw(10) << _num_slip_reversed
             << std::setw(10) << _num_modified
             << std::endl;

  _num_nl_iterations++;

  return solution_modified;
}

bool
FrictionalContactDamperProblem::limitSlip(NumericVector<Number>& vec_solution, NumericVector<Number>& ghosted_solution)
{
  NonlinearSystem & nonlinear_sys = getNonlinearSystem();
  unsigned int dim = nonlinear_sys.subproblem().mesh().dimension();

  MooseVariable * disp_x_var = &getVariable(0,_disp_x);
  MooseVariable * disp_y_var = &getVariable(0,_disp_y);
  MooseVariable * disp_z_var = NULL;
  if (dim == 3)
    disp_z_var = &getVariable(0,_disp_z);

  _it_slip_norm = 0.0;
  _inc_slip_norm = 0.0;
  TransientNonlinearImplicitSystem & system = getNonlinearSystem().sys();

  if (getDisplacedProblem() && _interactions.size() > 0)
  {
    _num_contact_nodes = 0;
    _num_slip_reversed = 0;
    _num_modified = 0;

    GeometricSearchData & displaced_geom_search_data = getDisplacedProblem()->geomSearchData();
    std::map<std::pair<unsigned int, unsigned int>, PenetrationLocator *> * penetration_locators = &displaced_geom_search_data._penetration_locators;

    AuxiliarySystem & aux_sys = getAuxiliarySystem();
    const NumericVector<Number> & aux_solution = *aux_sys.currentSolution();

    for (pl_iterator plit = penetration_locators->begin(); plit != penetration_locators->end(); ++plit)
    {
      PenetrationLocator & pen_loc = *plit->second;

      bool frictional_contact_this_interaction = false;

      std::set<std::pair<int,int> >::iterator ipit;
      std::pair<int,int> ms_pair(pen_loc._master_boundary,pen_loc._slave_boundary);
      ipit = _interactions.find(ms_pair);
      if (ipit != _interactions.end())
        frictional_contact_this_interaction = true;

      if (frictional_contact_this_interaction)
      {
        std::vector<dof_id_type> & slave_nodes = pen_loc._nearest_node._slave_nodes;

        for (unsigned int i=0; i<slave_nodes.size(); i++)
        {
          dof_id_type slave_node_num = slave_nodes[i];

          if (pen_loc._penetration_info[slave_node_num])
          {
            PenetrationInfo & info = *pen_loc._penetration_info[slave_node_num];
            const Node * node = info._node;

            if (node->processor_id() == processor_id())
            {


              if (info.isCaptured())
              {
                _num_contact_nodes++;
                
                RealVectorValue slip_correction = 0;
                RealVectorValue tangential_inc_slip_prev_iter = info._incremental_slip_prev_iter - (info._incremental_slip_prev_iter * info._normal) * info._normal;
                RealVectorValue tangential_inc_slip = info._incremental_slip - (info._incremental_slip * info._normal) * info._normal;
                RealVectorValue tangential_iterative_slip = tangential_inc_slip - tangential_inc_slip_prev_iter;
                RealVectorValue iterative_slip = info._incremental_slip - info._incremental_slip_prev_iter;
                if (_num_nl_iterations == 0)
                  info._slip_reversed = false;
                else
                {
                  //if (info._incremental_slip_prev_iter * info._incremental_slip < 0)
                  if (tangential_inc_slip_prev_iter * tangential_inc_slip < 0)
                  //if (false)
                  {
                    info._slip_reversed = true;
                    _num_slip_reversed++;
                    //slip_correction = -info._incremental_slip;
//                    slip_correction = -tangential_inc_slip;
//                    std::cout<<"BWS prev: "<<info._incremental_slip_prev_iter<<std::endl;
//                    std::cout<<"BWS curr: "<<info._incremental_slip<<std::endl;
//                    std::cout<<"BWS correction: "<<slip_correction<<std::endl;
                  }
                  else
                  {
                    Real iter_slip_mult = 1-(Real)_num_nl_iterations/50;
                    if (iter_slip_mult < 0)
                      iter_slip_mult = 0;
//                    std::cout<<"BWS iter_slip_mult: "<<iter_slip_mult<<std::endl;
                    //slip_correction = -iter_slip_mult*info._incremental_slip;
                    //slip_correction = -iter_slip_mult*iterative_slip;
                    //slip_correction = -0.5*iterative_slip;
                    //slip_correction = -iter_slip_mult*tangential_iterative_slip;
                  }
                  slip_correction = -0.1*iterative_slip;
                }
                info._incremental_slip_prev_iter = info._incremental_slip;

                VectorValue<dof_id_type> solution_dofs(node->dof_number(nonlinear_sys.number(), disp_x_var->number(), 0),
                                                       node->dof_number(nonlinear_sys.number(), disp_y_var->number(), 0),
                                                       (disp_z_var ? node->dof_number(nonlinear_sys.number(), disp_z_var->number(), 0) : 0));

                if (slip_correction.size() > 0)
                {
                  _num_modified++;
                  for (unsigned int i=0; i<dim; ++i)
                    vec_solution.add(solution_dofs(i), slip_correction(i));
                }
              }
            }
          }
        }
      }
    }
    _communicator.sum(_num_contact_nodes);
    _communicator.sum(_num_slip_reversed);
    _communicator.sum(_num_modified);
  }

  vec_solution.close();
  bool updated_solution = false;
  if (_num_modified > 0)
  {
    updated_solution = true;
    ghosted_solution = vec_solution;
    ghosted_solution.close();

//    updateContactPoints(ghosted_solution,false);

//    enforceRateConstraint(vec_solution, ghosted_solution);
  }
//    _communicator.sum(_it_slip_norm);
//    _it_slip_norm = std::sqrt(_it_slip_norm);

  return updated_solution;
}

unsigned int
FrictionalContactDamperProblem::numLocalFrictionalConstraints()
{
  GeometricSearchData & displaced_geom_search_data = getDisplacedProblem()->geomSearchData();
  std::map<std::pair<unsigned int, unsigned int>, PenetrationLocator *> * penetration_locators = &displaced_geom_search_data._penetration_locators;

  unsigned int num_constraints(0);

  for (pl_iterator plit = penetration_locators->begin(); plit != penetration_locators->end(); ++plit)
  {
    PenetrationLocator & pen_loc = *plit->second;

    bool frictional_contact_this_interaction = false;

    std::set<std::pair<int,int> >::iterator ipit;
    std::pair<int,int> ms_pair(pen_loc._master_boundary,pen_loc._slave_boundary);
    ipit = _interactions.find(ms_pair);
    if (ipit != _interactions.end())
      frictional_contact_this_interaction = true;

    if (frictional_contact_this_interaction)
    {
      std::vector<dof_id_type> & slave_nodes = pen_loc._nearest_node._slave_nodes;

      for (unsigned int i=0; i<slave_nodes.size(); i++)
      {
        dof_id_type slave_node_num = slave_nodes[i];

        PenetrationInfo * pinfo = pen_loc._penetration_info[slave_node_num];
        if (pinfo)
          if (pinfo->isCaptured())
            ++num_constraints;
      }
    }
  }
  return num_constraints;
}

void
FrictionalContactDamperProblem::updateContactPoints(NumericVector<Number>& ghosted_solution,
                                              bool update_incremental_slip)
{
  GeometricSearchData & displaced_geom_search_data = getDisplacedProblem()->geomSearchData();
  std::map<std::pair<unsigned int, unsigned int>, PenetrationLocator *> * penetration_locators = &displaced_geom_search_data._penetration_locators;

//  for (pl_iterator plit = penetration_locators->begin(); plit != penetration_locators->end(); ++plit)
//  {
//    PenetrationLocator & pen_loc = *plit->second;
//    pen_loc.setUpdate(true);
//  }

  //Do new contact search to update positions of slipped nodes
  _displaced_problem->updateMesh(ghosted_solution, *_aux.currentSolution());

//  for (pl_iterator plit = penetration_locators->begin(); plit != penetration_locators->end(); ++plit)
//  {
//    PenetrationLocator & pen_loc = *plit->second;
//    pen_loc.setUpdate(false);
//  }
}
