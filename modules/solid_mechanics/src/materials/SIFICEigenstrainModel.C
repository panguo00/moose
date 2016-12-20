/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "SIFICEigenstrainModel.h"
#include "SymmTensor.h"

template<>
InputParameters validParams<SIFICEigenstrainModel>()
{
  return validParams<Material>();
}

SIFICEigenstrainModel::SIFICEigenstrainModel(const InputParameters & parameters ):
  VolumetricModel( parameters )
{}

SIFICEigenstrainModel::~SIFICEigenstrainModel() {}

void
SIFICEigenstrainModel::modifyStrain(const unsigned int qp,
                                    const Real scale_factor,
                                    SymmTensor & strain_increment,
                                    SymmTensor & /*dstrain_increment_dT*/)
{
  Real ypos = _q_point[qp](1);
  strain_increment.xx() -= ypos * ypos / 1e4;
}
