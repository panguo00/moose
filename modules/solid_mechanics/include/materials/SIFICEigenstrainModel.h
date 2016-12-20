/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef SIFICEIGENSTRAINMODEL_H
#define SIFICEIGENSTRAINMODEL_H

#include "VolumetricModel.h"

class SymmTensor;

class SIFICEigenstrainModel;

template<>
InputParameters validParams<SIFICEigenstrainModel>();

class SIFICEigenstrainModel : public VolumetricModel
{
public:
  SIFICEigenstrainModel( const InputParameters & parameters );
  virtual ~SIFICEigenstrainModel();

  virtual void modifyStrain(const unsigned int qp,
                            const Real scale_factor,
                            SymmTensor & strain_increment,
                            SymmTensor & dstrain_increment_dT);

};

#endif // SIFICEIGENSTRAINMODEL_H
