#ifndef CHSPLIT2CHEMPOT_H
#define CHSPLIT2CHEMPOT_H

#include "Kernel.h"


//Forward Declarations
class CHSplit2ChemPot;

template<>
InputParameters validParams<CHSplit2ChemPot>();

class CHSplit2ChemPot : public Kernel
{
public:

  CHSplit2ChemPot(const std::string & name, InputParameters parameters);
  
protected:
  
  enum PFFunctionType
  {
    Residual,
    OffDiag
  };
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  virtual Real computeDFDC(PFFunctionType type);
  
  unsigned int _c_var;
  VariableValue & _c;

private:
};
#endif //CHSPLIT2CHEMPOT_H
