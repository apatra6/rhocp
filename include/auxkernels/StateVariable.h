#ifndef STATEVARIABLE_H
#define STATEVARIABLE_H

#include "AuxKernel.h"

class StateVariable : public AuxKernel
{
public:
  static InputParameters validParams();

  StateVariable(const InputParameters & parameters);
  virtual ~StateVariable() {}

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;

  unsigned int _sdv_id;

  const MaterialProperty<std::vector<Real> > & _state_var;
};
#endif //STATEVARIABLE_H
