#ifndef EFFECTIVEELASTICSTRAIN_H
#define EFFECTIVEELASTICSTRAIN_H

#include "AuxKernel.h"

class EffectiveElasticStrain : public AuxKernel
{
public:
  static InputParameters validParams();

  EffectiveElasticStrain(const InputParameters & parameters);
  virtual ~EffectiveElasticStrain() {}

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;

  const MaterialProperty<std::vector<Real> > & _state_var;
};
#endif //EFFECTIVEELASTICSTRAIN_H
