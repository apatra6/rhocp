#ifndef LATTICESTRAIN_H
#define LATTICESTRAIN_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class LatticeStrain : public AuxKernel
{
public:
  static InputParameters validParams();

  LatticeStrain(const InputParameters & parameters);
  virtual ~LatticeStrain() {}

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;

  Real _h, _k, _l, _g0, _g1, _g2;

  Real _tol;

  const MaterialProperty<std::vector<Real> > & _state_var;
  const MaterialProperty<Point> & _euler_ang;

  // The eigenstrains
  std::vector<MaterialPropertyName> _eigenstrain_names;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains;
};
#endif //LATTICESTRAIN_H
