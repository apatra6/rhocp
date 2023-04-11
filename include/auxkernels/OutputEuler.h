#ifndef OUTPUTEULER_H
#define OUTPUTEULER_H

#include "AuxKernel.h"

class OutputEuler : public AuxKernel
{
public:
  static InputParameters validParams();

  OutputEuler(const InputParameters & parameters);
  virtual ~OutputEuler() {}

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;

  unsigned int _angle_id;

  const MaterialProperty<Point> & _euler_ang;
};
#endif //OUTPUTEULER_H
