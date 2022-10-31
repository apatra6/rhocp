#ifndef GRAINAVERAGEVALUE_H
#define GRAINAVERAGEVALUE_H

#include "ElementIntegralVariablePostprocessor.h"
#include "EBSDMeshReader.h"

// Forward Declarations
class GrainAverageValue;

template <>
InputParameters validParams<GrainAverageValue>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class GrainAverageValue : public ElementIntegralVariablePostprocessor
{
public:
  GrainAverageValue(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  unsigned int _grain_id;
  const EBSDMeshReader * _EBSDFileReader;
  Real _volume;

};

#endif
