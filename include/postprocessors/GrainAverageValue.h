#ifndef GRAINAVERAGEVALUE_H
#define GRAINAVERAGEVALUE_H

#include "ElementIntegralVariablePostprocessor.h"
#include "EBSDMeshReader.h"

class GrainAverageValue : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  GrainAverageValue(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() const override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

protected:
  unsigned int _grain_id;
  const EBSDMeshReader * _EBSDFileReader;
  /// volume of the integration domain
  Real _volume;
  /// the integral value that is being accumulated
  Real _integral_value;

};

#endif
