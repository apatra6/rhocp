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
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  unsigned int _grain_id;
  const EBSDMeshReader * _EBSDFileReader;
  Real _volume;

};

#endif
