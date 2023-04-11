#ifndef ELEMENTNORMALIZEDVALUE_H
#define ELEMENTNORMALIZEDVALUE_H

#include "ElementIntegralVariablePostprocessor.h"

class ElementNormalizedValue : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  ElementNormalizedValue(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  Real _volume;
};

#endif
