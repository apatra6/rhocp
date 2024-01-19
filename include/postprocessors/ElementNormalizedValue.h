#pragma once

#include "ElementIntegralVariablePostprocessor.h"

class ElementNormalizedValue : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  ElementNormalizedValue(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() const override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  Real _volume;
};

