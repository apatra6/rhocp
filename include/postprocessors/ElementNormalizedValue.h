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
  virtual Real getValue() const override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

protected:
  /// volume of the integration domain
  Real _volume;
  /// the integral value that is being accumulated
  Real _integral_value;
};

#endif
