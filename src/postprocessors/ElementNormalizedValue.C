#include "ElementNormalizedValue.h"

registerMooseObject("RhocpApp", ElementNormalizedValue);

InputParameters
ElementNormalizedValue::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  return params;
}

ElementNormalizedValue::ElementNormalizedValue(const InputParameters & parameters)
  : ElementIntegralVariablePostprocessor(parameters), _volume(0)
{
}

void
ElementNormalizedValue::initialize()
{
  ElementIntegralVariablePostprocessor::initialize();
  _volume = 0;
}

void
ElementNormalizedValue::execute()
{
  ElementIntegralVariablePostprocessor::execute();

  if (abs(computeIntegral()) >= 1.0e-15) {
  _volume += _current_elem_volume;
  }
}

Real
ElementNormalizedValue::getValue()
{
  Real integral = ElementIntegralVariablePostprocessor::getValue();

  gatherSum(_volume);

  if (_volume == 0.0e0) {
      _volume = 1.0e0;
  }

  return integral / _volume;
}

void
ElementNormalizedValue::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  const ElementNormalizedValue & pps = static_cast<const ElementNormalizedValue &>(y);
  _volume += pps._volume;
}
