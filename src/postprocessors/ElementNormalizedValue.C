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
ElementNormalizedValue::getValue() const
{
  return _integral_value / _volume;
}

void
ElementNormalizedValue::finalize()
{
  gatherSum(_integral_value);
  gatherSum(_volume);
}

void
ElementNormalizedValue::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  const ElementNormalizedValue & pps = static_cast<const ElementNormalizedValue &>(y);

//   const auto & pps = static_cast<const ElementNormalizedValue &>(y);

  _integral_value += pps._integral_value;
  _volume += pps._volume;
}
