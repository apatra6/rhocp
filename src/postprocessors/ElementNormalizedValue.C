// This is a modified version of ElementAverageValue postprocessor from MOOSE
// This postprocessor calculates the volume average only for elements which have a non-zero value

#include "ElementNormalizedValue.h"

registerMooseObject("RhocpApp", ElementNormalizedValue);

InputParameters
ElementNormalizedValue::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  params.addClassDescription("Computes the volumetric average of a variable only for elements which have a non-zero value");
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

  if (abs(computeIntegral()) >= 1.e-15) {
    _volume += _current_elem_volume;
  }
}

Real
ElementNormalizedValue::getValue() const
{
  Real vol;
  if (_volume == 0.e0) {
    vol = 1.e6;
  }
  else {
    vol = _volume;
  }
  return _integral_value / vol;
}

void
ElementNormalizedValue::finalize()
{
  gatherSum(_volume);
  gatherSum(_integral_value);
}

void
ElementNormalizedValue::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  const auto & pps = static_cast<const ElementNormalizedValue &>(y);
  _volume += pps._volume;
}
