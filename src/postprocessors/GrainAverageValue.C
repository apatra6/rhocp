#include "GrainAverageValue.h"
#include "EBSDReaderPointDataAux.h"

registerMooseObject("RhocpApp", GrainAverageValue);

InputParameters
GrainAverageValue::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  params.addRequiredParam<unsigned int>("grain_id", "Grain number");
  params.addRequiredParam<UserObjectName>("EBSDFileReader", "Name of the EBSDReader UO");
  return params;
}

GrainAverageValue::GrainAverageValue(const InputParameters & parameters):
  ElementIntegralVariablePostprocessor(parameters),
  _grain_id(getParam<unsigned int>("grain_id")),
  _EBSDFileReader(&getUserObject<EBSDMeshReader>("EBSDFileReader")),
  _volume(0)
{
}

void
GrainAverageValue::initialize()
{
  ElementIntegralVariablePostprocessor::initialize();
  _volume = 0;
}

void
GrainAverageValue::execute()
{
  Point p = _current_elem->vertex_average();
  EBSDAccessFunctors::EBSDPointData data = _EBSDFileReader->getData(p);
  unsigned int current_gr_id = data._feature_id;

  if (current_gr_id == _grain_id) {
      ElementIntegralVariablePostprocessor::execute();

      _volume += _current_elem_volume;
  }
}

Real
GrainAverageValue::getValue() const
{
  return _integral_value / _volume;
}

void
GrainAverageValue::finalize()
{
  gatherSum(_integral_value);
  gatherSum(_volume);
}

void
GrainAverageValue::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  const GrainAverageValue & pps = static_cast<const GrainAverageValue &>(y);

//   const auto & pps = static_cast<const GrainAverageValue &>(y);

  _integral_value += pps._integral_value;
  _volume += pps._volume;
}
