#include "GrainAverageValue.h"
#include "EBSDReaderPointDataAux.h"

registerMooseObject("RhocpApp", GrainAverageValue);

template <>
InputParameters
validParams<GrainAverageValue>()
{
  InputParameters params = validParams<ElementIntegralVariablePostprocessor>();
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
  Point p = _current_elem->centroid();
  EBSDAccessFunctors::EBSDPointData data = _EBSDFileReader->getData(p);
  unsigned int current_gr_id = data._feature_id;

  if (current_gr_id == _grain_id) {
      ElementIntegralVariablePostprocessor::execute();

      _volume += _current_elem_volume;
  }
}

Real
GrainAverageValue::getValue()
{
  Real integral = ElementIntegralVariablePostprocessor::getValue();

  gatherSum(_volume);

  return integral / _volume;
}

void
GrainAverageValue::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  const GrainAverageValue & pps = static_cast<const GrainAverageValue &>(y);
  _volume += pps._volume;
}
