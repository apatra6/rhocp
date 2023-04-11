// This is a local copy of moose/modules/phase_field/src/aux_kernels/EBSDReaderPointDataAux.C
// This local class is used by EBSDMeshReader.C in rho-CP

#include "EBSDMeshReaderPointDataAux.h"

registerMooseObject("RhocpApp", EBSDMeshReaderPointDataAux);

InputParameters
EBSDMeshReaderPointDataAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>("ebsd_reader", "The EBSDReader GeneralUserObject");
  MooseEnum field_types = EBSDAccessFunctors::getPointDataFieldType();
  params.addRequiredParam<MooseEnum>(
      "data_name", field_types, "The data to be extracted from the EBSD data by this AuxKernel");
  return params;
}

EBSDMeshReaderPointDataAux::EBSDMeshReaderPointDataAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _ebsd_reader(getUserObject<EBSDMeshReader>("ebsd_reader")),
    _data_name(getParam<MooseEnum>("data_name")),
    _val(_ebsd_reader.getPointDataAccessFunctor(_data_name))
{
}

void
EBSDMeshReaderPointDataAux::precalculateValue()
{
  // EBSD data is defined at element centroids, so this only makes
  // sense as an Element AuxKernel
  Point p = _current_elem->vertex_average();

  _value = (*_val)(_ebsd_reader.getData(p));
}

Real
EBSDMeshReaderPointDataAux::computeValue()
{
  return _value;
}
