#include "OutputEuler.h"

registerMooseObject("RhocpApp", OutputEuler);

InputParameters
OutputEuler::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<unsigned int>("angle_id", "Euler angle id");
  return params;
}

OutputEuler::OutputEuler(const InputParameters & parameters) :
    AuxKernel(parameters),
    _angle_id(getParam<unsigned int>("angle_id")),
    _euler_ang(getMaterialProperty<Point>("euler_ang"))
{
}

Real OutputEuler::computeValue()
{
  return _euler_ang[_qp](_angle_id);
}
