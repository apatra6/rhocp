#include "StateVariable.h"

registerMooseObject("RhocpApp", StateVariable);

template<>
InputParameters validParams<StateVariable>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("sdv_id", "State variable number");
  return params;
}

StateVariable::StateVariable(const InputParameters & parameters) :
    AuxKernel(parameters),
    _sdv_id(getParam<unsigned int>("sdv_id")),
    _state_var(getMaterialProperty<std::vector<Real> >("state_var"))
{
}

Real StateVariable::computeValue()
{
  return _state_var[_qp][_sdv_id - 1];
}
