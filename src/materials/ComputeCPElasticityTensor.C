#include "ComputeCPElasticityTensor.h"

registerMooseObject("RhocpApp", ComputeCPElasticityTensor);

InputParameters
ComputeCPElasticityTensor::validParams()
{
  InputParameters params = ComputeElasticityTensorBase::validParams();
  params.addClassDescription("Compute CP elasticity tensor for crystal plasticity.");

  return params;
}

ComputeCPElasticityTensor::ComputeCPElasticityTensor(const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _Cel_cp(getMaterialProperty<RankFourTensor>("Cel_cp"))
{
  if (!isParamValid("elasticity_tensor_prefactor"))
    issueGuarantee(_elasticity_tensor_name, Guarantee::CONSTANT_IN_TIME);
}

void
ComputeCPElasticityTensor::computeQpElasticityTensor()
{
  _elasticity_tensor[_qp] = _Cel_cp[_qp];
}
