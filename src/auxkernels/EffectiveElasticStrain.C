#include "EffectiveElasticStrain.h"
#include "RankTwoTensor.h"

registerMooseObject("RhocpApp", EffectiveElasticStrain);

InputParameters
EffectiveElasticStrain::validParams()
{
  InputParameters params = AuxKernel::validParams();
  return params;
}

EffectiveElasticStrain::EffectiveElasticStrain(const InputParameters & parameters) :
    AuxKernel(parameters),
    _state_var(getMaterialProperty<std::vector<Real> >("state_var"))
{
}

Real EffectiveElasticStrain::computeValue()
{
  RankTwoTensor elF, elE;

  elF(0,0) = _state_var[_qp][3];
  elF(0,1) = _state_var[_qp][4];
  elF(0,2) = _state_var[_qp][5];
  elF(1,0) = _state_var[_qp][6];
  elF(1,1) = _state_var[_qp][7];
  elF(1,2) = _state_var[_qp][8];
  elF(2,0) = _state_var[_qp][9];
  elF(2,1) = _state_var[_qp][10];
  elF(2,2) = _state_var[_qp][11];

  elE = elF.transpose() * elF;
  elE(0,0) = elE(0,0) - 1.0;
  elE(1,1) = elE(1,1) - 1.0;
  elE(2,2) = elE(2,2) - 1.0;
  elE = 0.5 * elE;

  Real effElasticStrain;
  effElasticStrain = 0.0e0;

  effElasticStrain = (2./3./sqrt(2.))*sqrt((elE(0,0) - elE(1,1))*(elE(0,0) - elE(1,1)) + (elE(1,1) - elE(2,2))*(elE(1,1) - elE(2,2)) + (elE(2,2) - elE(0,0))*(elE(2,2) - elE(0,0)) + 6*elE(0,1)*elE(0,1) + 6*elE(1,2)*elE(1,2) + 6*elE(2,0)*elE(2,0));

  return effElasticStrain;
}
