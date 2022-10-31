#ifndef COMPUTECPELASTICITYTENSOR_H
#define COMPUTECPELASTICITYTENSOR_H

#include "ComputeElasticityTensorBase.h"

class ComputeCPElasticityTensor;

template <>
InputParameters validParams<ComputeCPElasticityTensor>();

class ComputeCPElasticityTensor : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();

  ComputeCPElasticityTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  const MaterialProperty<RankFourTensor> & _Cel_cp;
};
#endif //COMPUTECPELASTICITYTENSOR_H
