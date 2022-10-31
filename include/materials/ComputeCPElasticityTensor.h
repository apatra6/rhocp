#ifndef COMPUTECPELASTICITYTENSOR_H
#define COMPUTECPELASTICITYTENSOR_H

#pragma once

#include "ComputeElasticityTensorBase.h"

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
