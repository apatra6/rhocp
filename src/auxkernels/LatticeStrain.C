#include "LatticeStrain.h"
#include "RankTwoTensor.h"

registerMooseObject("RhocpApp", LatticeStrain);

InputParameters
LatticeStrain::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("AxuKernel for outputing lattice strain for a hkl plane and g vector");
  params.addRequiredParam<Real>("h", "Miller index: h");
  params.addRequiredParam<Real>("k", "Miller index: k");
  params.addRequiredParam<Real>("l", "Miller index: l");
  params.addRequiredParam<Real>("g0", "Diffraction vector: component 1");
  params.addRequiredParam<Real>("g1", "Diffraction vector: component 2");
  params.addRequiredParam<Real>("g2", "Diffraction vector: component 3");
  params.addParam<Real>("tol", 6.5 , "Tolerance angle");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", {}, "List of eigenstrains to be applied in this strain calculation");
  return params;
}

LatticeStrain::LatticeStrain(const InputParameters & parameters) :
    AuxKernel(parameters),
    _h(getParam<Real>("h")),
    _k(getParam<Real>("k")),
    _l(getParam<Real>("l")),
    _g0(getParam<Real>("g0")),
    _g1(getParam<Real>("g1")),
    _g2(getParam<Real>("g2")),
    _tol(getParam<Real>("tol")),
    _state_var(getMaterialProperty<std::vector<Real> >("state_var")),
    _euler_ang(getMaterialProperty<Point>("euler_ang")),
    _eigenstrain_names(getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
    _eigenstrains(_eigenstrain_names.size())
{
  for (auto i : make_range(_eigenstrain_names.size())) {
    _eigenstrains[i] = &getMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
  }
}

Real LatticeStrain::computeValue()
{

  RankTwoTensor R0, R, elR, elF, elE, Fp;

  // We can get the rotation tensor in two ways:

//   // Method 1: Creating the rotation tensor from the Euler angles
//   Point euler;
//   euler(0) = _euler_ang[_qp](0)*(22.0/7.0)/180.0;
//   euler(1) = _euler_ang[_qp](1)*(22.0/7.0)/180.0;
//   euler(2) = _euler_ang[_qp](2)*(22.0/7.0)/180.0;
//
//   Real s1 = sin(euler(0));
//   Real c1 = cos(euler(0));
//   Real s2 = sin(euler(1));
//   Real c2 = cos(euler(1));
//   Real s3 = sin(euler(2));
//   Real c3 = cos(euler(2));
//
//   // Bunge notation
//   R0(0,0) = c1*c3 - s1*s3*c2;
//   R0(1,0) = s1*c3 + c1*s3*c2;
//   R0(2,0) = s3*s2;
//   R0(0,1) = -c1*s3 - s1*c3*c2;
//   R0(1,1) = -s1*s3 + c1*c3*c2;
//   R0(2,1) = c3*s2;
//   R0(0,2) = s1*s2;
//   R0(1,2) = -c1*s2;
//   R0(2,2) = c2;

  // Method 2: getting the direction cosine tensor from DDCPStressUpdate
  R0(0,0) = _state_var[_qp][0];
  R0(0,1) = _state_var[_qp][1];
  R0(0,2) = _state_var[_qp][2];
  R0(1,0) = _state_var[_qp][3];
  R0(1,1) = _state_var[_qp][4];
  R0(1,2) = _state_var[_qp][5];
  R0(2,0) = _state_var[_qp][6];
  R0(2,1) = _state_var[_qp][7];
  R0(2,2) = _state_var[_qp][8];

  elF(0,0) = _state_var[_qp][18];
  elF(0,1) = _state_var[_qp][19];
  elF(0,2) = _state_var[_qp][20];
  elF(1,0) = _state_var[_qp][21];
  elF(1,1) = _state_var[_qp][22];
  elF(1,2) = _state_var[_qp][23];
  elF(2,0) = _state_var[_qp][24];
  elF(2,1) = _state_var[_qp][25];
  elF(2,2) = _state_var[_qp][26];

  elE = elF.transpose() * elF;
  elE(0,0) = elE(0,0) - 1.0;
  elE(1,1) = elE(1,1) - 1.0;
  elE(2,2) = elE(2,2) - 1.0;
  elE = 0.5 * elE;

  RankTwoTensor total_eigenstrain;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      total_eigenstrain(i,j) = 0.e0;
    }
  }

  for (auto i : make_range(_eigenstrain_names.size())) {
    total_eigenstrain += (*_eigenstrains[i])[_qp];
  }
  elE = elE + total_eigenstrain;

  std::vector<Real> e_value(3);
  RankTwoTensor e_vector, N1, N2, N3;

  RankTwoTensor temp1 = elF.transpose() * elF;
  temp1.symmetricEigenvaluesEigenvectors(e_value, e_vector);

  const Real lambda1 = std::sqrt(e_value[0]);
  const Real lambda2 = std::sqrt(e_value[1]);
  const Real lambda3 = std::sqrt(e_value[2]);

  N1.vectorOuterProduct(e_vector.column(0), e_vector.column(0));
  N2.vectorOuterProduct(e_vector.column(1), e_vector.column(1));
  N3.vectorOuterProduct(e_vector.column(2), e_vector.column(2));

  RankTwoTensor elU =  N1 * lambda1 + N2 * lambda2 + N3 * lambda3;
  RankTwoTensor invelU(elU.inverse());

  elR = elF * invelU;

  R =  elR * R0;

  Real N0[24][3], n[24][3], nr[24][3];

  N0[0][0] = _h; N0[0][1] = _k; N0[0][2] = _l;
  N0[1][0] = (-1.0)*_h; N0[1][1] = _k; N0[1][2] = _l;
  N0[2][0] = _h; N0[2][1] = (-1.0)*_k; N0[2][2] = _l;
  N0[3][0] = _h; N0[3][1] = _k; N0[3][2] = (-1.0)*_l;

  N0[4][0] = _h; N0[4][1] = _l; N0[4][2] = _k;
  N0[5][0] = (-1.0)*_h; N0[5][1] = _l; N0[5][2] = _k;
  N0[6][0] = _h; N0[6][1] = (-1.0)*_l; N0[6][2] = _k;
  N0[7][0] = _h; N0[7][1] = _l; N0[7][2] = (-1.0)*_k;

  N0[8][0] = _k; N0[8][1] = _h; N0[8][2] = _l;
  N0[9][0] = (-1.0)*_k; N0[9][1] = _h; N0[9][2] = _l;
  N0[10][0] = _k; N0[10][1] = (-1.0)*_h; N0[10][2] = _l;
  N0[11][0] = _k; N0[11][1] = _h; N0[11][2] = (-1.0)*_l;

  N0[12][0] = _k; N0[12][1] = _l; N0[12][2] = _h;
  N0[13][0] = (-1.0)*_k; N0[13][1] = _l; N0[13][2] = _h;
  N0[14][0] = _k; N0[14][1] = (-1.0)*_l; N0[14][2] = _h;
  N0[15][0] = _k; N0[15][1] = _l; N0[15][2] = (-1.0)*_h;

  N0[16][0] = _l; N0[16][1] = _h; N0[16][2] = _k;
  N0[17][0] = (-1.0)*_l; N0[17][1] = _h; N0[17][2] = _k;
  N0[18][0] = _l; N0[18][1] = (-1.0)*_h; N0[18][2] = _k;
  N0[19][0] = _l; N0[19][1] = _h; N0[19][2] = (-1.0)*_k;

  N0[20][0] = _l; N0[20][1] = _k; N0[20][2] = _h;
  N0[21][0] = (-1.0)*_l; N0[21][1] = _k; N0[21][2] = _h;
  N0[22][0] = _l; N0[22][1] = (-1.0)*_k; N0[22][2] = _h;
  N0[23][0] = _l; N0[23][1] = _k; N0[23][2] = (-1.0)*_h;

  Real val;

//---------------------------Removing repetition-------------------------

  for (unsigned int l=0; l<24; l++)
  {
    for (unsigned int k=0; k<24; k++)
    {
        if (l!=k){
          if (N0[l][0]==N0[k][0] && N0[l][1]==N0[k][1] && N0[l][2]==N0[k][2])
          {
            N0[k][0] = 0.0; N0[k][1] = 0.0; N0[k][2] = 0.0;
          }
        }
    }
  }

//--------------------------Normalizing vectors----------------------------

  for (unsigned int l=0; l<24; l++)
  {
    val = N0[l][0]*N0[l][0] + N0[l][1]*N0[l][1] + N0[l][2]*N0[l][2];
    if (val < 1.e-15)
    {
      N0[l][0] = 0.0;
      N0[l][1] = 0.0;
      N0[l][2] = 0.0;
    }
    else
    {
      N0[l][0] = N0[l][0]/sqrt(val);
      N0[l][1] = N0[l][1]/sqrt(val);
      N0[l][2] = N0[l][2]/sqrt(val);
    }
  }


  for (unsigned int l = 0; l < 24; l++)
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      n[l][i] = 0.0;
      for (unsigned int j = 0; j < 3; j++)
      {
          n[l][i] = n[l][i] + R(i,j) * N0[l][j];
      }
    }
  }


  Real latticeStrain, cosangle;
  latticeStrain = 0.0e0;

  for (unsigned int l = 0; l < 24; l++)
  {
    if ((n[l][0]*n[l][0] + n[l][1]*n[l][1] + n[l][2]*n[l][2])>1e-15)
    {
      cosangle = ((n[l][0]*_g0 + n[l][1]*_g1 + n[l][2]*_g2) / sqrt(n[l][0]*n[l][0] + n[l][1]*n[l][1] + n[l][2]*n[l][2]) / sqrt(_g0*_g0 + _g1*_g1 + _g2*_g2));

      if (abs(acos(cosangle)* 180.0/(22.0/7.0)) <= _tol)
      {
        for (unsigned int i = 0; i < 3; i++)
        {
            for (unsigned int j = 0; j < 3; j++)
            {
              latticeStrain = latticeStrain + (n[l][i]*elE(i,j)*n[l][j]);
            }
        }
        break;
      }
    }
  }
  return latticeStrain;

}
