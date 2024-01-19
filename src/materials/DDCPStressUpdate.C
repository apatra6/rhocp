#include "DDCPStressUpdate.h"
#include "MatrixTools.h"
#include "MooseRandom.h"
#include "MooseException.h"
#include <algorithm>
#include <dlfcn.h>
#include <fstream>
#define QUOTE(macro) stringifyName(macro)

registerMooseObject("RhocpApp", DDCPStressUpdate);

InputParameters
DDCPStressUpdate::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Stress calculation using the Crystal Plasticity material model for FCC and BCC crystals");
  params.addRequiredParam<FileName>("propsFile", "The file with the material parameters");
  params.addRequiredParam<FileName>("slipSysFile", "The file with the crystallography of slip systems");
  params.addRequiredParam<unsigned int>("num_props", "The number of material properties this UMAT is going to use");
  params.addRequiredParam<unsigned int>("num_slip_sys", "The number of slip systems");
  params.addRequiredParam<unsigned int>("num_state_vars", "The number of state variables this UMAT is going to use");
  params.addParam<Real>("tol", 1.0e-6, "Tolerance");
  params.addCoupledVar("temp", 300, "Temperature");
  params.addParam<UserObjectName>("EulerAngFileReader", "Name of the EulerAngleReader UO");
  params.addParam<UserObjectName>("EBSDFileReader", "Name of the EBSDReader UO");
  params.addParam<UserObjectName>("GrainAreaSize", "Name of the GrainAreaSize UO");
  params.addParam<int>("isEulerRadian", 0, "Are Euler angles specified in radians");
  params.addParam<int>("isEulerBunge", 0, "Are Euler angles specified in Bunge notation");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", {}, "List of eigenstrains to be applied in this strain calculation");
  return params;
}

DDCPStressUpdate::DDCPStressUpdate(const InputParameters & parameters) :
    ComputeStressBase(parameters),
    _propsFile(getParam<FileName>("propsFile")),
    _slipSysFile(getParam<FileName>("slipSysFile")),
    _num_props(getParam<unsigned int>("num_props")),
    _num_slip_sys(getParam<unsigned int>("num_slip_sys")),
    _num_state_vars(getParam<unsigned int>("num_state_vars")),
    _tol(getParam<Real>("tol")),
    _temp(coupledValue("temp")),
    _EulerAngFileReader(isParamValid("EulerAngFileReader")
                               ? &getUserObject<EulerAngleReader>("EulerAngFileReader")
                               : NULL),
    _EBSDFileReader(isParamValid("EBSDFileReader")
                               ? &getUserObject<EBSDMeshReader>("EBSDFileReader")
                               : NULL),
    _GrainAreaSize(isParamValid("GrainAreaSize")
                               ? &getUserObject<GrainAreaSize>("GrainAreaSize")
                               : NULL),
    _isEulerRadian(getParam<int>("isEulerRadian")),
    _isEulerBunge(getParam<int>("isEulerBunge")),
    _euler_ang(declareProperty<Point>("euler_ang")),
    _eigenstrain_names(getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
    _eigenstrains(_eigenstrain_names.size()),
    _eigenstrains_old(_eigenstrain_names.size()),
    _state_var(declareProperty<std::vector<Real> >("state_var")),
    _state_var_old(getMaterialPropertyOld<std::vector<Real> >("state_var")),
    _properties(declareProperty<std::vector<Real> >("properties")),
    _properties_old(getMaterialPropertyOld<std::vector<Real> >("properties")),
    _y(declareProperty<std::vector<std::vector<Real> > >("slip plane normals")),
    _y_old(getMaterialPropertyOld<std::vector<std::vector<Real> > >("slip plane normals")),
    _z(declareProperty<std::vector<std::vector<Real> > >("slip plane directions")),
    _z_old(getMaterialPropertyOld<std::vector<std::vector<Real> > >("slip plane directions")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _rotation_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "rotation_increment")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _Cel_cp(declareProperty<RankFourTensor>("Cel_cp"))
{
  for (auto i : make_range(_eigenstrain_names.size())) {
    _eigenstrains[i] = &getMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
    _eigenstrains_old[i] = &getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_names[i]);
  }
}

DDCPStressUpdate::~DDCPStressUpdate()
{
}

void DDCPStressUpdate::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();

  //Initialize state variable vector
  _state_var[_qp].resize(_num_state_vars);

  // Initialize _properties vector
  _properties[_qp].resize(_num_props);

  // Initialize _y and _z vectors
  _y[_qp].resize(3, std::vector<Real>(_num_slip_sys));
  _z[_qp].resize(3, std::vector<Real>(_num_slip_sys));

} // End initQpStatefulProperties()

void DDCPStressUpdate::computeQpStress()
{
  Real PI = 4.0*atan(1.0);

  Real psi[3]; // Euler angles
  RankTwoTensor E_el; // Elastic Green stran tensor
  RankTwoTensor E_tot; // Green stran tensor E_el + E_p
  RankTwoTensor F0; // F at beginning of sub increment
  RankTwoTensor F1; // F at end of sub increment
  RankTwoTensor F_el; //  Elastic part of F
  RankTwoTensor F_el_inv; // Inverse of elastic part of F
  RankTwoTensor F_p_inv_0; // Inverse of plastic part of F at beginning of step
  RankTwoTensor F_p_inv; // Inverse of plastic part of F at end of step
  RankTwoTensor E_p; // Plastic strain tensor
  Real E_eff;
  Real E_p_eff;
  Real E_p_eff_cum;
  // Real gamma_dot_g_avg0;
  Real rho_m_avg0, rho_m_avg;
  Real rho_i_avg0, rho_i_avg;

  Real dir_cos[3][3]; // Rotation matrix
  Real dir_cos0[3][3]; // Rotation matrix at the beginning of simulation
  Real C0[3][3][3][3]; // Elastic stiffness tensor
  Real C[3][3][3][3];
  Real C_avg[3][3][3][3];

  RankTwoTensor sig_avg;
  RankTwoTensor sig;
  RankTwoTensor Spk2;

  Real ddpdsig[3][3][3][3];
  Real ddsdde_4th[3][3][3][3];

  // slip system dependent arrays
  std::vector<std::vector<Real>> xs0(3, std::vector<Real>(_num_slip_sys)); // intermediate config slip directions in global coords
  std::vector<std::vector<Real>> xs(3, std::vector<Real>(_num_slip_sys)); // current config slip directions in global coords
  std::vector<std::vector<Real>> xm0(3, std::vector<Real>(_num_slip_sys)); // intermediate config plane normals in global coords
  std::vector<std::vector<Real>> xm(3, std::vector<Real>(_num_slip_sys)); // current config plane normals in global coords

  std::vector<Real> rho_m0(_num_slip_sys); // mobile dislocation density at beginning of step
  std::vector<Real> rho_m(_num_slip_sys); // mobile dislocation density at end of step
  std::vector<Real> rho_i0(_num_slip_sys); // immobile dislocation density at beginning of step
  std::vector<Real> rho_i(_num_slip_sys); // immobile dislocation density at end of step
  std::vector<Real> bstress0(_num_slip_sys); // backstress at beginning of step
  std::vector<Real> bstress(_num_slip_sys); // backstress at end of step
  std::vector<Real> tau(_num_slip_sys); // resolved shear stress
  std::vector<Real> tau_eff(_num_slip_sys); // effective resolved shear stress
  std::vector<Real> s_a(_num_slip_sys); // athermal slip resistance
  std::vector<Real> s_t(_num_slip_sys); // thermal slip resistance
  std::vector<Real> gamma_dot(_num_slip_sys); // shear strain rate
  std::vector<Real> gamma_dot_g(_num_slip_sys); // shear strain rate due to glide
  std::vector<Real> gamma_try(_num_slip_sys); // trial shear strain rate
  std::vector<Real> gamma_try_g(_num_slip_sys); // trial shear strain rate due to glide

  std::vector<Real> residual(_num_slip_sys); // residual during Newton-Raphson iteration

  std::vector<std::vector<Real>> dTaudGd(_num_slip_sys, std::vector<Real>(_num_slip_sys)); // d(tau)/d(gamma_dot)

  // interaction coefficient martix between different slip systems
  std::vector<std::vector<Real>> A(_num_slip_sys, std::vector<Real>(_num_slip_sys));
  std::vector<std::vector<Real>> H(_num_slip_sys, std::vector<Real>(_num_slip_sys));

  // Our model is based on V.R decomposition as opposed to MOOSE's R.U decomposition of the deformation gradient
  // F = V.R = R.U
  // V = R.U.R_T
  // Our model needs everything in the current configuration
  RankTwoTensor delV;
  delV = _rotation_increment[_qp] * _strain_increment[_qp] * _rotation_increment[_qp].transpose();

  // Rotate stress to current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress_old[_qp] * _rotation_increment[_qp].transpose();

  // Recover "old" state variables
  std::copy(_state_var_old[_qp].begin(), _state_var_old[_qp].end(), _state_var[_qp].begin());
  std::copy(_properties_old[_qp].begin(), _properties_old[_qp].end(), _properties[_qp].begin());

  // Recover "old" miller indices
  std::copy(_y_old[_qp].begin(), _y_old[_qp].end(), _y[_qp].begin());
  std::copy(_z_old[_qp].begin(), _z_old[_qp].end(), _z[_qp].begin());

  if (_t_step <= 1) {
  // read in Euler angles
  if (_EBSDFileReader){
      Point p = _current_elem->vertex_average();
      EBSDAccessFunctors::EBSDPointData data = _EBSDFileReader->getData(p);
      _grainid = data._feature_id;

      psi[0] = data._phi1;
      psi[1] = data._Phi;
      psi[2] = data._phi2;
  }
  else if(_EulerAngFileReader){
      _grainid = _current_elem->subdomain_id();

      psi[0] = _EulerAngFileReader->getData(_grainid,0);
      psi[1] = _EulerAngFileReader->getData(_grainid,1);
      psi[2] = _EulerAngFileReader->getData(_grainid,2);
  }
  else{
      mooseError("Euler angle data not found.");
  }

  // read grain size, available
  if (_GrainAreaSize){
    _grain_size = _GrainAreaSize->getGrainSize(_grainid);
  }
  else{
    _grain_size = 1.e10;
  }
  }

  // Read in material parameters
  if (_t_step <= 1){
      readPropsFile();
  }
  assignProperties();

  // Initialize Kronecker delta tensor
  Real del[3][3];
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      del[i][j] = 0;
    }
    del[i][i] = 1;
  }

  // Read in slip system data
  if (_t_step <= 1){
    MooseUtils::checkFileReadable(_slipSysFile);
    std::ifstream file_slip_sys;
    file_slip_sys.open(_slipSysFile.c_str());

    // Assign slip system normals and slip directions
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      for(unsigned int j = 0; j < 3; j++){
        //y[j][i] = reader.getData(c++)[0];
        file_slip_sys >> _y[_qp][j][i];
      }
      for (unsigned int j = 0; j < 3; j++){
        //z[j][i] = reader.getData(c++)[0];
        file_slip_sys >> _z[_qp][j][i];
      }
    }
    file_slip_sys.close();

    // Normalize the Miller indices to unit vectors
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      normalize_vector(&_y[_qp][0][i],&_y[_qp][1][i],&_y[_qp][2][i]);
      normalize_vector(&_z[_qp][0][i],&_z[_qp][1][i],&_z[_qp][2][i]);
    }

    for (unsigned int k = 0; k < _num_slip_sys; k++) {
      Real sum1 = 0.0;
      for (unsigned int i = 0; i<3; i++) {
        sum1 = sum1 + _y[_qp][i][k]*_z[_qp][i][k];
      }
      if (abs(sum1) > 1.0e-8) {
        mooseError("The Miller indices are WRONG for slip system:",(k + 1));
        break;
      }
    }
  }

  // Initialize interaction coefficient matrix
  for (unsigned int i = 0; i < _num_slip_sys; i++) {
    for (unsigned int j = 0; j < _num_slip_sys; j++) {
      A[i][j] = Alatent;
      H[i][j] = 0.0;
    }
    A[i][i] = 1.0;
    H[i][i] = 1.0;
  }

  if (_t_step == 1) {
    // Convert Euler angles to radians
    psi[0] = psi[0]*PI/180;
    psi[1] = psi[1]*PI/180;
    psi[2] = psi[2]*PI/180;

    // Initialize F_p_inv_0
    for (unsigned int i = 0; i<3; i++) {
      for(unsigned int j = 0; j<3; j++) {
        F_p_inv_0(i,j) = 0.0;
      }
      F_p_inv_0(i,i) = 1.0;
    }

    // Initialize E_p.
    for (unsigned int i = 0; i < 3; i++) {
      for(unsigned int j = 0; j < 3; j++) {
        E_p(i,j) = 0;
      }
    }

    // Initialize E_eff
    E_eff = 0;

    // Initialize E_p_eff
    E_p_eff = 0;

    // Initialize E_p_eff_cum
    E_p_eff_cum = 0;

    // Initialize average values of defect densities and the corresponding loop sizes
    rho_m_avg0 = rho_m_zero;
    rho_i_avg0 = rho_i_zero;

    for (unsigned int n = 0; n<_num_slip_sys; n++){
      rho_m0[n] = rho_m_zero;
    }

    for (unsigned int n = 0; n<_num_slip_sys; n++){
      rho_i0[n] = rho_i_zero;
    }

    for (unsigned int n = 0; n<_num_slip_sys; n++){
      bstress0[n] = 0.0;
    }
  }
  //  End of initializations.  Read in internal variables

  if(_t_step > 1) {
    // checkpoint("Initializing ISVs, nonzero time step")
    int n = -1;
    // Read in dir_cos0 SDV 1-9
    for (unsigned int i = 0; i<3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        dir_cos0[i][j] = _state_var[_qp][n];
      }
    }

    // Read in dir_cos SDV 10-18
    for (unsigned int i = 0; i<3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        dir_cos[i][j] = _state_var[_qp][n];
      }
    }

    // Read the elastic part of F SDV 19-27
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        // F_el[i][j] = _state_var[_qp][n]; // Not needed. This is only for post-processing
      }
    }

    // Read inverse of the plastic part of F SDV 28-36
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        F_p_inv_0(i,j) = _state_var[_qp][n];
      }
    }

    // Read E_p SDV 37-45
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        E_p(i,j) = _state_var[_qp][n];
      }
    }

    // Read E_eff SDV 46
    n = n + 1;
    E_eff = _state_var[_qp][n];

    // Read E_p_eff SDV 47
    n = n + 1;
    E_p_eff = _state_var[_qp][n];

    // Read average values of defect densities and loop sizes SDV 48-49
    n = n + 1;
    rho_m_avg0 = _state_var[_qp][n];

    n = n + 1;
    rho_i_avg0 = _state_var[_qp][n];

    // Read average dislocation glide rates  SDV 50
    n = n + 1;
    // gamma_dot_g_avg0 = _state_var[_qp][n];

    // Read mobile dislocation density values SDV 50 + _num_slip_sys
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      rho_m0[i] = _state_var[_qp][n];
    }

    // Read immobile dislocation density values SDV 50 + _num_slip_sys*2
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      rho_i0[i] = _state_var[_qp][n];
    }

    // Read back stress values SDV 50 + _num_slip_sys*3
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      bstress0[i] = _state_var[_qp][n];
    }

    // std::cout << "\nTotal no. of state vars read:" << n;
  } // End of initializations

  if(_t_step == 1) {
    Real s1 = sin(psi[0]);
    Real c1 = cos(psi[0]);
    Real s2 = sin(psi[1]);
    Real c2 = cos(psi[1]);
    Real s3 = sin(psi[2]);
    Real c3 = cos(psi[2]);

    dir_cos0[0][0] = c1*c3-s1*s3*c2;
    dir_cos0[0][1] = s1*c3+c1*s3*c2;
    dir_cos0[0][2] = s3*s2;
    dir_cos0[1][0] = -c1*s3-s1*c3*c2;
    dir_cos0[1][1] = -s1*s3+c1*c3*c2;
    dir_cos0[1][2] = c3*s2;
    dir_cos0[2][0] = s1*s2;
    dir_cos0[2][1] = -c1*s2;
    dir_cos0[2][2] = c2;

    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        dir_cos[i][j] = dir_cos0[i][j];
      }
    }
  }

  // Initialize ANISOTROPIC 4th rank elastic stiffness tensor
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          C0[i][j][k][l] = C12*del[i][j]*del[k][l] + C44*(del[i][k]*del[j][l]+del[i][l]*del[k][j]);
        }
      }
    }
  }
  C0[0][0][0][0] = C11;
  C0[1][1][1][1] = C11;
  C0[2][2][2][2] = C11;

  // Initialize arrays for averaging over grains
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig_avg(i,j) = 0.0;
    }
  }

  // Begin calculations over the element
  rotate_4th(dir_cos, C0, C);

  // Convert Miller Indices to global coordinates
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    for (unsigned int i = 0; i < 3; i++) {
      xs0[i][n] = 0.0;
      xm0[i][n] = 0.0;
      for (unsigned int j = 0; j < 3; j++) {
        xs0[i][n] = xs0[i][n] + dir_cos0[i][j]*_z[_qp][j][n];
        xm0[i][n] = xm0[i][n] + dir_cos0[i][j]*_y[_qp][j][n];
      }
    }
  }

  // Initialize number of sub-increments.  Note that the subincrement initially equals the total increment. This remains the case unless the process starts to diverge.
  unsigned int N_incr, N_incr_total, iNR, N_ctr;
  Real dt_incr;
  RankTwoTensor dfgrd0, dfgrd1, Ftheta, Ftheta_old;
  RankTwoTensor total_eigenstrain, total_eigenstrain_old;
  RankTwoTensor array1, array2;
  Real array1_matrix[3][3], array2_matrix[3][3];

  bool converged, improved;

  // Initialize deformation gradients for beginning and end of subincrement.
  for (auto i : make_range(_eigenstrain_names.size())) {
    total_eigenstrain += (*_eigenstrains[i])[_qp];
    total_eigenstrain_old += (*_eigenstrains_old[i])[_qp];
  }

  Ftheta.zero();
  Ftheta_old.zero();
  for (int i=0; i<3; ++i)
  {
      Ftheta(i,i) = sqrt(1.0 + 2.0*total_eigenstrain(i,i));
      Ftheta_old(i,i) = sqrt(1.0 + 2.0*total_eigenstrain_old(i,i));
  }

  dfgrd0 = _deformation_gradient_old[_qp] * Ftheta_old.inverse();
  dfgrd1 = _deformation_gradient[_qp] * Ftheta.inverse();
  F0 = dfgrd0;
  F1 = dfgrd1;

  // counters for time step sub-increment
  N_incr = 1;
  N_incr_total = 1;

  Real sse_old, sse_ref;

  // Begin time step subincrement
  do {
  dt_incr = _dt/N_incr_total;

  F0 = dfgrd0 + (dfgrd1 - dfgrd0)*(N_incr - 1)/N_incr_total;
  F1 = dfgrd0 + (dfgrd1 - dfgrd0)*N_incr/N_incr_total;

  // Initialize mobile dislocation density
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    rho_m[n] = rho_m0[n];
  }

  // Initialize immobile dislocation density
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    rho_i[n] = rho_i0[n];
  }

  // Initialize backstress
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    bstress[n] = bstress0[n];
  }

  // Initialize average defect densities and loop sizes
  rho_m_avg = rho_m_avg0;
  rho_i_avg = rho_i_avg0;


  // Multiply F() by F_p_inv() to get F_el()
  F_el = F0 * F_p_inv_0;
  if (F_el.det() == 0.0e0) {
    F_el(0,0) = 1.0e0;
    F_el(1,1) = 1.0e0;
    F_el(2,2) = 1.0e0;
  }

  F_el_inv = F_el.inverse();

  // Rotate xs0 and xm0 to current coordinates, called xs and xm
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    for (unsigned int i = 0; i < 3; i++) {
      xs[i][n] = 0.0;
      xm[i][n] = 0.0;
      for (unsigned int j = 0; j < 3; j++) {
        xs[i][n] = xs[i][n] + F_el(i,j)*xs0[j][n];
        xm[i][n] = xm[i][n] + xm0[j][n]*F_el_inv(j,i);
      }
    }
  }

  // Calculate elastic Green strain
  E_el = F_el.transpose() * F_el;
  E_el(0,0) = E_el(0,0) - 1.0;
  E_el(1,1) = E_el(1,1) - 1.0;
  E_el(2,2) = E_el(2,2) - 1.0;
  E_el = 0.5 * E_el;

  // Multiply the anisotropic stiffness tensor by the Green strain to get the 2nd Piola Kirkhhoff stress
  Real E_el_matrix[3][3], Spk2_matrix[3][3];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      E_el_matrix[i][j] = E_el(i,j);
    }
  }
  aaaa_dot_dot_bb(C,E_el_matrix,Spk2_matrix);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Spk2(i,j) = Spk2_matrix[i][j];
    }
  }

  // Convert from PK2 stress to Cauchy stress
  Real det = F_el.det();
  sig =  F_el * Spk2 * F_el.transpose();
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig(i,j) = sig(i,j)/det;
    }
  }

  // Calculate resolved shear stress for each slip system
  for (unsigned int k = 0; k < _num_slip_sys; k++) {
    tau[k] = 0.0;
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        tau[k] = tau[k] + xs[i][k]*xm[j][k]*sig(i,j);
      }
    }
  }

  // Calculate effective shear stress for each slip system
  for (unsigned int k = 0; k < _num_slip_sys; k++) {
    tau_eff[k] = tau[k] - bstress[k];
  }

  // Calculate athermal slip resistance for each slip system.
  for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    Real sum1 = 0.0;
    for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {
      sum1 = sum1 + p0*A[ia][ib]*(rho_m[ib] + rho_i[ib]);
    }
    s_a[ia] = tau0 + hp_coeff/sqrt(grain_size) + k_rho*G*b_mag*sqrt(sum1);
  }

  // Calculate reference shear stress for each slip system.
  for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    s_t[ia] = frictional_stress;
  }

  // Calculate 1st estimate of gamma_dot for each slip system.
  std::vector<Real> d_disl(_num_slip_sys); // dislocation spacing

  for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    Real sum1 = 0.0;
    for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {
      sum1 = sum1 + H[ia][ib]*(rho_m[ib] + rho_i[ib]);
    }
    d_disl[ia] = 1.e0/sqrt(sum1);
  }

  for(unsigned int k = 0; k < _num_slip_sys; k++) {
    Real refgdg = gammadot0g;
    if(abs(tau_eff[k]) > s_a[k]) {
      gamma_dot_g[k] = refgdg*exp((-delF0/B_k/_temp[_qp])*power(1 - power((abs(tau_eff[k]) - s_a[k])/s_t[k], p), q))*sgn(tau_eff[k]);
    }
    else {
      gamma_dot_g[k] = 1e-18;
    }
    gamma_dot[k] = gamma_dot_g[k];
  }

  // Calculate d(Tau)/d(Gamma_dot)
  for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
      dTaudGd[ia][ib] = 0.0;

      for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 3; j++) {
          array1_matrix[i][j] = xs0[i][ib]*xm0[j][ib];
        }
      }
      aaaa_dot_dot_bb(C0,array1_matrix,array2_matrix);

      for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 3; j++) {
          array1_matrix[i][j] = xs0[i][ia]*xm0[j][ia];
        }
      }
      dTaudGd[ia][ib] = aa_dot_dot_bb(array1_matrix,array2_matrix);
      dTaudGd[ia][ib] = -1.0*dTaudGd[ia][ib]*dt_incr;
    }
  }

  // Begin Newton-Raphson iterative loops.
  iNR = 0;

  converged = false;

  do {
    iNR = iNR + 1;
    // converged = true;

    NR_residual (_num_slip_sys, xs0, xm0, _temp[_qp], dt_incr, gamma_dot, F1, F_el, F_p_inv, F_p_inv_0, C, rho_m0, rho_m, rho_i0, rho_i, bstress0, bstress, sig, tau, tau_eff, s_a, s_t, A, H, residual, sse);

    sse_ref = sse;

//       if (sse > 0.0e0) {
//         std::cout << "\n iNR:" << iNR;
//         std::cout << "\n sse:" << sse << "\n";
//       }

    // Begin calculation of the partial derivatives needed for the Newton-Raphson step

    // Calculate derivative of the mobile dislocation density, rho_m, immobile dislocation density, rho_i, backstress, bstress, w.r.t. gamma-dot-beta
    std::vector<std::vector<Real>> temp1(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> temp2(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> drhomdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> drhoidgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dbsdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dldrhom(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dsadgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dstdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        temp1[ia][ib] = 0.0;
        temp2[ia][ib] = 0.0;
        drhomdgb[ia][ib] = 0.0;
        drhoidgb[ia][ib] = 0.0;
        dbsdgb[ia][ib] = 0.0;
        dldrhom[ia][ib] = 0.0;
      }
    }

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      Real sum1 = 0.0;
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        sum1 = sum1 + H[ia][ib]*(rho_m[ib] + rho_i[ib]);
      }
      dldrhom[ia][ia] = - 0.5*pow(d_disl[ia],2.0)*H[ia][ia]/sqrt(sum1);

      temp2[ia][ia] = temp2[ia][ia] + (k_M/b_mag/d_disl[ia] - k_I/b_mag/d_disl[ia]
                        - (k_ann*R_c*rho_m[ia]/b_mag))*(sgn(gamma_dot[ia])*dt_incr);

      temp1[ia][ia] = temp1[ia][ia] + 1.0 + (k_M/b_mag/pow(d_disl[ia], 2))*dldrhom[ia][ia]*abs(gamma_dot[ia])
                        *dt_incr + k_ann*R_c*abs(gamma_dot[ia])*dt_incr/b_mag
                        + 0.5*k_I*H[ia][ia]*abs(gamma_dot[ia])*dt_incr/sqrt(sum1)/b_mag;
    }

    std::vector<std::vector<Real>> temp1_inv(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    // Method 2: MatrixTools inverse
    MatrixTools::inverse (temp1,temp1_inv);

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        drhomdgb[ia][ib] = 0.0;
        for(unsigned int ic = 0; ic < _num_slip_sys; ic++) {
          drhomdgb[ia][ib] = drhomdgb[ia][ib] + temp1_inv[ia][ic]*temp2[ic][ib];
        }
      }
    }

    // immobile dislocations
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      Real sum1 = 0.0;
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        drhoidgb[ia][ib] = 0.0;
        sum1 = sum1 + H[ia][ib]*(rho_m[ib] + rho_i[ib]);
      }
      Real t3 = k_I*sqrt(sum1)/b_mag - k_D*rho_i[ia];
      Real t4 = 1.0e0 - dt_incr*abs(gamma_dot[ia])*0.5*k_I*H[ia][ia]/b_mag/sqrt(sum1) + k_D*dt_incr*abs(gamma_dot[ia]);

      drhoidgb[ia][ia] = t3*sgn(gamma_dot[ia])*dt_incr/t4;
    }

    // backstress
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        dbsdgb[ia][ib] = dbsdgb[ia][ib] + k_bs1*(0.5*(k_rho*G*b_mag)/sqrt(rho_m[ia] + rho_i[ia]))*dt_incr*(drhomdgb[ia][ib] + drhoidgb[ia][ib])*sgn(tau[ia] -bstress[ia])*abs(gamma_dot[ia]);
      }

      dbsdgb[ia][ia] = dbsdgb[ia][ia] + k_bs1*k_rho*G*b_mag*dt_incr*sqrt(rho_m[ia] + rho_i[ia])*sgn(tau[ia] - bstress[ia])*sgn(gamma_dot[ia]) - k_bs2*dt_incr*  bstress[ia]*sgn(gamma_dot[ia]);
    }

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        dbsdgb[ia][ib] = dbsdgb[ia][ib]/(1.e0 + k_bs2*abs(gamma_dot[ia])*dt_incr);
      }
    }


    // Calculate derivative of the thermal slip resistance, s_t w.r.t. gamma-dot-beta
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        dstdgb[ia][ib] = 0.0;
      }
    }

    // Calculate derivative of the athermal slip resistance, s_a w.r.t. gamma-dot-beta.
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      Real sum1 = 0.0;
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        sum1 = sum1 + p0*A[ia][ib]*(rho_m[ib] + rho_i[ib]);

        dsadgb[ia][ib] = 0.5*k_rho*G*b_mag*(p0*A[ia][ib]*(drhomdgb[ia][ib]
          + drhoidgb[ia][ib])/sqrt(sum1));
      }
    }

    // Form "A-matrix" of derivatives wrt d_gamma_beta.
    std::vector<std::vector<Real>> array3(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      tau_eff[ia] = tau[ia] - bstress[ia];

      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        array3[ia][ib] = -gamma_dot_g[ia]*(p*q*(delF0/B_k/_temp[_qp])*power(1.e0 - power((abs(tau_eff[ia])
            - s_a[ia])/s_t[ia], p), q - 1)*power((abs(tau_eff[ia]) - s_a[ia])/s_t[ia], p - 1))*
          (((dTaudGd[ia][ib] - dbsdgb[ia][ib])*sgn(tau_eff[ia]) - dsadgb[ia][ib])/s_t[ia] - (abs(tau_eff[ia]) - s_a[ia])*dstdgb[ia][ib]/(s_t[ia]*s_t[ia]));
      }
      array3[ia][ia] = array3[ia][ia] + 1;
    }

    // Calculate the gradient of sse wrt gamma_dot(). Will be used later to ensure that line search is in the correct direction
    std::vector<Real> gradient(_num_slip_sys);

    for(unsigned int j = 0; j < _num_slip_sys; j++) {
      gradient[j] = 0.0;
      for(unsigned int i = 0; i < _num_slip_sys; i++) {
        gradient[j] = gradient[j] + residual[i]*array3[i][j];
      }
      gradient[j] = 2.0*gradient[j];
    }

    // Solve for increment of gamma_dot
    std::vector<Real> d_gamma_dot(_num_slip_sys);

    std::vector<std::vector<Real>> array3_inv(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    MatrixTools::inverse(array3,array3_inv);

    for(unsigned int j = 0; j < _num_slip_sys; j++) {
      d_gamma_dot[j] = 0.0;
      for(unsigned int i = 0; i < _num_slip_sys; i++) {
        d_gamma_dot[j] = d_gamma_dot[j] - array3_inv[j][i]*residual[i];
      }
    }

    // Check to make sure that N-R step leads 'down hill' the sse surface
    Real sum1 = 0.0;
    for(unsigned int i = 0; i < _num_slip_sys; i++) {
      sum1 = sum1 - gradient[i]*d_gamma_dot[i];
    }

    if(sum1 > 0.0) {
      for(unsigned int i = 0; i < _num_slip_sys; i++) {
        d_gamma_dot[i] = -1*d_gamma_dot[i];
      }
    }

    // Multiply step size by two 'cause next loop will divide it by 2
    for(unsigned int k = 0; k < _num_slip_sys; k++) {
      d_gamma_dot[k] = d_gamma_dot[k]*2;
    }

    // Begin line search
    improved = false;

    N_ctr = 0;

    do {
      sse_old = sse;

      // Divide step size by 2
      for (unsigned int k = 0; k < _num_slip_sys; k++) {
        d_gamma_dot[k] = d_gamma_dot[k]/2.0;

        gamma_try_g[k] = gamma_dot_g[k] + d_gamma_dot[k];
        gamma_try[k] = gamma_try_g[k];
      }

      NR_residual (_num_slip_sys, xs0, xm0, _temp[_qp], dt_incr, gamma_try, F1, F_el, F_p_inv, F_p_inv_0, C, rho_m0, rho_m, rho_i0, rho_i, bstress0, bstress, sig, tau, tau_eff, s_a, s_t, A, H, residual, sse);

      if ((sse_old <= sse_ref) && (sse >= sse_old) && (N_ctr > 0)) {
        improved = true;

        break;
      }

      N_ctr++;

    } while ((improved == false) && (N_ctr < max_loops)); // End line search

    // add d_gamma_dot to gamma_dot to get new values for this iteration
    for (unsigned int k = 0; k < _num_slip_sys; k++) {
      gamma_dot[k] = gamma_dot[k] + d_gamma_dot[k]*2.0;
      gamma_dot_g[k] = gamma_dot_g[k] + d_gamma_dot[k]*2.0;
    }

    // if (sse_old > tolerance) then this step has not converged
    if (sse_old <= _tol) {
      converged = true;

      // break;
    }

    // if (sse_old > sse_ref/2) then convergence is too slow and increment is divided into two subincrements
    if ((sse_old > (sse_ref/2.0)) && (converged == false)) {
      N_incr = 2*N_incr - 1;
      N_incr_total = 2*N_incr_total;

      // std::cout << "\n Total sub-increments:" << N_incr_total << "\n";
      break;
    }

  } while((iNR < 10000) && (converged == false));

  if (iNR == 10000) {
    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        _stress[_qp](i,j) = MooseRandom::rand()*2.0e3;
        for (unsigned int k = 0; k < 3; k++){
          for (unsigned int l = 0; l < 3; l++){
            _Jacobian_mult[_qp](i,j,k,l) = MooseRandom::rand()*2.0e3;
          }
        }
      }
    }

    throw MooseException("Crystal plasticity Newton Raphson did not converge");
    return;
  }

  // if another subincrement remains to be done, then reinitialize F_p_inv_0 and state variables.  F0 and F1 gets reinitialized back at the top of this loop
  if (N_incr < N_incr_total) {
    if (N_incr == (N_incr/2)*2) { // even
      N_incr = N_incr/2 + 1;
      N_incr_total = N_incr_total/2;
    }
    else { // odd
      N_incr = N_incr + 1;
    }

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        F_p_inv_0(i,j) = F_p_inv(i,j);
      }
    }

    for (unsigned int k = 0; k < _num_slip_sys; k++) {
      rho_m0[k] = rho_m[k];
      rho_i0[k] = rho_i[k];
      bstress0[k] = bstress[k];
    }
  }
  else if (N_incr == N_incr_total) {
    break;
  }
  } while (N_incr_total < 10000); // End time step subincrementation

  if (N_incr_total > 10000) {
    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        _stress[_qp](i,j) = MooseRandom::rand()*2.0e3;
        for (unsigned int k = 0; k < 3; k++){
          for (unsigned int l = 0; l < 3; l++){
            _Jacobian_mult[_qp](i,j,k,l) = MooseRandom::rand()*2.0e3;
          }
        }
      }
    }

    throw MooseException("Crystal plasticity time step subincrement did not converge");
    return;
  }

  // Calculate the average Cauchy stress
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig_avg(i,j) = sig_avg(i,j) + sig(i,j);
    }
  }

  Real F_el_matrix[3][3];
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++) {
      F_el_matrix[i][j] = F_el(i,j);
    }
  }

  // Write out Euler Angles in Bunge Convention
  aa_dot_bb(F_el_matrix,dir_cos0,dir_cos);
  bunge_angles(dir_cos,psi);

  // Calculate average elasticity tensor
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          C_avg[i][j][k][l] = C[i][j][k][l];
        }
      }
    }
  }

  // Rotate xs0 and xm0 to current coordinates, called xs and xm
  F_el_inv = F_el.inverse();
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    for (unsigned int i = 0; i < 3; i++) {
      xs[i][n] = 0.0;
      xm[i][n] = 0.0;
      for (unsigned int j = 0; j < 3; j++) {
        xs[i][n] = xs[i][n] + F_el(i,j)*xs0[j][n];
        xm[i][n] = xm[i][n] + xm0[j][n]*F_el_inv(j,i);
      }
    }
  }

  // Calculate the derivative of the plastic part of the rate of deformation tensor in the current configuration wrt sigma
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          ddpdsig[i][j][k][l] = 0.0e0;
        }
      }
    }
  }

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          for(unsigned int n = 0; n < _num_slip_sys; n++) {
            tau_eff[n] = tau[n] - bstress[n];

            ddpdsig[i][j][k][l] = ddpdsig[i][j][k][l] + (xs[i][n]*xm[j][n] + xm[i][n]*xs[j][n])*(((xs[k][n]*xm[l][n] + xm[k][n]*xs[l][n])*(gamma_dot_g[n]*(p*q*(delF0/B_k/_temp[_qp])*power(1.0 - power((abs(tau_eff[n]) - s_a[n])/s_t[n], p), q - 1.0)*power((abs(tau_eff[n]) - s_a[n])/s_t[n], p - 1.0))*sgn(tau_eff[n])/s_t[n])));
          }
        }
      }
    }
  }

  // End calculations over the integration point

  // Calculate Green strain
  E_tot = F1.transpose() * F1;
  E_tot(0,0) = E_tot(0,0) - 1.0;
  E_tot(1,1) = E_tot(1,1) - 1.0;
  E_tot(2,2) = E_tot(2,2) - 1.0;
  E_tot = 0.5 * E_tot;

  // Begin calculation of the Jacobian (the tangent stiffness matrix)

  // Store sig_avg() into sig()
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig(i,j) = sig_avg(i,j);
    }
  }

  // Output stress, strain, etc. for post processing
  Real sum1;
  sum1 = 0.0;
  for(unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sum1 = sum1 + sig(i,j) * sig(i,j);
    }
  }

  // Calculate the inverse of F_el
  F_el_inv = F_el.inverse();

  // Scale by appropriate constants and divide by num_grains to get the average value. ALSO multiply by 'dtime' which is d(sig)/d(sig_dot)
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          ddpdsig[i][j][k][l] = ddpdsig[i][j][k][l]*_dt/4.0;
        }
      }
    }
  }

  // Multiply the 4th rank elastic stiffness tensor by the derivative of the plastic part of the rate of deformation tensor wrt sig_dot
  Real array6[3][3][3][3], array4[6][6], array4_inv[6][6];
  aaaa_dot_dot_bbbb(C_avg,ddpdsig,array6);

  std::vector<std::vector<Real>> array4_matrix(6, std::vector<Real>(6));
  std::vector<std::vector<Real>> array4_inv_matrix(6, std::vector<Real>(6));

  // Add 4th rank identity tensor to array6()
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          array6[i][j][k][l] = array6[i][j][k][l] + 0.5*(del[i][k]*del[j][l] + del[i][l]*del[j][k]);
        }
      }
    }
  }

  forth_to_Voigt(array6,array4);

  for (unsigned int i = 0; i < 6; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      array4_matrix[i][j] = array4[i][j];
    }
  }

  // Method 2: MatrixTools::inverse
  MatrixTools::inverse(array4_matrix,array4_inv_matrix);

  for (unsigned int i = 0; i < 6; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      array4_inv[i][j] = array4_inv_matrix[i][j];
    }
  }

  Voigt_to_forth(array4_inv,array6);

//   // Method 3: inverse of 4th order symmetric tensor
//   // Need to take the inverse of Array4.  Since it relates two 2nd rank tensors that are both symmetric, Array4 can be xformed to Voigt notation style in order to do the inverse, then xformed back.
//   RankFourTensor array6RFT, array6RFTinv;
//
//   for (unsigned int i = 0; i < 3; i++) {
//       for (unsigned int j = 0; j < 3; j++) {
//           for (unsigned int k = 0; k < 3; k++) {
//               for (unsigned int l = 0; l < 3; l++) {
//                   array6RFT(i,j,k,l) = array6[i][j][k][l];
//               }
//           }
//       }
//   }
//
//   array6RFTinv = array6RFT.invSymm();
//
//   for (unsigned int i = 0; i < 3; i++) {
//       for (unsigned int j = 0; j < 3; j++) {
//           for (unsigned int k = 0; k < 3; k++) {
//               for (unsigned int l = 0; l < 3; l++) {
//                   array6[i][j][k][l] = array6RFTinv(i,j,k,l);
//               }
//           }
//       }
//   }

  //  Multiply Array6 by C, the elastic stiffness matrix to finally get the Jacobian, but as a 4th rank tensor.
  aaaa_dot_dot_bbbb (array6,C_avg,ddsdde_4th);

  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++){
          _Jacobian_mult[_qp](i,j,k,l) = ddsdde_4th[i][j][k][l];
        }
      }
    }
  }

  _stress[_qp] = sig;

  _euler_ang[_qp](0) = psi[0];
  _euler_ang[_qp](1) = psi[1];
  _euler_ang[_qp](2) = psi[2];

  // checkpoint("Storing ISVs")
  int n = -1;
  // Store dir_cos0 SDV 1-9
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = dir_cos0[i][j];
    }
  }

  // Store dir_cos SDV 10-18
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = dir_cos[i][j];
    }
  }

  // Store the elastic part of F SDV 19-27
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = F_el(i,j);
    }
  }

  // Store inverse of the plastic part of F SDV 28-36
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = F_p_inv(i,j);
    }
  }

  // Plastic strain calculations
  RankTwoTensor F_p;
  F_p = F_p_inv.inverse();

  E_p = F_p.transpose() * F_p;
  E_p(0,0) = E_p(0,0) - 1.0;
  E_p(1,1) = E_p(1,1) - 1.0;
  E_p(2,2) = E_p(2,2) - 1.0;
  E_p = 0.5 * E_p;

  // Store E_p SDV 37-45
  for (unsigned int i = 0; i<3; i++) {
    for (unsigned int j = 0; j<3; j++) {
      n = n + 1;
      _state_var[_qp][n] = E_p(i,j);
    }
  }

  // Effective strain calculations
  Real s1 = 0.0; Real s2 = 0.0;

  for (unsigned int i = 0; i<3; i++) {
    for (unsigned int j = 0; j<3; j++) {
      s1 = s1 + E_tot(i,j)*E_tot(i,j);
      s2 = s2 + E_p(i,j)*E_p(i,j);
    }
  }

  // Effective total strain
  E_eff = sqrt((2./3.) * s1);

  // Effective plastic strain
  E_p_eff = sqrt((2./3.) * s2);

  // Store E_eff SDV 46
  n = n + 1;
  _state_var[_qp][n] = E_eff;

  // Store E_p_eff SDV 47
  n = n + 1;
  _state_var[_qp][n] = E_p_eff;

  // Store average values of defect densities and loop sizes SDV 48-49
  rho_m_avg = 0.0;
  rho_i_avg = 0.0;
  Real gamma_dot_avg = 0.0;
  for (unsigned int i = 0; i < _num_slip_sys; i++) {
    rho_m_avg = rho_m_avg + rho_m[i];
    rho_i_avg = rho_i_avg + rho_i[i];
    gamma_dot_avg = gamma_dot_avg + gamma_dot[i];
  }
  rho_m_avg = rho_m_avg/_num_slip_sys;
  rho_i_avg = rho_i_avg/_num_slip_sys;
  gamma_dot_avg = gamma_dot_avg/_num_slip_sys;

  n = n + 1;
  _state_var[_qp][n] = rho_m_avg;

  n = n + 1;
  _state_var[_qp][n] = rho_i_avg;

  // Store average dislocation glide rates SDV 50
  n = n + 1;
  _state_var[_qp][n] = gamma_dot_avg;

  // Read mobile dislocation density values SDV 50 + _num_slip_sys
  for (unsigned int i = 0; i < _num_slip_sys; i++) {
    n = n + 1;
    _state_var[_qp][n] = rho_m[i];
  }

  // Read immobile dislocation density values SDV 50 + _num_slip_sys*2
  for (unsigned int i = 0; i < _num_slip_sys; i++) {
    n = n + 1;
    _state_var[_qp][n] = rho_i[i];
  }

  // Read backstress values SDV 50 + _num_slip_sys*3
  for (unsigned int i = 0; i < _num_slip_sys; i++) {
    n = n + 1;
    _state_var[_qp][n] = bstress[i];
  }

  // elasticity tensor
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      for (unsigned int k = 0; k < 3; ++k){
	    for (unsigned int l = 0; l < 3; ++l){
	      // Double check
	      _Cel_cp[_qp](i,j,k,l) = C[i][j][k][l];
	    }
      }
    }
  }

} // End computeStress()

// NR_residual()
void DDCPStressUpdate::NR_residual (unsigned int num_slip_sys, std::vector<std::vector<Real>> &xs0, std::vector<std::vector<Real>> &xm0, Real temp, Real dt, std::vector<Real> gdt, RankTwoTensor F1, RankTwoTensor &F_el, RankTwoTensor &F_p_inv, RankTwoTensor F_p_inv_0, Real C[3][3][3][3], std::vector<Real> &rho_m0, std::vector<Real> &rho_m, std::vector<Real> &rho_i0, std::vector<Real> &rho_i, std::vector<Real> &bstress0, std::vector<Real> &bstress, RankTwoTensor &sig, std::vector<Real> &tau, std::vector<Real> &tau_eff, std::vector<Real> &s_a, std::vector<Real> &s_t, std::vector<std::vector<Real>> A, std::vector<std::vector<Real>> H, std::vector<Real> &residual, Real &sse){

  Real xL_p_inter[3][3];
  std::vector<Real> shr_g(num_slip_sys);

  // RankTwoTensor F_el;
  RankTwoTensor E_el;
  RankTwoTensor Spk2;

  std::vector<std::vector<Real>> xs(3, std::vector<Real>(num_slip_sys));
  std::vector<std::vector<Real>> xm(3, std::vector<Real>(num_slip_sys));

  // Calculate the plastic part of the velocity gradient in the intermediate configuration
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      xL_p_inter[i][j] = 0.0;
      for (unsigned int k = 0; k < num_slip_sys; k++) {
        xL_p_inter[i][j] = xL_p_inter[i][j] + xs0[i][k]*xm0[j][k]*gdt[k];
      }
    }
  }

  // Begin calculation process of F_p_n+1 = exp(xL_p_inter*dt).F_p_n
  Real array1[3][3], array2[3][3];
  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      array1[i][j] = xL_p_inter[i][j]*dt;
    }
  }

  // Caculate Omega
  Real sum1 = 0;
  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      sum1 = sum1 + array1[i][j]*array1[i][j];
    }
  }
  Real omega = sqrt(0.5 * sum1);

  // Continue calculating intermediate stuff needed for F_p_n+1
  aa_dot_bb(array1,array1,array2);

  for(unsigned int i = 0; i < 3; i++) {
    if(omega > 0.0) { // if omega==0 then no need
    for(unsigned int j = 0; j < 3; j++) {
      array1[i][j] = array1[i][j]*sin(omega)/omega + array2[i][j]*(1.0-cos(omega))/power(omega, 2);
      }
    }
    array1[i][i] = 1.0 + array1[i][i];
  }

  // Finally multiply arrays to get F_p_inv at end of time step.
  RankTwoTensor array1_tensor, array2_tensor;

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      array1_tensor(i,j) = array1[i][j];
    }
  }

  array2_tensor = array1_tensor.inverse();
  F_p_inv = F_p_inv_0 * array2_tensor;

  // Multiply F() by F_p_inv() to get F_el()
  F_el = F1 * F_p_inv;
  if (F_el.det() == 0.0e0) {
    F_el(0,0) = 1.0e0;
    F_el(1,1) = 1.0e0;
    F_el(2,2) = 1.0e0;
  }
  RankTwoTensor F_el_inv = F_el.inverse();

  // Rotate director vectors from intermediate configuration to the current configuration.
  for(unsigned int n = 0; n < num_slip_sys; n++) {
    for(unsigned int i = 0; i < 3; i++) {
      xs[i][n] = 0.0;
      xm[i][n] = 0.0;
      for(unsigned int j = 0; j < 3; j++) {
        xs[i][n] = xs[i][n] + F_el(i,j)*xs0[j][n];
        xm[i][n] = xm[i][n] + xm0[j][n]*F_el_inv(j,i);
      }
    }
  }

  // Calculate elastic Green Strain
  E_el = F_el.transpose() * F_el;
  E_el(0,0) = E_el(0,0) - 1.0;
  E_el(1,1) = E_el(1,1) - 1.0;
  E_el(2,2) = E_el(2,2) - 1.0;
  E_el = 0.5 * E_el;

  // Multiply the stiffness tensor by the Green strain to get the 2nd Piola Kirkhhoff stress
  Real E_el_matrix[3][3], Spk2_matrix[3][3];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      E_el_matrix[i][j] = E_el(i,j);
    }
  }
  aaaa_dot_dot_bb(C,E_el_matrix,Spk2_matrix);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Spk2(i,j) = Spk2_matrix[i][j];
    }
  }

  // Convert from PK2 stress to Cauchy stress
  Real det = F_el.det();
  sig =  F_el * Spk2 * F_el.transpose();
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig(i,j) = sig(i,j)/det;
    }
  }

  // Calculate resolved shear stress for each slip system.
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    tau[k] = 0.0e0;
    for(unsigned int i = 0; i < 3; i++) {
      for(unsigned int j = 0; j < 3; j++) {
        tau[k] = tau[k] + xs[i][k]*xm[j][k]*sig(i,j);
      }
    }
  }

  // Calculate dislocation densities, rho_m, rho_i
  std::vector<Real> d_disl(num_slip_sys);
  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    Real sum1 = 0.0e0;
    for(unsigned int ib = 0; ib < num_slip_sys; ib++) {
      sum1 = sum1 + H[ia][ib]*(rho_m0[ib] + rho_i0[ib]);
    }
    d_disl[ia] = 1.0e0/sqrt(sum1);
  }

  std::vector<Real> drhomdt(num_slip_sys);
  std::vector<Real> drhoidt(num_slip_sys);
  std::vector<Real> dbstressdt(num_slip_sys);

  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    Real sum1 = 0.0e0;
    for(unsigned int ib = 0; ib < num_slip_sys; ib++) {
      sum1 = sum1 + H[ia][ib]*(rho_m0[ib] + rho_i0[ib]);
    }
    drhomdt[ia] = (k_M/b_mag/d_disl[ia] - k_ann*R_c*rho_m0[ia]/b_mag
      - (k_I*sqrt(sum1))/b_mag)*abs(gdt[ia]);

    drhoidt[ia] = ((k_I*sqrt(sum1))/b_mag - k_D*rho_i0[ia])*abs(gdt[ia]);
  }

  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    rho_m[ia] = rho_m0[ia] + drhomdt[ia]*dt;
    if(rho_m[ia] < 1.0e2) {
      rho_m[ia] = 1.0e2;
    }
    rho_i[ia] = rho_i0[ia] + drhoidt[ia]*dt;
    if(rho_i[ia] < 1.0e2) {
      rho_i[ia] = 1.0e2;
    }
  }

  // Calculate the back stress for each slip system.
  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    dbstressdt[ia] = (k_bs1*k_rho*G*b_mag*sqrt(rho_m[ia] + rho_i[ia])*sgn(tau[ia] - bstress0[ia]) - k_bs2*bstress0[ia])*abs(gdt[ia]);

    bstress[ia] = bstress0[ia] + dbstressdt[ia]*dt;
  }

  // Calculate the athermal slip resistance for each slip system.
  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    Real sum1 = 0.0;
    for(unsigned int ib = 0; ib < num_slip_sys; ib++) {
      sum1 = sum1 + p0*A[ia][ib]*(rho_m[ib] + rho_i[ib]);
    }
    s_a[ia] = tau0 + hp_coeff/sqrt(grain_size) + k_rho*G*b_mag*sqrt(sum1);
  }

  // Calculate thermal slip resistance for each slip system.
  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    s_t[ia] = frictional_stress;
  }

  // Calculate new slip rate for each slip system.
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    Real refgdg = gammadot0g;

    tau_eff[k] = tau[k] - bstress[k];

    if(abs(tau_eff[k]) > s_a[k]) {
      shr_g[k] = refgdg*exp((-delF0/B_k/temp)*power(1.0 - power((abs(tau_eff[k]) - s_a[k])/s_t[k], p), q))*sgn(tau_eff[k]);
    } else {
      shr_g[k] = 1e-18;
    }
  }

  // Calculate residual values
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    residual[k] = gdt[k] - shr_g[k];

    if (residual[k] >= 1.0e100) {
      residual[k] = 1.0e100;
    }
    if (residual[k] <= -1.0e100) {
      residual[k] = -1.0e100;
    }
  }

  // Calculate root mean square error that is to be minimized to zero
  sse = 0.0e0;
  Real gamma_dot_max = abs(gdt[0]);
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    sse = sse + abs(gdt[k])*(residual[k]*residual[k]);

    if(abs(gdt[k]) > gamma_dot_max) {
      gamma_dot_max = abs(gdt[k]);
    }
  }
  if (gamma_dot_max > 1.e-18) {
    sse = sqrt(sse/gamma_dot_max)/num_slip_sys;
  }

  // std::cout << "\n sse:" << sse;
} // end NR_residual

void DDCPStressUpdate::readPropsFile() {
  // MooseUtils::DelimitedFileReader reader(_propsFile);

  MooseUtils::checkFileReadable(_propsFile);
  std::ifstream file_prop;
  file_prop.open(_propsFile.c_str());

  // Assign slip system normals and slip directions for a BCC material
  for (unsigned int i = 0; i < _num_props; i++) {
    // file_prop >> _properties[_qp][i];
    if (!(file_prop >> _properties[_qp][i])) {
      mooseError("Required number of material parameters not supplied for crystal plasticity model");
    }
  }

  file_prop.close();
}

void DDCPStressUpdate::assignProperties(){

  // Elastic constants
  C11 = _properties[_qp][0];
  C11_perK = _properties[_qp][1];
  C12 = _properties[_qp][2];
  C12_perK = _properties[_qp][3];
  C44 = _properties[_qp][4];
  C44_perK = _properties[_qp][5];
  G = _properties[_qp][6];
  G_perK = _properties[_qp][7];

  C11 = C11 - C11_perK*_temp[_qp];
  C12 = C12 - C12_perK*_temp[_qp];
  C44 = C44 - C44_perK*_temp[_qp];
  G = G - G_perK*_temp[_qp];

  // use the appropriate one
  // G = 0.5e0*(C11 - C12);
  // G = sqrt(C44*0.5e0*(C11 - C12));
  // G = C44;

  b_mag = _properties[_qp][8];  // Burgers vector magnitude

  // Flow parameters
  gammadot0g = _properties[_qp][9];  // reference strain rate
  enthalpy_const = _properties[_qp][10];  // Multiplication constant for enthalpy term
  p = _properties[_qp][11];  // Shape parameter (dislocation glide)
  q = _properties[_qp][12];  // Shape parameter (dislocation glide)

  // Hardening parameters
  tau0 = _properties[_qp][13];  // Athermal slip resistance constant
  hp_coeff = _properties[_qp][14];  // Hall-Petch coefficient <TBA>
//   if (_grain_size >= 1.e10) {
  if (! _GrainAreaSize){
    grain_size = _properties[_qp][15];  // Grain size (for Hall-Petch term) <TBA>
  }
  else {
    // if the GrainAreaSize object exists, then grain_size is taken directly from the mesh
    grain_size = _grain_size;
  }
  frictional_stress = _properties[_qp][16]; // Lattice frictional resistance
  p0 = _properties[_qp][17];  // Dislocation barrier strength
  k_rho = _properties[_qp][18];  // Taylor hardening parameter
  k_I = _properties[_qp][19];  // Mean free path constant (dislocation network)

  Alatent = _properties[_qp][20];  // Latent hardening coefficient

  // Dislocation evolution parameters
  rho_m_zero = _properties[_qp][21];  // Initial mobile dislocation density
  rho_i_zero = _properties[_qp][22];  // Initial immobile dislocation density
  d_disl_zero = _properties[_qp][23];  // Initial dislocation line length
  k_M = _properties[_qp][24]; // Dislocation line generation constant
  R_c = _properties[_qp][25];  // Critical capture radius for dislocation annihilation
  k_ann = _properties[_qp][26];  // Dislocation evolution annihilation constant
  k_D = _properties[_qp][27];  // Immobile dislocation evolution dynamic recovery constant
  k_bs1 = _properties[_qp][28];  // Backstress evolution constant 1
  k_bs2 = _properties[_qp][29];  // Backstress evolution constant 2

  B_k = 1.3806503e-23;  // Boltzmann constant
  freq = 1e13;  // Debye frequency

  // adjusted to change units of product with stress to SI units
  act_vol = (pow(b_mag,3))*1.0e-3;

  // delF0 scaled to SI units
  delF0 = (enthalpy_const*G*1.0e6)*(pow(b_mag,3))*1.0e-9;
}

// normalize to unit vector
void DDCPStressUpdate::normalize_vector(Real *x, Real *y, Real *z) {
  Real length = sqrt(*x * *x + *y * *y + *z * *z);

  *x = *x/length;
  *y = *y/length;
  *z = *z/length;
}

// Rotate 4th order stiffness tensor
void DDCPStressUpdate::rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]) {
  Real d[3][3][3][3];

  for(unsigned int m = 0; m < 3; m++) {
    for(unsigned int n = 0; n < 3; n++) {
      for(unsigned int k = 0; k < 3; k++) {
        for(unsigned int l = 0; l < 3; l++) {
          d[m][n][k][l] = a[k][0]*(a[l][0]*b[m][n][0][0] + a[l][1]*b[m][n][0][1] + a[l][2]*b[m][n][0][2]) +
          a[k][1]*(a[l][0]*b[m][n][1][0] +a[l][1]*b[m][n][1][1] + a[l][2]*b[m][n][1][2]) +
          a[k][2]*(a[l][0]*b[m][n][2][0] +a[l][1]*b[m][n][2][1] + a[l][2]*b[m][n][2][2]);
        }
      }
    }
  }

  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      for(unsigned int k = 0; k < 3; k++) {
        for(unsigned int l = 0; l < 3; l++) {
          c[i][j][k][l] = a[i][0]*(a[j][0]*d[0][0][k][l] + a[j][1]*d[0][1][k][l] + a[j][2]*d[0][2][k][l]) +
          a[i][1]*(a[j][0]*d[1][0][k][l] + a[j][1]*d[1][1][k][l] + a[j][2]*d[1][2][k][l]) +
          a[i][2]*(a[j][0]*d[2][0][k][l] + a[j][1]*d[2][1][k][l] + a[j][2]*d[2][2][k][l]);
        }
      }
    }
  }
}

// 4th order tensor to Voigt notation
void DDCPStressUpdate::forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]) {
  for(unsigned int i = 1; i <= 3; i++) {
    for(unsigned int j = i; j <= 3; j++) {
      unsigned int ia = i;
      if(i != j) {
        ia = 9 - i - j;
      }
      for(unsigned int k = 1; k <= 3; k++) {
        for(unsigned int l = k; l <= 3; l++) {
          unsigned int ib = k;
          if(k != l) {
            ib = 9 - k - l;
          }
          b[ia-1][ib-1] = a[i-1][j-1][k-1][l-1];
        }
      }
    }
  }
}

// Voigt tensor to 4th order tensor
void DDCPStressUpdate::Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]) {
  for(unsigned int i = 1; i <= 3; i++) {
    for(unsigned int j = 1; j <= 3; j++) {
      unsigned int ia = i;
      if(i != j) {
        ia = 9 - i - j;
      }
      for(unsigned int k = 1; k <= 3; k++) {
        for(unsigned int l = 1; l <= 3; l++) {
          unsigned int ib = k;
          if(k != l) {
            ib = 9 - k - l;
          }
          a[i-1][j-1][k-1][l-1] = b[ia-1][ib-1];
          if(ia > 3) {
            a[i-1][j-1][k-1][l-1] = a[i-1][j-1][k-1][l-1]/2;
          }
          if(ib > 3) {
            a[i-1][j-1][k-1][l-1] = a[i-1][j-1][k-1][l-1]/2;
          }
        }
      }
    }
  }
}

void DDCPStressUpdate::aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          product[i][j][k][l] = 0.0;
          for (unsigned int m1 = 0; m1 < 3; m1++) {
            for (unsigned int m2 = 0; m2 < 3; m2++) {
              product[i][j][k][l] = product[i][j][k][l] + a[i][j][m1][m2]*b[m1][m2][k][l];
            }
          }
        }
      }
    }
  }
}

void DDCPStressUpdate::aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product[i][j] = 0.0;
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          product[i][j] = product[i][j] + a[i][j][k][l]*b[k][l];
        }
      }
    }
  }
}

void DDCPStressUpdate::aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product[i][j] = 0.0;
      for (unsigned int k = 0; k < 3; k++) {
        product[i][j] = product[i][j] + a[i][k]*b[k][j];
      }
    }
  }
}

Real DDCPStressUpdate::aa_dot_dot_bb(Real a[3][3], Real b[3][3]) {

  Real product = 0.0e0;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product = product + a[i][j]*b[i][j];
    }
  }
  return product;
}

// calculate Euler angles in Bunge notation
void DDCPStressUpdate::bunge_angles(Real (&array1)[3][3], Real (&psi0)[3]) {

  Real PI = 4.0*atan(1.0);
  Real phi_1, Phi, phi_2, sth;

  // From VPSC
  if (array1[2][2] >= 1.0e0) array1[2][2] = 1.0e0;
  if (array1[2][2] <= -1.0e0) array1[2][2] = -1.0e0;

  Phi = acos(array1[2][2]);

  if (abs(array1[2][2]) >= 0.99999) {
    phi_2 = 0.0;
    phi_1 = atan2(array1[0][1],array1[0][0]); // atan2(array1[0][1],array1[0][0]);
  }
  else {
    sth = sin(Phi);
    phi_2 = atan2(array1[0][2]/sth,array1[1][2]/sth);
    phi_1 = atan2(array1[2][0]/sth,-array1[2][1]/sth);
  }

  psi0[0] = 180.0e0*phi_1/PI;
  psi0[1] = 180.0e0*Phi/PI;
  psi0[2] = 180.0e0*phi_2/PI;
}

// return power (with sign)
Real DDCPStressUpdate::power(Real x, Real y) {

  Real ans;

  if (x == 0.0e0) {
    if (y > 0.0e0) {
      ans = 0.0e0;
    }
    else if (y < 0.0e0) {
      ans = 1.0e300;
    }
    else {
      ans = 1.0e0;
    }
  }
  else {
    ans = y*log10(abs(x));
    if (ans > 300) {
      ans = 1.0e300;
    }
    else {
      ans = pow(10.0e0,ans);
    }
    if (x < 0.0e0) {
      ans = -ans;
    }
  }

  return ans;
}

// return signum of a real number
Real DDCPStressUpdate::sgn(Real x) {

  Real ans;
  ans = 1.0;

  if (x < 0.0e0) {
    ans = -1.0;
  }

  return ans;
}

// return maximum of two scalars
Real DDCPStressUpdate::max_val(Real a, Real b)
{
  if (a>b)
  {
    return a;
  }
  else
  {
    return b;
  }
}
