#include "DDJ2StressUpdate.h"
#include "MatrixTools.h"
#include "MooseRandom.h"
#include "MooseException.h"
#include <algorithm>
#include <dlfcn.h>
#include <fstream>
#define QUOTE(macro) stringifyName(macro)

registerMooseObject("RhocpApp", DDJ2StressUpdate);

InputParameters
DDJ2StressUpdate::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Stress calculation using the J2 plasticity material model");
  params.addRequiredParam<FileName>("propsFile", "The file with the material parameters");
  params.addParam<UserObjectName>("GrainAreaSize", "Name of the GrainAreaSize UO");
  params.addRequiredParam<unsigned int>("num_props", "The number of material properties this UMAT is going to use");
  params.addRequiredParam<unsigned int>("num_state_vars", "The number of state variables this UMAT is going to use");
  params.addParam<Real>("tol", 1.0e-6, "Tolerance");
  params.addCoupledVar("temp", 300, "Temperature");
  return params;
}

DDJ2StressUpdate::DDJ2StressUpdate(const InputParameters & parameters) :
    ComputeStressBase(parameters),
    _propsFile(getParam<FileName>("propsFile")),
    _GrainAreaSize(isParamValid("GrainAreaSize")
                               ? &getUserObject<GrainAreaSize>("GrainAreaSize")
                               : NULL),
    _num_props(getParam<unsigned int>("num_props")),
    _num_state_vars(getParam<unsigned int>("num_state_vars")),
    _tol(getParam<Real>("tol")),
    _temp(coupledValue("temp")),
    _state_var(declareProperty<std::vector<Real> >("state_var")),
    _state_var_old(getMaterialPropertyOld<std::vector<Real> >("state_var")),
    _properties(declareProperty<std::vector<Real> >("properties")),
    _properties_old(getMaterialPropertyOld<std::vector<Real> >("properties")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _rotation_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "rotation_increment")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _Cel_cp(declareProperty<RankFourTensor>("Cel_cp"))
{
}

DDJ2StressUpdate::~DDJ2StressUpdate()
{
}

void DDJ2StressUpdate::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();

  //Initialize state variable vector
  _state_var[_qp].resize(_num_state_vars);

  // Initialize _properties vector
  _properties[_qp].resize(_num_props);

} // End initQpStatefulProperties()

void DDJ2StressUpdate::computeQpStress()
{
  Real PI = 4.0*atan(1.0);

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
  Real sig_vm, sig_vm_star;

  Real C0[3][3][3][3]; // Elastic stiffness tensor
  Real C[3][3][3][3];
  Real C_avg[3][3][3][3];
  Real dir_cos0[3][3];

  RankTwoTensor sig_avg;
  RankTwoTensor sig;
  RankTwoTensor Spk2;
  RankTwoTensor sdev, sdev_star;
  RankTwoTensor Np, Np_star;

  Real ddpdsig[3][3][3][3];
  Real ddsdde_4th[3][3][3][3];

  // ISVs
  Real rho_m0; // mobile dislocation density at beginning of step
  Real rho_m; // mobile dislocation density at end of step
  Real rho_i0; // immobile dislocation density at beginning of step
  Real rho_i; // immobile dislocation density at end of step
  RankTwoTensor bstress0; // backstress at beginning of step
  RankTwoTensor bstress; // backstress at end of step
  Real s_a; // athermal slip resistance
  Real s_t; // thermal slip resistance

  Real gamma_dot_trial; // shear strain rate
  Real gamma_dot_g; // shear strain rate due to glide

  Real residual; // residual during Newton-Raphson iteration
  Real del_epsilon, del_epsilon_old; // current and previous values of gamma_dot

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

  // Read grain size, if available
  _grainid = _current_elem->subdomain_id();

  if (_GrainAreaSize){
    _grain_size = _GrainAreaSize->getGrainSize(_grainid);
  }
  else{
    _grain_size = 1.e10;
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

  if (_t_step == 1) {
  // Initialize F_p_inv_0
  for (unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
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

  rho_m0 = rho_m_zero;
  rho_i0 = rho_i_zero;

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      bstress0(i,j) = 0.0e0;
    }
  }
  }
  //  End of initializations.  Read in internal variables

  if(_t_step > 1) {
    // checkpoint("Initializing ISVs, nonzero time step")
    int n = -1;
    // Read the elastic part of F SDV 1-9
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        // F_el[i][j] = _state_var[_qp][n]; // Not needed. This is only for post-processing
      }
    }

    // Read inverse of the plastic part of F SDV 10-18
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        F_p_inv_0(i,j) = _state_var[_qp][n];
      }
    }

    // Read E_p SDV 19-27
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        E_p(i,j) = _state_var[_qp][n];
      }
    }

    // Read E_eff SDV 28
    n = n + 1;
    E_eff = _state_var[_qp][n];

    // Read E_p_eff SDV 29
    n = n + 1;
    E_p_eff = _state_var[_qp][n];

    // Read Np_star SDV 30-38
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        Np_star(i,j) = _state_var[_qp][n];
      }
    }

    // Read del Epdot Effective value SDV 39
    n = n + 1;
    del_epsilon = _state_var[_qp][n];

    // Read average dislocation glide rates  SDV 40
    n = n + 1;
    // gamma_dot_g_avg0 = _state_var[_qp][n];

    // Read mobile dislocation density values SDV 41
    n = n + 1;
    rho_m0 = _state_var[_qp][n];

    // Read immobile dislocation density values SDV 42
    n = n + 1;
    rho_i0 = _state_var[_qp][n];

    // Read back stress values SDV 43-51
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        bstress0(i,j) = _state_var[_qp][n];
      }
    }

    // std::cout << "\nTotal no. of state vars read:" << n;
  } // End of initializations

  // Initialize ISOTROPIC 4th rank elastic stiffness tensor
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          C0[i][j][k][l] = (YM*poisson_ratio/(1.e0 + poisson_ratio)/(1.e0 - 2.e0*poisson_ratio))*del[i][j]*del[k][l] + (YM/2.e0/(1.e0 + poisson_ratio))*(del[i][k]*del[j][l]+del[i][l]*del[k][j]);
        }
      }
    }
  }

  // Initialize arrays for averaging over grains
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig_avg(i,j) = 0.0;
    }
  }

  // Begin calculations over the element
  // This is not used at the moment, but can used to rotate the elastic stiffness tensor
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      dir_cos0[i][j] = 0.e0;
    }
    dir_cos0[i][i] = 1.e0;
  }

  rotate_4th(dir_cos0, C0, C);

  // Initialize number of sub-increments.  Note that the subincrement initially equals the total increment. This remains the case unless the process starts to diverge.
  unsigned int N_incr, N_incr_total, iNR, N_ctr;
  Real dt_incr, depdse;
  RankTwoTensor dfgrd0, dfgrd1;
  Real array1[3][3], array2[3][3];
  Real array1_matrix[3][3], array2_matrix[3][3];

  bool converged, improved;

  Real isubincr = 0;

  // Initialize deformation gradients for beginning and end of subincrement.
  dfgrd0 = _deformation_gradient_old[_qp];
  dfgrd1 = _deformation_gradient[_qp];
  F0 = dfgrd0;
  F1 = dfgrd1;

  // counters for time step sub-increment
  N_incr = 1;
  N_incr_total = 1;

  // Begin time step subincrement
  do {
  dt_incr = _dt/N_incr_total;

  F0 = dfgrd0 + (dfgrd1 - dfgrd0)*(N_incr - 1)/N_incr_total;
  F1 = dfgrd0 + (dfgrd1 - dfgrd0)*N_incr/N_incr_total;

  // Initialize mobile dislocation density
  rho_m = rho_m0;

  // Initialize immobile dislocation density
  rho_i = rho_i0;

  // Initialize backstress
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      bstress(i,j) = bstress0(i,j);
  }

  // Multiply F() by F_p_inv() to get F_el()
  F_el = F0 * F_p_inv_0;
  if (F_el.det() == 0.0e0) {
    F_el(0,0) = 1.0e0;
    F_el(1,1) = 1.0e0;
    F_el(2,2) = 1.0e0;
  }

  F_el_inv = F_el.inverse();

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

  // Calculate deviatoric stress tensor
  Real trace_sig = (sig(0,0)+sig(1,1)+sig(2,2))/3.0;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sdev(i,j) = sig(i,j);
    }
    sdev(i,i) = sdev(i,i) - trace_sig;
  }

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sdev_star(i,j) = sig(i,j) - bstress(i,j);
    }
    sdev_star(i,i) = sdev_star(i,i) - trace_sig;
  }

  // Calculate von Mises effective stress
  Real s1 = 0.0; Real s2 = 0.0;

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      s1 = s1 + sdev(i,j)*sdev(i,j);
      s2 = s2 + sdev_star(i,j)*sdev_star(i,j);
    }
  }

  sig_vm = sqrt(s1*3.0/2.0);
  sig_vm_star = sqrt(s2*3.0/2.0);

  if (sig_vm <= 1.e-20) sig_vm = 1.e-20;
  if (sig_vm_star <= 1.e-20) sig_vm_star = 1.e-20;

  // Calculate normal to the Cauchy stress tensor
  Real v;
  v = sqrt(3.0/2.0)*(1.0/sig_vm);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Np(i,j) = v*sdev(i,j);
    }
  }

  v = sqrt(3.0/2.0)*(1.0/sig_vm_star);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Np_star(i,j) = v*sdev_star(i,j);
    }
  }

  // Calculate athermal slip resistance
  Real sum1 = 0.0;
  sum1 = p0*(rho_m + rho_i);
  s_a = tau0 + hp_coeff/sqrt(grain_size) + k_rho*G*b_mag*sqrt(sum1);

  // Calculate reference shear stress
  s_t = frictional_stress;

  // Calculate 1st estimate of gamma_dot
  Real refgdg = gammadot0g;
  if(sig_vm_star > s_a) {
      gamma_dot_g = refgdg*exp((-delF0/B_k/_temp[_qp])*power(1 - power((sig_vm_star - s_a)/s_t, p), q));
    }
    else {
      gamma_dot_g = 1e-18;
    }
    gamma_dot_trial = gamma_dot_g;
  }

  // Begin Newton-Raphson iterative loops.
  iNR = 0;

  converged = false;

  do {
    iNR = iNR + 1;

    if (iNR == 1) {
        del_epsilon_old = 1.e0;
    }
    else {
      del_epsilon_old = del_epsilon;
    }

    Real drhomded, drhoided, dsadrho, dsaded, dstded, dfded, dvmded, xlambda, d_disl;

    // TODO
    NR_residual_J2 (_temp[_qp], dt_incr, F1, F_p_inv, F_p_inv_0, C, Np_star, rho_m0, rho_m, rho_i0, rho_i, bstress0, bstress, sig, sig_vm, sig_vm_star, s_a, s_t, gamma_dot_trial, residual);

//       if (sse > 0.0e0) {
//         std::cout << "\n iNR:" << iNR;
//         std::cout << "\n sse:" << sse << "\n";
//       }

    // Begin calculation of the partial derivatives needed for the Newton-Raphson step

    // Calculate d(sig_vm)/d(Epsilon_dot)
    Real Np_star_matrix[3][3];
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        Np_star_matrix[i][j] = Np_star(i,j);
      }
    }

    aa_dot_dot_bbbb (Np_star_matrix,C,array1);
    dvmded = aa_dot_dot_bb (array1,Np_star_matrix);
    dvmded = -1.5*dvmded;

    // Calculate d(s_t)/d(Epsilon_dot)
    dstded = 0.e0;

    // Calculate d(s_a)/d(Epsilon_dot)
    dsaded = 0.e0;
    drhomded = 0.e0;
    drhoided = 0.e0;

    Real sum1 = 0.e0;
    sum1 = rho_m + rho_i;
    d_disl = 1.0e0/sqrt(rho_m + rho_i);

    drhomded = (k_M/b_mag/d_disl - k_ann*R_c*rho_m/b_mag - (k_I*sqrt(sum1))/b_mag)*dt_incr;

    drhoided = (k_I*sqrt(sum1)/b_mag - k_D*rho_i)*dt_incr;

    dsadrho = 0.5*k_rho*G*b_mag*sqrt(p0/sum1);

    dsaded = dsadrho*(drhoided + drhomded);

    // Calculate the value of d(function)/d(epsilon_dot)
    xlambda =(dvmded - dsaded)/s_t - dstded*(sig_vm_star - s_a)/(s_t*s_t);

    dfded = 1.0 - gamma_dot_trial*(delF0/B_k/_temp[_qp])*q*power(1.0 - power((abs(sig_vm_star) - s_a)/s_t,p), q - 1)*p*power((abs(sig_vm_star) - s_a)/s_t, p - 1)*xlambda;

    // Calculate the del_epsilon
    del_epsilon = -residual/dfded;

    //  Update gamma_dot_trial
    gamma_dot_trial = gamma_dot_trial + del_epsilon;

    // Check whether the del_epsilon is less than tolerance
    E_tot = F1.transpose() * F1;
    E_tot(0,0) = E_tot(0,0) - 1.0;
    E_tot(1,1) = E_tot(1,1) - 1.0;
    E_tot(2,2) = E_tot(2,2) - 1.0;
    E_tot = 0.5 * E_tot;

    Real s1 = 0.0;

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        s1 = s1 + E_tot(i,j)*E_tot(i,j);
      }
    }
    E_eff = sqrt((2./3.) * s1);

    if (E_eff < 1.e-20) E_eff = 1.e-20;

    if ((del_epsilon*dt_incr/E_eff) <= _tol) {
      converged = true;
    }
    else if ((del_epsilon > 0.7*del_epsilon_old) && isubincr == 1) {
      N_incr = 2*N_incr - 1;
      N_incr_total = 2*N_incr_total;

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

    throw MooseException("J2 plasticity Newton Raphson did not converge");
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

    rho_m0 = rho_m;
    rho_i0 = rho_i;
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        bstress0(i,j) = bstress(i,j);
      }
    }
  }
  else if (N_incr == N_incr_total) {
    break;
  }
  } while (N_incr_total < 10000); // End time step subincrementation

  if (N_incr_total > 10000) {
    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        _stress[_qp](i,j) = MooseRandom::rand()*5.0e3;
        for (unsigned int k = 0; k < 3; k++){
          for (unsigned int l = 0; l < 3; l++){
            _Jacobian_mult[_qp](i,j,k,l) = MooseRandom::rand()*5.0e3;
          }
        }
      }
    }

    throw MooseException("J2 plasticity time step subincrement did not converge");
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

  Real v;
  if (sig_vm_star <= 1.e-20) sig_vm_star = 1.e-20;

  v = sqrt(1.5)*(1.0/sig_vm_star);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Np_star(i,j) = v*sdev_star(i,j);
    }
  }

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          depdse = gamma_dot_trial*(delF0/B_k/_temp[_qp])*q*(p/s_t)*   power(1.e0 - power((abs(sig_vm_star)-s_a)/s_t,p),    q-1)*power((abs(sig_vm_star)-s_a)/s_t,p-1);

          ddpdsig[i][j][k][l] = ddpdsig[i][j][k][l] + 1.5e0*Np_star(i,j)*    Np_star(k,l)*depdse;
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

  // checkpoint("Storing ISVs")
  int n = -1;
  // Store the elastic part of F SDV 1-9
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = F_el(i,j);
    }
  }

  // Store inverse of the plastic part of F SDV 10-18
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

  // Store E_p SDV 19-27
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = E_p(i,j);
    }
  }

  // Effective strain calculations
  Real s1 = 0.0; Real s2 = 0.0;

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      s1 = s1 + E_tot(i,j)*E_tot(i,j);
      s2 = s2 + gamma_dot_trial*gamma_dot_trial*Np_star(i,j)*Np_star(i,j)*_dt*_dt;
    }
  }

  // Effective total strain
  E_eff = sqrt((2./3.) * s1);

  // Effective plastic strain
  E_p_eff = E_p_eff + sqrt((2./3.) * s2);

  // Store E_eff SDV 28
  n = n + 1;
  _state_var[_qp][n] = E_eff;

  // Store E_p_eff SDV 29
  n = n + 1;
  _state_var[_qp][n] = E_p_eff;

  // Store Np_star  SDV 30-38
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = Np_star(i,j);
    }
  }

  // Store del Epdot Effective value SDV 39
  n = n + 1;
  _state_var[_qp][n] = del_epsilon;

  // Store average dislocation glide rates SDV 40
  n = n + 1;
  _state_var[_qp][n] = gamma_dot_trial;

  // Store dislocation densities SDV 41-42
  n = n + 1;
  _state_var[_qp][n] = rho_m;

  n = n + 1;
  _state_var[_qp][n] = rho_i;

  // Read backstress values SDV 43-51
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = bstress(i,j);
    }
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

// NR_residual_J2()
void DDJ2StressUpdate::NR_residual_J2 (Real temp, Real dt, RankTwoTensor F1, RankTwoTensor &F_p_inv, RankTwoTensor F_p_inv_0, Real C[3][3][3][3], RankTwoTensor Np_star, Real &rho_m0, Real &rho_m, Real &rho_i0, Real &rho_i, RankTwoTensor &bstress0, RankTwoTensor &bstress, RankTwoTensor &sig, Real &sig_vm, Real &sig_vm_star, Real &s_a, Real &s_t, Real &gamma_dot_trial, Real &residual){

  Real xL_p_inter[3][3];

  RankTwoTensor F_el;
  RankTwoTensor E_el;
  RankTwoTensor Spk2, sdev, sdev_star;
  RankTwoTensor Np;

  Real gamma_dot_NR;

  // Calculate the plastic part of the velocity gradient in the intermediate configuration
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      xL_p_inter[i][j] = sqrt(1.5e0)*gamma_dot_trial*Np_star(i,j);
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

  // Calculate deviatoric stress tensor
  Real trace_sig = (sig(0,0)+sig(1,1)+sig(2,2))/3.0;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sdev(i,j) = sig(i,j);
    }
    sdev(i,i) = sdev(i,i) - trace_sig;
  }

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sdev_star(i,j) = sig(i,j) - bstress(i,j);
    }
    sdev_star(i,i) = sdev_star(i,i) - trace_sig;
  }

  // Calculate von Mises effective stress
  Real s1 = 0.0; Real s2 = 0.0;

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      s1 = s1 + sdev(i,j)*sdev(i,j);
      s2 = s2 + sdev_star(i,j)*sdev_star(i,j);
    }
  }

  sig_vm = sqrt(s1*3.0/2.0);
  sig_vm_star = sqrt(s2*3.0/2.0);

  if (sig_vm <= 1.e-20) sig_vm = 1.e-20;
  if (sig_vm_star <= 1.e-20) sig_vm_star = 1.e-20;

  // Calculate normal to the Cauchy stress tensor
  Real v;
  v = sqrt(3.0/2.0)*(1.0/sig_vm);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Np(i,j) = v*sdev(i,j);
    }
  }

  v = sqrt(3.0/2.0)*(1.0/sig_vm_star);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Np_star(i,j) = v*sdev_star(i,j);
    }
  }

  // Calculate dislocation densities, rho_m, rho_i
  Real d_disl;
  sum1 = rho_m0 + rho_i0;
  d_disl = 1.0e0/sqrt(sum1);

  Real drhomdt, drhoidt;
  RankTwoTensor dbstressdt;

  drhomdt = (k_M/b_mag/d_disl - k_ann*R_c*rho_m0/b_mag - (k_I*sqrt(sum1))/b_mag)*abs(gamma_dot_trial);

  drhoidt = ((k_I*sqrt(sum1))/b_mag - k_D*rho_i0)*abs(gamma_dot_trial);

  rho_m = rho_m0 + drhomdt*dt;
  if(rho_m < 1.0e2) {
    rho_m = 1.0e2;
  }
  rho_i = rho_i0 + drhoidt*dt;
  if(rho_i < 1.0e2) {
    rho_i = 1.0e2;
  }

  // Calculate the back stress
  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      dbstressdt(i,j) = (k_bs1*k_rho*G*b_mag*sqrt(rho_m0 + rho_i0)*Np_star(i,j) - k_bs2*bstress0(i,j))*abs(gamma_dot_trial);

      bstress(i,j) = bstress0(i,j) + dbstressdt(i,j)*dt;
    }
  }

  // Calculate the athermal slip resistance
  sum1 = p0*(rho_m + rho_i);
  s_a = tau0 + hp_coeff/sqrt(grain_size) + k_rho*G*b_mag*sqrt(sum1);

  // Calculate thermal slip resistance
  s_t = frictional_stress;

  // Calculate new slip rate
  Real refgdg = gammadot0g;
  if(sig_vm_star > s_a) {
    gamma_dot_NR = refgdg*exp((-delF0/B_k/_temp[_qp])*power(1 - power((sig_vm_star - s_a)/s_t, p), q));
  }
  else {
    gamma_dot_NR = 1e-18;
  }

  // Calculate residual
  residual = gamma_dot_trial - gamma_dot_NR;

} // end NR_residual_J2

void DDJ2StressUpdate::readPropsFile() {
  // MooseUtils::DelimitedFileReader reader(_propsFile);

  MooseUtils::checkFileReadable(_propsFile);
  std::ifstream file_prop;
  file_prop.open(_propsFile.c_str());

  for (unsigned int i = 0; i < _num_props; i++) {
    // file_prop >> _properties[_qp][i];
    if (!(file_prop >> _properties[_qp][i])) {
      mooseError("Required number of material parameters not supplied for J2 plasticity model");
    }
  }

  file_prop.close();
}

void DDJ2StressUpdate::assignProperties(){

  // Elastic constants
  YM = _properties[_qp][0];
  YM_perK = _properties[_qp][1];
  poisson_ratio = _properties[_qp][2];
  poisson_ratio_perK = _properties[_qp][3];
  G = _properties[_qp][4];
  G_perK = _properties[_qp][5];

  YM = YM - YM_perK*_temp[_qp];
  poisson_ratio = poisson_ratio - poisson_ratio_perK*_temp[_qp];
  G = G - G_perK*_temp[_qp];

  b_mag = _properties[_qp][6];  // Burgers vector magnitude

  // Flow parameters
  gammadot0g = _properties[_qp][7];  // reference strain rate
  enthalpy_const = _properties[_qp][8];  // Multiplication constant for enthalpy term
  p = _properties[_qp][9];  // Shape parameter (dislocation glide)
  q = _properties[_qp][10];  // Shape parameter (dislocation glide)

  // Hardening parameters
  tau0 = _properties[_qp][11];  // Athermal slip resistance constant
  hp_coeff = _properties[_qp][12];  // Hall-Petch coefficient <TBA>
//   if (_grain_size >= 1.e10) {
  if (! _GrainAreaSize){
    grain_size = _properties[_qp][13];  // Grain size (for Hall-Petch term) <TBA>
  }
  else {
    // if the GrainAreaSize object exists, then grain_size is taken directly from the mesh
    grain_size = _grain_size;
  }
  frictional_stress = _properties[_qp][14]; // Lattice frictional resistance
  p0 = _properties[_qp][15];  // Dislocation barrier strength
  k_rho = _properties[_qp][16];  // Taylor hardening parameter
  k_I = _properties[_qp][17];  // Mean free path constant (dislocation network)

  // Dislocation evolution parameters
  rho_m_zero = _properties[_qp][18];  // Initial mobile dislocation density
  rho_i_zero = _properties[_qp][19];  // Initial immobile dislocation density
  d_disl_zero = _properties[_qp][20];  // Initial dislocation line length
  k_M = _properties[_qp][21]; // Dislocation line generation constant
  R_c = _properties[_qp][22];  // Critical capture radius for dislocation annihilation
  k_ann = _properties[_qp][23];  // Dislocation evolution annihilation constant
  k_D = _properties[_qp][24];  // Immobile dislocation evolution dynamic recovery constant
  k_bs1 = _properties[_qp][25];  // Backstress evolution constant 1
  k_bs2 = _properties[_qp][26];  // Backstress evolution constant 2

  B_k = 1.3806503e-23;  // Boltzmann constant
  freq = 1e13;  // Debye frequency

  // adjusted to change units of product with stress to SI units
  act_vol = (pow(b_mag,3))*1.0e-3;

  // delF0 scaled to SI units
  delF0 = (enthalpy_const*G*1.0e6)*(pow(b_mag,3))*1.0e-9;
}

// normalize to unit vector
void DDJ2StressUpdate::normalize_vector(Real *x, Real *y, Real *z) {
  Real length = sqrt(*x * *x + *y * *y + *z * *z);

  *x = *x/length;
  *y = *y/length;
  *z = *z/length;
}

// Rotate 4th order stiffness tensor
void DDJ2StressUpdate::rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]) {
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
void DDJ2StressUpdate::forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]) {
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
void DDJ2StressUpdate::Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]) {
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

void DDJ2StressUpdate::aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]) {

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

void DDJ2StressUpdate::aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]) {

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

void DDJ2StressUpdate::aa_dot_dot_bbbb(Real a[3][3], Real b[3][3][3][3], Real (&product)[3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product[i][j] = 0.0;
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          product[i][j] = product[i][j] + a[k][l]*b[k][l][i][j];
        }
      }
    }
  }
}


void DDJ2StressUpdate::aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product[i][j] = 0.0;
      for (unsigned int k = 0; k < 3; k++) {
        product[i][j] = product[i][j] + a[i][k]*b[k][j];
      }
    }
  }
}

Real DDJ2StressUpdate::aa_dot_dot_bb(Real a[3][3], Real b[3][3]) {

  Real product = 0.0e0;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product = product + a[i][j]*b[i][j];
    }
  }
  return product;
}

// return power (with sign)
Real DDJ2StressUpdate::power(Real x, Real y) {

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
Real DDJ2StressUpdate::sgn(Real x) {

  Real ans;
  ans = 1.0;

  if (x < 0.0e0) {
    ans = -1.0;
  }

  return ans;
}

// return maximum of two scalars
Real DDJ2StressUpdate::max_val(Real a, Real b)
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
