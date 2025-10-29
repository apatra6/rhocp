// Developed by Vikram Roy

#include "ThermalIrradiationCPUpdate.h"
#include "MatrixTools.h"
#include "MooseRandom.h"
#include "MooseException.h"
#include <algorithm>
#include <dlfcn.h>
#include <fstream>
#define QUOTE(macro) stringifyName(macro) 

registerMooseObject("RhocpApp", ThermalIrradiationCPUpdate);

InputParameters
ThermalIrradiationCPUpdate::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Stress calculation using the crystal plasticity material model for thermal and irradiation hardening and creep");
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
  params.addParam<bool>("climbmodel", false, "Consider Climb Strains in Analysis");
  params.addParam<bool>("Transition_slip", false, "Transition Slip System Activity");
  params.addParam<FileName>("intMatFile","", "The file with the interaction matrix coefficients");
  params.addParam<bool>("modify_intmat", false, "Modify the Intmat Values Based on Change in Dislocation Density");
  params.addParam<std::vector<Real>>("coefficients",{0}, "Coefficients a,b & length conversion factor to m, for modifying intmat with change in dislocation density");
  params.addParam<unsigned int>("sliplaw", 1, "Slip Law Formulation for Calculations");
  params.addParam<unsigned int>("crsslaw", 1, "CRSS Law Formulation for Calculations");
  params.addParam<bool>("irad_defects", false, "Are Irradiation Defects present in calculations ?");
  params.addParam<bool>("deltaH_eV", false, "deltaH in Props file given in eV instead of multiplier.");
  return params;
}

ThermalIrradiationCPUpdate::ThermalIrradiationCPUpdate(const InputParameters & parameters) :
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
    _climbmodel(getParam<bool>("climbmodel")),
    _transition_slip_on(getParam<bool>("Transition_slip")),
    _intMatfile(getParam<FileName>("intMatFile")),
    _modify_intmat(getParam<bool>("modify_intmat")),
    _coefficients(getParam<std::vector<Real>>("coefficients")),    
    _sliplaw(getParam<unsigned int>("sliplaw")),
    _crsslaw(getParam<unsigned int>("crsslaw") ),
    _irad_defects(getParam<bool>("irad_defects")),
    _deltaH_eV(getParam<bool>("deltaH_eV")),
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
  
  if(!_intMatfile.empty()){
    readfile(_intMatfile);
  }
}

ThermalIrradiationCPUpdate::~ThermalIrradiationCPUpdate()
{
}

void ThermalIrradiationCPUpdate::initQpStatefulProperties()
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

void ThermalIrradiationCPUpdate::computeQpStress()
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
  std::vector<Real> gamma_dot_c (_num_slip_sys); // Shear Strain Rate due to climb
  std::vector<Real> gamma_try(_num_slip_sys); // trial shear strain rate
  std::vector<Real> gamma_try_g(_num_slip_sys); // trial shear strain rate due to glide
  std::vector<Real> gamma_try_c(_num_slip_sys); // trial shear strain rate due to climb

  std::vector<Real> residual(_num_slip_sys); // residual during Newton-Raphson iteration
  std::vector<std::vector<Real>> dTaudGd(_num_slip_sys, std::vector<Real>(_num_slip_sys)); // d(tau)/d(gammadot_g)

  // (Vik) New Vars for Climb
  Real interdensity, vacdensity; // Density of Interstitials and Vacancies
  Real Di, Dv; // Interstial and Vacancy Diffusivity at the Temp of Calculation
  Real vacdensity_th; // Zero stress thermal equilibrium density of vacancies
  Real lg; // Mean Free Path of Dislocation Climb 
  std::vector<Real> climbstress(_num_slip_sys);
  std::vector<Real> zi(_num_slip_sys); // Finite Stress Interstitial Capture Efficency
  std::vector<Real> zv(_num_slip_sys); // Finite Stress Vacancy Capture Efficency
  std::vector<Real> Kis(_num_slip_sys); // Insterstial sink (Dislocation) reaction rate coefficient
  std::vector<Real> Kvs(_num_slip_sys); // Vacancy sink (Dislocation) Reaction rate coefficient
  std::vector<Real> vel_climb(_num_slip_sys); // Climb Velocity
  std::vector<Real> vacdensity_0(_num_slip_sys); // Density of Vacancies on each slip system 
  std::vector<std::vector<Real>> dEtadGc(_num_slip_sys, std::vector<Real>(_num_slip_sys)); // d(climbstress)/d(gammadot_climb)   

  // Variables for Irradiation Defects
  std::vector<Real> Nloop0(_num_slip_sys);  // Number Density of Irradiation Loops - Initial Value Read from Input File/SDV
  std::vector<Real> dloop0(_num_slip_sys);  // Diameter of Irradiation Loops - Initial Value Read from Input File/SDV
  std::vector<Real> Nloop(_num_slip_sys);  // Number Density of Irradiation Loops - Final Value to be stored into SDV
  std::vector<Real> dloop(_num_slip_sys);  // Diameter of Irradiation Loops - Final Value to be stored into SDV
  
  


  // interaction coefficient martix between different slip systems
  std::vector<std::vector<Real>> A(_num_slip_sys, std::vector<Real>(_num_slip_sys));
  std::vector<std::vector<Real>> H(_num_slip_sys, std::vector<Real>(_num_slip_sys));

  // Resize Some Vectors if Climb is considered
  if(_climbmodel){
    residual.resize(2*_num_slip_sys);
  }
  else {
    interdensity = 0;
    vacdensity = 0;
    Di = 0;
    Dv = 0;
    vacdensity_th = 0;
    for(unsigned int k = 0; k < _num_slip_sys; k++){
      gamma_dot_c[k] = 0;
      gamma_try_c[k] = 0;
      climbstress[k] = 0;
      zi[k] = 0;
      zv[k] = 0;
      Kis[k] = 0;
      Kvs[k] = 0;
      vel_climb[k] = 0;
      vacdensity_0[k] = 0;
    }
  }


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

  //Read SintMat Data, if available
  if(!_intMatfile.empty()){   

   if(_modify_intmat){
    if(_coefficients.size() != 3){
      mooseError("variable coefficient should contain only three Real values when modify_intmat = true");
    }
   }
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
  if(!_intMatfile.empty()){
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      for (unsigned int j = 0; j < _num_slip_sys; j++) {
        A[i][j] = sintmat[i][j];
        H[i][j] = 1.0;
      }
      H[i][i] = 1.0;
    }
  }
  else {
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      for (unsigned int j = 0; j < _num_slip_sys; j++) {
        A[i][j] = p0*Alatent;
        H[i][j] = 1.0;
      }
      A[i][i] = p0*1.0;
      H[i][i] = 1.0;
    }    
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

    if(_irad_defects){
    for (unsigned int n = 0; n <_num_slip_sys; n++){
      Nloop0[n] = Nloop_zero/(_num_slip_sys);
      dloop0[n] = dloop_zero;
    }
    }
    else {
    for (unsigned int n = 0; n <_num_slip_sys; n++){
      Nloop0[n] = 1e-18;
      dloop0[n] = 1e-18;
      beta_loop = 1;
    }      
    }

    if(_transition_slip_on && _temp[_qp] < transition_temp){
      enthalpy_const = enthalpy_const2;
      p = p2;
      q = q2;
   if(_deltaH_eV){ delF0 = enthalpy_const*1.6e-19;  }
    else {    // delF0 scaled to SI units
    delF0 = (enthalpy_const*G*1.0e6)*(pow(b_mag,3))*1.0e-9;
    }      
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

    // Read Dislocation Loop Densities from SDV 50 + _num_slip_sys*4
    if(_irad_defects){
      for (unsigned int i = 0; i < _num_slip_sys; i++) {
        n = n + 1;
        Nloop0[i] = _state_var[_qp][n];
      }

    //  Read Dislocation Loop Densities from SDV 50 + _num_slip_sys*5
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
        n = n + 1;
        dloop0[i] = _state_var[_qp][n];
      }   
    }
    else {
      for (unsigned int i = 0; i < _num_slip_sys; i++){
        Nloop0[i] = 1e-18;
        dloop0[i] = 1e-18;
      }
    }
    
    if(_transition_slip_on && _temp[_qp] < transition_temp){
      enthalpy_const = enthalpy_const2;
      p = p2;
      q = q2;
   if(_deltaH_eV){ delF0 = enthalpy_const*1.6e-19;  }
    else {    // delF0 scaled to SI units
    delF0 = (enthalpy_const*G*1.0e6)*(pow(b_mag,3))*1.0e-9;
    }  

   C11 = C11P - C11_perK*_temp[_qp];
   C12 = C12P - C12_perK*_temp[_qp];
   C44 = C44P - C44_perK*_temp[_qp];
   G = GP - G_perK*_temp[_qp];          
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

  // // Modify the Interaction Coefficients if specified in the inputfile
  // if(_modify_intmat && !_intMatfile.empty()){
  //   std::vector<std::vector<Real>> sintmatm(_num_slip_sys, std::vector<Real>(_num_slip_sys));
  //   Real rhoref{1e12};
  //   for(unsigned int n = 0; n < _num_slip_sys; n++) {
  //     for(unsigned int i = 0; i < _num_slip_sys ; i++){
  //       Real rhoforest{0};
  //       for(unsigned int j = 0; j < _num_slip_sys; j++){
  //         if(i == j) 
  //           {continue;}
  //         else {
  //           rhoforest += rho_m0[j];
  //         }
  //       }
  //       rhoforest = rhoforest/(_coefficients[2] * _coefficients[2]); // Convert Rhoforest to SI units

  //       Real temp1 = log(0.35 * b_mag * _coefficients[2] * sqrt(rhoref));
  //       Real temp2 = log(0.35 * b_mag * _coefficients[2] * sqrt(rhoforest));
  //       Real temp3 = temp2/temp1;
  //       Real temp4 = (_coefficients[0] + _coefficients[1]*temp3);
  //       Real cfac = temp4 * temp4;

  //       sintmatm[n][i] = sintmat[n][i]*cfac;
  //       A[n][i] = sintmatm[n][i];

  // }
  // }
  // }

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

  // Initialize number of sub-increments.  Note that the subincrement initially equals the total increment. 
  // This remains the case unless the process starts to diverge.
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
      Ftheta(i,i) = sqrt(1.0 + 2.0* total_eigenstrain(i,i));
      Ftheta_old(i,i) = sqrt(1.0 + 2.0* total_eigenstrain_old(i,i));
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

  // Initialize Number Density and Size of Dislocation Loops
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    Nloop[n] = Nloop0[n];
    dloop[n] = dloop0[n];
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
  calc_crss(_num_slip_sys, rho_m, rho_i, Nloop, dloop, s_a, A);


  // Calculate reference shear stress for each slip system.
  for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    s_t[ia] = frictional_stress;
  }

  // Calculate 1st estimate of gamma_dot for each slip system.
  std::vector<Real> d_disl(_num_slip_sys); // dislocation spacing
  std::vector<Real> d_disl4(_num_slip_sys); // Effective Trap Mean Free Path
  calc_mfp(_num_slip_sys, rho_m, rho_i, Nloop, dloop, H, d_disl);
    d_disl4 = d_disl;



  calc_sliprate(_num_slip_sys, tau_eff, s_a, s_t, gamma_dot_g, gamma_dot);


  // Calculate d(Tau)/d(gammadot_g)
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


  // (Vik) Calculate the resolved climb stresses on each of the slip system
  if(_climbmodel){
    for (unsigned int k = 0; k < _num_slip_sys; k++) {
        climbstress[k] = 0.0;
        for (unsigned int j = 0; j< 3; j++) {
            for (unsigned int i = 0; i < 3 ; i++) {
              climbstress[k] = climbstress[k] + xs[i][k]*xs[j][k]*sig(i,j);
            }
        }
    }

    // (Vik) Variables for Calculation
    Real array3_matrix[3][3], array4_matrix[3][3];
    Real dummy1, dummy2;

    // (Vik) Calculate the Di & Dv
    Di = Di0 * exp(-((Em_i*Ecorr)/(B_k * _temp[_qp])));
    Dv = Dv0 * exp(-((Em_v*Ecorr)/(B_k * _temp[_qp])));

    // (Vik) Calculate the zero stress thermal equilibrium Density of vacancies 
    vacdensity_th = Nv * exp((-G0*Ecorr) / (B_k * _temp[_qp]));
    
    // (Vik) Calculate The ClimbRate
    calc_climbrate(_num_slip_sys, Di, Dv, vacdensity_th, interdensity, vacdensity, climbstress,  zi, zv, Kis, Kvs, vacdensity_0, vel_climb, rho_m, rho_i, gamma_dot_c);


  // (Vik)Calculating d(climbstress)/d(gammadot_climb)
  for (unsigned int ia = 0; ia <_num_slip_sys; ia++){
    for (unsigned int ib = 0; ib< _num_slip_sys; ib++) {
      dEtadGc[ia][ib] = 0.0;

      for(unsigned int i = 0; i< 3; i++) {
        for(unsigned int j = 0; j< 3; j++){ 
        array3_matrix[i][j] = xs0[i][ib]*xs0[j][ib];
        }
      }
      aaaa_dot_dot_bb(C0, array3_matrix, array4_matrix);

      for (unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 3; j++) {
          array3_matrix[i][j] = xs0[i][ia]*xs0[i][ia];
        }
      }
      dEtadGc[ia][ib] = aa_dot_dot_bb(array3_matrix, array4_matrix);
      dEtadGc[ia][ib] = -1.0*dEtadGc[ia][ib]* dt_incr;
    }
  }
  }
  else {
      Di = 0;
      Dv = 0;
      vacdensity_th = 0;
      interdensity = 0;
      vacdensity = 0;
    for (unsigned int i = 0; i < _num_slip_sys; i++){
      climbstress[i] = 0;
      zi[i] = 0;
      zv[i] = 0;
      Kis[i] = 0;
      Kvs[i] = 0;
      vacdensity_0[i] = 0;
      vel_climb[i] = 1e-25;
      gamma_dot_c[i] = 0;
      for(unsigned int j = 0; j < _num_slip_sys; j++){
        dEtadGc[i][j] = 0;
      }

    }
  }


    //k_D2 = k_D2_0;

  

  // Begin Newton-Raphson iterative loops.
  iNR = 0;

  converged = false;

  do {
    iNR = iNR + 1;
    // converged = true;

    NR_residual (_num_slip_sys, xs0, xm0, _temp[_qp], dt_incr, gamma_dot, gamma_dot_c, F1, F_el, F_p_inv, F_p_inv_0, C, rho_m0, rho_m, rho_i0, rho_i, bstress0, bstress, Nloop0, Nloop, dloop0, dloop, sig, tau, climbstress, tau_eff, s_a, s_t, A, H, interdensity, vacdensity, vacdensity_th, Di, Dv, residual, sse);

    sse_ref = sse;

//       if (sse > 0.0e0) {
//         std::cout << "\n iNR:" << iNR;
//         std::cout << "\n sse:" << sse << "\n";
//       }

    // Begin calculation of the partial derivatives needed for the Newton-Raphson step

    // Calculate derivative of the mobile dislocation density, rho_m, immobile dislocation density, rho_i, backstress, bstress, w.r.t. gamma-dot-beta

    // Variables for derivative with respect to glide component
    std::vector<std::vector<Real>> temp1(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> temp2(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> temp3(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> temp4(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> drhomdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> drhoidgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dbsdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dldrhom(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dsadgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dsa1dgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dsa2dgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dstdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> dNldgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> ddldgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> temp5(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> temp6(_num_slip_sys, std::vector<Real>(_num_slip_sys));


    // Initialize Small Number to various arrays
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        temp1[ia][ib]      = 1e-18;
        temp2[ia][ib]      = 1e-18;
        temp3[ia][ib]      = 1e-18;
        temp4[ia][ib]      = 1e-18;
        drhomdgb[ia][ib]   = 1e-18;
        drhoidgb[ia][ib]   = 1e-18;
        dbsdgb[ia][ib]     = 1e-18;
        dldrhom[ia][ib]    = 1e-18;
        dNldgb[ia][ib]     = 1e-18;
        ddldgb[ia][ib]     = 1e-18;
        temp5[ia][ib]      = 1e-18;
        temp6[ia][ib]      = 1e-18;
        dsa1dgb[ia][ib]    = 1e-18;
        dsa2dgb[ia][ib]    = 1e-18;

        }
    }

    // Update Effective Mean Free Paths
    calc_mfp(_num_slip_sys, rho_m, rho_i, Nloop, dloop, H, d_disl);
      d_disl4 = d_disl;

    //drhomdgb
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      temp2[ia][ia] = temp2[ia][ia] + (k_M/b_mag/d_disl[ia] - k_I/b_mag/d_disl4[ia]
                        - (k_ann*R_c*rho_m[ia]/b_mag))*(sgn(gamma_dot[ia])*dt_incr);

      temp1[ia][ia] = temp1[ia][ia] + 1.0 - ( (0.5*k_M*d_disl[ia]/b_mag)
                                          +   (k_ann*R_c/b_mag)
                                          +   (0.5*k_I*d_disl4[ia]/b_mag) ) * abs(gamma_dot[ia])*dt_incr;
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
      
      Real t1 = (k_I*d_disl4[ia])/(2*b_mag);
      Real t2 = k_D;
      Real t3 = 0.0;// sqrt(Nloop[ia]*dloop[ia])*dloop[ia]/(4*b_mag*sqrt(rho_i[ia]));
      Real t4 = abs(gamma_dot[ia])*dt_incr;
      temp3[ia][ia] = 1 + (-t1 + t2 )*t4;

      Real t5 = (k_I/(b_mag*d_disl4[ia]));
      Real t6 = k_D*rho_i[ia];
      //Real t7 = sqrt(Nloop[ia]*dloop[ia])*sqrt(rho_i[ia])*dloop[ia]/(2*b_mag);
      Real t8 = dt_incr*sgn(gamma_dot[ia]);
      temp4[ia][ia] = (t5 - t6)*t8;
      }
    std::vector<std::vector<Real>> temp3_inv(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    // Method 2: MatrixTools inverse
    MatrixTools::inverse (temp3,temp3_inv);    

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        drhoidgb[ia][ib] = 0.0;
        for(unsigned int ic = 0; ic < _num_slip_sys; ic++) {
          drhoidgb[ia][ib] = drhoidgb[ia][ib] + temp3_inv[ia][ic]*temp4[ic][ib];
        }
      }
    }

    //dNloopdgb
    if(_irad_defects) {
      for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
          Real t1 = sqrt(Nloop[ia]*dloop[ia]);
          Real t2 = sqrt(rho_m[ia]);
          Real t3 = -(dt_incr/(2*b_mag));
          temp5[ia][ia] = (1 + (t2*abs(gamma_dot[ia]*dt_incr))/(4*b_mag*t1));

        //   for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        //     temp6[ia][ib] = t3*(1/(2*t2))*t1*abs(gamma_dot[ia])*drhoidgb[ia][ib];
        
        // }
        temp6[ia][ia] = temp6[ia][ia] + t1 * t2*t3*sgn(gamma_dot[ia]);

        } 

      std::vector<std::vector<Real>> temp5_inv(_num_slip_sys, std::vector<Real>(_num_slip_sys));

      // Method 2: MatrixTools inverse
      MatrixTools::inverse (temp5,temp5_inv);

      for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
        for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
          dNldgb[ia][ib] = 0.0;
          for(unsigned int ic = 0; ic < _num_slip_sys; ic++) {
            dNldgb[ia][ib] = dNldgb[ia][ib] + temp5_inv[ia][ic]*temp6[ic][ib];
          }
        }
      }
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


    // Calculate derivative of the thermal slip resistance, s_t w.r.t. gamma-dot-beta (glide component)
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        dstdgb[ia][ib] = 0.0;
      }
    }

    // Calculate derivative of the athermal slip resistance, s_a w.r.t. gamma-dot-beta. (glide component)
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++){
        Real sum1 = 0.0;
        for(unsigned int ib1 = 0; ib1 < _num_slip_sys; ib1++) {
          //sum1 = sum1 + p0*A[ia][ib]*(rho_m[ib] + rho_i[ib]);
          sum1 = sum1 + A[ia][ib1]*(rho_m[ib1] + rho_i[ib1]);
        }

      Real s_a1{1e-18}; // Accounting for Hardening Due to Dislocation Densities
      s_a1 = k_rho*G*b_mag*sqrt(sum1);        
      
      //dsadgb[ia][ib] = 0.5*k_rho*G*b_mag*(p0*A[ia][ib]*(drhomdgb[ia][ib] + drhoidgb[ia][ib])/sqrt(sum1));
      dsa1dgb[ia][ib] = 0.5*k_rho*G*b_mag*(A[ia][ib]*(drhomdgb[ia][ib] + drhoidgb[ia][ib])/sqrt(sum1));
      

      Real s_a2{1e-18}; // Acounting for Hardening Due to Dislocation Loops
      if(_irad_defects){
      Real t1 = sqrt(Nloop[ia]*dloop[ia]);
      s_a2 = 0.5*G*b_mag*t1;
      dsa2dgb[ia][ib] = 0.5*G*b_mag*(dNldgb[ia][ib] + ddldgb[ia][ib])/t1;
      }

      Real s_a_res{0}; // Resultant s_a due to s_a1 and s_a2
      s_a_res = sqrt(s_a1*s_a1 + s_a2*s_a2);

      Real t2 = s_a1/s_a_res;
      Real t3 = s_a2/s_a_res;

      dsadgb[ia][ib] = t2*dsa1dgb[ia][ib] + t3*dsa2dgb[ia][ib];            
      }
    }

   // (Vik) Derivatives w.r.t. gamma-dot-beta (climb component)
   // drhomdgc = Derivative of Mobile Dislocation wrt to climb slip rate
   std::vector<std::vector<Real>> drhomdgc(_num_slip_sys, std::vector<Real>(_num_slip_sys)); 
   // drhoidgc = Derivative of Immobile Dislocation wrt to climb slip rate
   std::vector<std::vector<Real>> drhoidgc(_num_slip_sys, std::vector<Real>(_num_slip_sys)); 
   // dldgc = Derivative of Mean Free Path wrt to climb slip rate
   std::vector<std::vector<Real>> dldgc(_num_slip_sys, std::vector<Real>(_num_slip_sys)); 
   // dvcdgc = Derivative of Climb Velocity wrt to climb slip rate
   std::vector<std::vector<Real>> dvcdgc(_num_slip_sys, std::vector<Real>(_num_slip_sys)); 
   //dzidgc, dzvdgc = Derivative of Interstitial and Vacancy Capture Efficencies wrt to climb slip rate
   std::vector<std::vector<Real>> dzidgc(_num_slip_sys, std::vector<Real>(_num_slip_sys));    
   std::vector<std::vector<Real>> dzvdgc(_num_slip_sys, std::vector<Real>(_num_slip_sys));    
   //dcvth0dgc = Derivative of Thermal Vacancy Density wrt to climb slip rate
   std::vector<std::vector<Real>> dcvth0dgc(_num_slip_sys, std::vector<Real>(_num_slip_sys));    
   // temp11  and temp12 = Temporary Variables for Calculation
   std::vector<std::vector<Real>> temp11(_num_slip_sys, std::vector<Real>(_num_slip_sys));
   std::vector<std::vector<Real>> temp12(_num_slip_sys, std::vector<Real>(_num_slip_sys));
   std::vector<std::vector<Real>> temp11_inv(_num_slip_sys, std::vector<Real>(_num_slip_sys));

   
   // Initialize All Newly Defined Varibles to Zero
   for(unsigned int ia = 0; ia < _num_slip_sys; ia++){
    for(unsigned int ib = 0; ib <_num_slip_sys; ib++) {
      temp11[ia][ib]      = 1e-18;
      temp12[ia][ib]      = 1e-18;
      drhomdgc[ia][ib]    = 1e-18;
      drhoidgc[ia][ib]    = 1e-18;
      dvcdgc[ia][ib]      = 1e-18;
      dzidgc[ia][ib]      = 1e-18;
      dzvdgc[ia][ib]      = 1e-18;
      dcvth0dgc[ia][ib]   = 1e-18;
      dldgc[ia][ib]       = 1e-18;
    }
   }

  if(_climbmodel){

      // Calculate drhomdgc
      // for (unsigned int ia = 0; ia < _num_slip_sys; ia++ ) {
      //   Real sum1 = 0.0;
      //   Real sum2 = 0.0;
      //   Real sum3 = 0.0;
      //   for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {
      //     sum1 = sum1 + H[ia][ib]*(rho_m[ib] + rho_i[ib]);
      //     sum2 = sum2 + H[ia][ib]*(rho_m[ib]);
      //     sum3 = sum3 + H[ia][ib]*(rho_i[ib]);
      //   }
      //   d_disl2[ia] = pow( (sqrt(sum2)), -1); // Mean Free Path
      //   d_disl3[ia] = 1/sqrt(sum3);
      //   Real d_disl {1/sqrt(sum1)};
      //   Real t1 {k_M*d_disl/2/b_mag};
      //   Real t2_0 {pow(rho_m[ia], -0.5)/(2* b_mag)};
      //   Real t2 = t1*(t2_0)*abs(gamma_dot_c[ia])*dt_incr;
      //   Real t3 = d_disl2[ia]/(2*b_mag)*abs(gamma_dot_c[ia])*dt_incr;

      //   Real t4 = k_M/b_mag/d_disl;
      //   Real t5 = 1/b_mag/d_disl2[ia];
      //   Real t6 = k_D2*rho_i[ia];

      //   temp12[ia][ia] = temp12[ia][ia] + (t4 - t5 + t1*(t5 - t6))*(sgn(gamma_dot_c[ia]) * dt_incr);

      //   temp11[ia][ia] = temp11[ia][ia] + 1 - t1 - t2  + t3;                   

      // }


    // MatrixTools::inverse(temp11, temp11_inv);

    for (unsigned int ia = 0; ia <_num_slip_sys; ia++) {
      // for (unsigned int ib = 0; ib < _num_slip_sys; ib++ ) {
      //   drhomdgc[ia][ib] = 0.0;
      //   for(unsigned int ic = 0; ic < _num_slip_sys; ic++) {
          drhomdgc[ia][ia] = -k_c*rho_m[ia]*sgn(gamma_dot_c[ia])*dt_incr;
      //  }
     // }
    }

    // drhoidgc 
    for (unsigned int ia = 0; ia < _num_slip_sys; ia++){
      // for(unsigned int ib = 0; ib < _num_slip_sys; ib++){
      //   Real t1 = k_c*drhomdgc[ia][ib]*abs(gamma_dot_c[ia]);
      //   drhoidgc[ia][ib] = t1*dt_incr;
      // }
      drhoidgc[ia][ia] = drhoidgc[ia][ia] + (- k_D*rho_i[ia])*sgn(gamma_dot_c[ia])*dt_incr;
    }

    // Calculate dvcdgc
    // Calculate dependent terms 
    for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        dzidgc[ia][ib] = (Zsi/G) * dEtadGc[ia][ib];
        dzvdgc[ia][ib] = 0.0;
        Real x1 = (Ecorr * atomvol) / (B_k * _temp[_qp]);
        dcvth0dgc[ia][ib] = vacdensity_0[ia] * dEtadGc[ia][ib] * x1;
      }
    }

    // dvcdgc
    for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for (unsigned int ib = 0; ib< _num_slip_sys; ib++) {
        dvcdgc[ia][ib] = atomvol/b_mag * ( (dzidgc[ia][ib] * Di * interdensity) -
        - (zv[ia]*Dv*(0- dcvth0dgc[ia][ib])) );
      }
    }

    //dldgc
    
    // for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    //   Real x1 = pow(d_disl2[ia], 2);
    //   Real x2 = pow(rho_m[ia], -0.5);
    //   Real x3 = -0.5*x1*x2;      
    //   for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {     
    //   dldgc[ia][ib] = x3*drhomdgc[ia][ib] ;
    //   }
    // }
  } // end if (climbmodel)   

    // (Vik) Form "A-matrix" of derivatives wrt d_gamma_beta.
    unsigned int totalsys;
    if(_climbmodel){
      totalsys = 2*_num_slip_sys;
    }
    else {
      totalsys = _num_slip_sys;
    }
     
    std::vector<std::vector<Real>> array3(totalsys, std::vector<Real>(totalsys));

    // Initialize all components to zero
    for (unsigned int ia = 0; ia< (totalsys) ; ia++){
      for (unsigned int ib = 0; ib<(totalsys) ; ib++){
        array3[ia][ib] = 0;
      }
    }

    // Slip Systems
    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      tau_eff[ia] = tau[ia] - bstress[ia];

     
      if(_sliplaw == 2) {  // * For a Slip Model which contains both athermal term s-a & thermal term s-t
        for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
          array3[ia][ib] = -gamma_dot_g[ia]*(p*q*(delF0/B_k/_temp[_qp])*power(1.e0 - power((abs(tau_eff[ia])
              - s_a[ia])/s_t[ia], p), q - 1)*power((abs(tau_eff[ia]) - s_a[ia])/s_t[ia], p - 1))*
            (((dTaudGd[ia][ib] - dbsdgb[ia][ib])*sgn(tau_eff[ia]) - dsadgb[ia][ib])/s_t[ia] - (abs(tau_eff[ia]) - s_a[ia])*dstdgb[ia][ib]/(s_t[ia]*s_t[ia]));
        }
      }

      else if(_sliplaw ==1) {

      // *  For Viscoplastic Slip Model Containing only s-a term 
        for(unsigned int ib = 0; ib < _num_slip_sys; ib++){
          Real x1 = -delF0/(B_k*_temp[_qp]);
          Real x2 = abs(tau_eff[ia])/s_a[ia];
          if (x2 == 0) { x2 = 1e-18;}
          Real x3 = pow(x2, p);
          Real x4 = pow(x2, p-1);
          Real x5 = q* pow(1-x3, q-1);
          Real x6 = -p*x4;
          Real x7 = s_a[ia]*dTaudGd[ia][ib]*sgn(tau_eff[ia]) - tau_eff[ia]*dsadgb[ia][ib];
          Real x8 = pow(s_a[ia],2);
          Real x9 = x7/x8;
          Real x10 = drhomdgb[ia][ib]/rho_m[ia]*0;
          array3[ia][ib] = -gamma_dot_g[ia]*(x1*x5*x6*x9 - x10);
        }
      }
        array3[ia][ia] = array3[ia][ia] + 1;
      }
    
    if(_climbmodel){
      //Climb Systems
      for(unsigned int ia = _num_slip_sys; ia< (totalsys); ia++){
        int ij = ia - _num_slip_sys;
        for(unsigned int ib= _num_slip_sys; ib< (totalsys); ib++){
          int ik = ib - _num_slip_sys;
          array3[ia][ib] = - (   (drhomdgc[ij][ik]*b_mag*vel_climb[ij])
                              + (rho_m[ij]*b_mag*dvcdgc[ij][ik]) );
        }
        array3[ia][ia] = array3[ia][ia] + 1;
      }
    }

    // Calculate the gradient of sse wrt gamma_dot(). Will be used later to ensure that line search is in the correct direction
    std::vector<Real> gradient(totalsys);

    for(unsigned int j = 0; j < (totalsys); j++) {
      gradient[j] = 0.0;
      for(unsigned int i = 0; i < (totalsys); i++) {
        gradient[j] = gradient[j] + residual[i]*array3[i][j];
      }
      gradient[j] = 2.0*gradient[j];
    }

    // Solve for increment of gamma_dot
    std::vector<Real> d_gamma_dot(totalsys);

    std::vector<std::vector<Real>> array3_inv(totalsys, std::vector<Real>(totalsys));

    MatrixTools::inverse(array3,array3_inv);

    for(unsigned int j = 0; j < (totalsys); j++) {
      d_gamma_dot[j] = 0.0;
      for(unsigned int i = 0; i < (totalsys); i++) {
        d_gamma_dot[j] = d_gamma_dot[j] - array3_inv[j][i]*residual[i];
      }
    }

    // Check to make sure that N-R step leads 'down hill' the sse surface
    Real sum1 = 0.0;
    for(unsigned int i = 0; i < (totalsys); i++) {
      sum1 = sum1 - gradient[i]*d_gamma_dot[i];
    }

    if(sum1 > 0.0) {
      for(unsigned int i = 0; i < (totalsys); i++) {
        d_gamma_dot[i] = -1*d_gamma_dot[i];
      }
    }

    // Multiply step size by two 'cause next loop will divide it by 2
    for(unsigned int k = 0; k < (totalsys); k++) {
      d_gamma_dot[k] = d_gamma_dot[k]*2;
    }

    // Begin line search
    improved = false;

    N_ctr = 0;

    do {
      sse_old = sse;

      // Divide step size by 2
      for (unsigned int k = 0; k < (_num_slip_sys); k++) {
        d_gamma_dot[k] = d_gamma_dot[k]/2.0;
        gamma_try_g[k] = gamma_dot_g[k] + d_gamma_dot[k];
      }

      if(_climbmodel){
        for (unsigned int k = 0; k < (_num_slip_sys); k++){
          int j = k + _num_slip_sys;
          d_gamma_dot[j] = d_gamma_dot[j]/2.0;
          gamma_try_c[k] = gamma_dot_c[k] + d_gamma_dot[j];
        }
      }
      else {
        for (unsigned int k = 0; k < (_num_slip_sys); k++){
          gamma_try_c[k] = 0;
        }
      }        
      
      NR_residual (_num_slip_sys, xs0, xm0, _temp[_qp], dt_incr, gamma_try_g,   gamma_try_c, F1, F_el, F_p_inv, F_p_inv_0, C, rho_m0, rho_m, rho_i0,     rho_i, bstress0, bstress, Nloop0, Nloop, dloop0, dloop, sig, tau,    climbstress, tau_eff, s_a, s_t, A, H, interdensity, vacdensity,     vacdensity_th, Di, Dv, residual, sse);

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

    if(_climbmodel){
      for (unsigned int k = 0; k < _num_slip_sys; k++) {
        int j = k + _num_slip_sys;
        gamma_dot_c[k] = gamma_dot_c[k] + d_gamma_dot[j]*2.0;
      }      
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

      //std::cout << "\n Total sub-increments:" << N_incr_total << "\n";
      break;
    }

  } while((iNR < 10000) && (converged == false));


    if (iNR == 10000) {

      std::cout << "\n SSE: \n"; 
      std::cout << sse;
      std::cout << "\n SSE_OLD ";
      std::cout << sse_old;
      std::cout << "\n SSE_REF: ";
      std::cout << sse_ref;

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


      throw MooseException("Crystal plasticity Newton Raphson did not converge\n for grain no no " + std::to_string(_grainid) + 
      "\n  Qp no: " + std::to_string(_qp) +
      "\n Increment no:  " + std::to_string(N_incr_total)); 
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

    throw MooseException("Crystal plasticity time step subincrement did not converge for element no " + std::to_string(_grainid) + " ");
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

  // Calculate Mean Free Path for Further Calculations
  std::vector<Real> d_disl(_num_slip_sys); // dislocation spacing
  calc_mfp(_num_slip_sys, rho_m, rho_i, Nloop, dloop, H, d_disl);



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
          // Slip Systems

          if(_sliplaw == 2) {// * For Slip Law Containing s-a & s-t
            for(unsigned int n = 0; n < _num_slip_sys; n++) {
              tau_eff[n] = tau[n] - bstress[n];

              ddpdsig[i][j][k][l] = 
              ddpdsig[i][j][k][l] + 
              (xs[i][n]*xm[j][n] + xm[i][n]*xs[j][n])*(((xs[k][n]*xm[l][n] + xm[k][n]*xs[l][n])*
              (gamma_dot_g[n]*(p*q*(delF0/B_k/_temp[_qp])*power(1.0 - power((abs(tau_eff[n]) - s_a[n])/s_t[n], p), q - 1.0)*
              power((abs(tau_eff[n]) - s_a[n])/s_t[n], p - 1.0))*sgn(tau_eff[n])/s_t[n])));
            }
          }
          else if (_sliplaw == 1){ // * For ViscoPlastic Slip Law
            for (unsigned int n = 0; n < _num_slip_sys; n++){
              tau_eff[n] = tau[n] - bstress[n];
              Real x1 = (delF0/(B_k*_temp[_qp]));
              Real x2 = abs(tau[n])/s_a[n]; if (x2 == 0) { x2 = 1e-18;}
              Real x3 = pow(x2, p);
              Real x4 = pow(x2, p-1);
              Real x5 = pow(1 - x3, q-1);

              ddpdsig[i][j][k][l] = 
              ddpdsig[i][j][k][l] +
              (xs[i][n]*xm[j][n] + xm[i][n]*xs[j][n])*
              gamma_dot_g[n]*(x1*q*x5*p*x4)*sgn(tau_eff[n])/s_a[n];
              (xs[k][n]*xm[l][n] + xm[k][n]*xs[l][n]);
            }
          }

          // Climb Systems
          if(_climbmodel) {
          for (unsigned int n = _num_slip_sys; n< (2*_num_slip_sys); n++) {
            int m = n - _num_slip_sys;
            Real x1 = (Ecorr * atomvol)/(B_k*_temp[_qp]);
            Real x2 = Zsi/G;
            ddpdsig[i][j][k][l] = 
            ddpdsig[i][j][k][l] +
            (xs[i][m]*xs[j][m] + xs[i][m]*xs[j][m])  *  (xs[k][m]*xs[l][m]+xs[k][m]*xs[l][m]) * 
            (rho_m[m]*b_mag)* ((atomvol/b_mag)*
            (x2*Di*interdensity + zv[m]*Dv*vacdensity_0[m]*x1));
          }
          }
        }
      }
    }
  }

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
  Real gamma_dot_climb_avg;
  Real climbfrac1 = 0, climbfrac2=0;
  for (unsigned int i = 0; i < _num_slip_sys; i++) {
    rho_m_avg = rho_m_avg + rho_m[i];
    rho_i_avg = rho_i_avg + rho_i[i];
    gamma_dot_avg = gamma_dot_avg + abs(gamma_dot[i]);
    gamma_dot_climb_avg = gamma_dot_climb_avg + abs(gamma_dot_c[i]);
  }
  rho_m_avg = rho_m_avg/_num_slip_sys;
  rho_i_avg = rho_i_avg/_num_slip_sys;
  climbfrac1 = (gamma_dot_climb_avg)/ (gamma_dot_avg + gamma_dot_climb_avg);

  gamma_dot_avg = gamma_dot_avg/_num_slip_sys;
  gamma_dot_climb_avg = gamma_dot_climb_avg/_num_slip_sys;

 

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

  if(_irad_defects) {
    // Store Nloop values SDV 50 + _num_slip_sys*4
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      _state_var[_qp][n] = Nloop[i];
    }

    // Store dloop values SDV 50 + _num_slip_sys*5
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      _state_var[_qp][n] = dloop[i];
    }
  }
  
  // Read Dislocation Glide Rate SDV50 + _num_slip_sys*6
  for(unsigned int i = 0; i < _num_slip_sys; i++) {
    n = n+1;
    _state_var[_qp][n] = gamma_dot[i];
  }

  if(_climbmodel){
    // Read Dislocation Climb Rates SDV50 + _num_slip_sys*7
    for(unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n+1;
      _state_var[_qp][n] = gamma_dot_c[i];
    } 

    // Store average climb rates SDV50 + _num_slip_sys*7 + 1
    n = n+1;
    _state_var[_qp][n] = gamma_dot_climb_avg;

    // Store Effect of Climb rates SDV50 + _num_slip_sys*7 + 3
    n = n+1;
    _state_var[_qp][n] = climbfrac1;

  // Climbfrac calculations by method 2
  Real Lp_g[3][3], Lp_c[3][3];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Lp_g[i][j] = 0.0;
      Lp_c[i][j] = 0.0;
      for (unsigned int k = 0; k < _num_slip_sys; k++) {
        Lp_g[i][j] = Lp_g[i][j] +   xs0[i][k]*xm0[j][k]* gamma_dot[k];
        Lp_c[i][j] = Lp_c[i][j] +   xs0[i][k]*xs0[j][k]* gamma_dot_c[k];
      }
    }
  }

  Real Ep_g = 0, Ep_c = 0, sg= 0, sc = 0;
  for (unsigned int i = 0; i<3; i++) {
    for (unsigned int j = 0; j<3; j++) {
      sg = sg + Lp_g[i][j]*Lp_g[i][j];
      sc = sc + Lp_c[i][j]*Lp_c[i][j];
    }
  }
  Ep_g = sqrt((2./3.) * sg);    
  Ep_c = sqrt((2./3.) * sc);
  climbfrac2 = Ep_c/(Ep_c + Ep_g);       
  
  n = n+1;
  _state_var[_qp][n] = climbfrac2;

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
void ThermalIrradiationCPUpdate::NR_residual (
  unsigned int num_slip_sys, // Number of Slip Systems
  std::vector<std::vector<Real>> &xs0, // Slip Direction Matrix
  std::vector<std::vector<Real>> &xm0, // Slip Normal Matrix
  Real temp,                           // Temperature
  Real dt,                             // Time Increment
  std::vector<Real> gdt,               // gammadot_glide Trial Component
  std::vector<Real> gdc,               // gammadot_Climb Trial Component 
  RankTwoTensor F1,                    // Final Deformation Gradient of Current SubIncrement
  RankTwoTensor &F_el,                 // Elastic Part of F
  RankTwoTensor &F_p_inv,              // Inverse of Plastic Part of F at end of step
  RankTwoTensor F_p_inv_0,             // Inverse of Plastic Part of F at beginging of step
  Real C[3][3][3][3],                  // Elastic Stiffness Tensor
  std::vector<Real> &rho_m0,           // Mobile Dislocation Density at the Begining of Step
  std::vector<Real> &rho_m,            // Mobile Dislocation Density at the End of Step
  std::vector<Real> &rho_i0,           // Immobile Dislocation Density at the Begining of Step
  std::vector<Real> &rho_i,            // Immobile Dislocation Density at the End of Step
  std::vector<Real> &bstress0,         // BackStress at the Begining of Step
  std::vector<Real> &bstress,          // BackStress at the End of Step
  std::vector<Real> &Nloop0,           // Density of Dislocation Loops at the Begining of Step
  std::vector<Real> &Nloop,            // Density of Dislocation Loops at the End of Step
  std::vector<Real> &dloop0,           // Diameter of Dislocation Loops at the Begining of Step 
  std::vector<Real> &dloop,            // Diameter of Dislocation Loops at the End of Step 
  RankTwoTensor &sig,                  // Cauchy Stress
  std::vector<Real> &tau,              // Resolved Shear Stress on Each Slip System
  std::vector<Real> &climbstress,      // Resolved Climb Stress on Each Slip System
  std::vector<Real> &tau_eff,          // Effective Shear Stress on Each Slip System
  std::vector<Real> &s_a,              // athermal slip resistance of Each Slip System
  std::vector<Real> &s_t,              // thermal slip resistance of Each Slip System
  std::vector<std::vector<Real>> A,    // Slip Interaction Coefficient Matrix
  std::vector<std::vector<Real>> H,    // Hardening Interaction Coefficient Matrix
  Real interdensity,                   // Density of Interstitials
  Real vacdensity,                     // Density of Vacancies
  Real vacdensity_th,                  // Density of Vacancies which are in Thermal Equilibrium with material
  Real Di,                             // Diffusivity of Interstitials 
  Real Dv,                             // Diffusivity of Vacancies 
  std::vector<Real> &residual,         // Difference Between the Trial Slip Value and Slip Values calculated by New Fp
  Real &sse){                          // RMS of Residual



  Real xL_p_inter[3][3];
  std::vector<Real> shr_g(num_slip_sys);
  std::vector<Real> shr(num_slip_sys);
  std::vector<Real> shr_c(num_slip_sys);

  // RankTwoTensor F_el;
  RankTwoTensor E_el;
  RankTwoTensor Spk2;

  std::vector<std::vector<Real>> xs(3, std::vector<Real>(num_slip_sys));
  std::vector<std::vector<Real>> xm(3, std::vector<Real>(num_slip_sys));

  // (Vik) Calculate the plastic part of the velocity gradient in the intermediate configuration - Contribution of Glide & Climb Slip
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      xL_p_inter[i][j] = 0.0;
      for (unsigned int k = 0; k < num_slip_sys; k++) {
        xL_p_inter[i][j] = xL_p_inter[i][j] + 
        xs0[i][k]*xm0[j][k]* gdt[k] + 
        xs0[i][k]*xs0[j][k]* gdc[k];
      }
    }
  }

  // (Vik) Calculate the Strain Rate from the velocity gradient
  Real xL_p_inter_s[3][3];
  Real edot{0.e0};
  for (unsigned int i = 0; i <3; i++){
    for (unsigned int j = 0; j <3 ; j++){
      xL_p_inter_s[i][j] = 0.5*(xL_p_inter[i][j] + xL_p_inter[j][i]);
    }
  }

  for (unsigned int i = 0; i <3; i++){
    for (unsigned int j = 0; j <3 ; j++){
      edot = edot + xL_p_inter_s[i][j]*xL_p_inter_s[i][j];
    }
  } 
  edot = sqrt(edot); 
  edot = sqrt(0.667e0)*edot;
  if (abs(edot) < 1e-18){
        edot = 1e-18;
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

  // (Vik) Calculate resolved climb stress for each slip system
  if(_climbmodel){
    for (unsigned int k = 0; k < num_slip_sys; k++) {
      climbstress[k] = 0.0e0;
      for (unsigned int i = 0; i< 3; i++) {
        for(unsigned int j = 0; j < 3; j++) {
          climbstress[k] = climbstress[k] + xs[i][k]*xs[j][k]*sig(i,j);
        }
      }
    }
  }

  // (Vik) Calculate dislocation densities, rho_m, rho_i
  std::vector<Real> d_disl(num_slip_sys);
  std::vector<Real> d_disl4(num_slip_sys);
  calc_mfp(_num_slip_sys, rho_m0, rho_i0, Nloop0, dloop0, H, d_disl);

  d_disl4 = d_disl;
  // for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
  //   Real sum1 = 1e-18;
  //   Real sum2 = 1e-18;
  //   for(unsigned int ib = 0; ib < num_slip_sys; ib++) {
  //     sum1 = sum1 + H[ia][ib]*(rho_m0[ib] + rho_i0[ib]);
  //     sum2 = sum2 + H[ia][ib]*(rho_m0[ib]);
  //   }
  //   d_disl[ia] = 1.0e0/sqrt(sum1);
  // }

  std::vector<Real> drhomdt(num_slip_sys);
  std::vector<Real> drhoidt(num_slip_sys);
  std::vector<Real> dbstressdt(num_slip_sys);
  std::vector<Real> dNloopdt(num_slip_sys);
  std::vector<Real> ddloopdt(num_slip_sys);

  // Effective Mean Free Path for Dislocation Trap
  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    // Real sum1 = 1e-18;
    // for(unsigned int ib = 0; ib < num_slip_sys; ib++) {
    //   sum1 = sum1 + H[ia][ib]*(rho_m0[ib] + rho_i0[ib]);
    // }
    // Real sum2{1e-18};
    // Real sum3{1e-18};
    // if(_irad_defects) {
    // sum2 = beta_loop*sqrt(Nloop[ia]*dloop[ia]);        
    // sum3 = 1.0/( (1.0/sum1) + (1.0/sum2));
    // }
    // else{
    //   sum3 = sum1;
    // }

    
    //Real k_c = 1e4;

    drhomdt[ia] = ((k_M/b_mag/d_disl[ia]                          // Dislocation Multiplication Component
                 - k_ann*R_c*rho_m0[ia]/b_mag                     // Mutual Annhilation
                 - k_I/b_mag/d_disl4[ia])*abs(gdt[ia])           // Trapped Dislocations
                 - k_c * abs(gdc[ia]) * rho_m0[ia]);              // Climb Component Annhilation   
                 //+ k_M/b_mag/d_disl[ia]*abs(gdc[ia]);             // Dislocation Multiplication by Climb


    
    if(_climbmodel){
      k_D2 = k_D2_0*pow((edot0/edot), 1/n0);
      k_D = k_D2;
    }
    else{
      k_D2 = 0;
    }

    drhoidt[ia] = (k_I/b_mag/d_disl4[ia]                       // Trapped Dislocation - Glide Component
                  - k_D*rho_i0[ia])*abs(gdt[ia])                // Dynamic Recovery of Dislocations
                  - k_D2 * rho_i0[ia]* abs(gdc[ia]);            // Relaxation of Immobile Dislocations due to Climb
                //  +  k_c * abs(gdc[ia]) * rho_m0[ia];            // Climb Component Addition


    if (_irad_defects) {

    dNloopdt[ia] = k_c*rho_m0[ia]*abs(gdc[ia])/(pi*dloop0[ia]) // Increase in Number Density Dislocation Loops Due to Climb
      - (1.2/(2*b_mag))*sqrt(Nloop0[ia]*dloop0[ia])*sqrt(rho_m0[ia])*abs(gdt[ia]); // Decrease in Number Density of Dislocation Loops due to Annhilation
    
    Real zi1 = zi0 + (Zsi/G)*climbstress[ia];
    ddloopdt[ia] = (2.0/pi)*(Ai0_by_i0)*(zi1* Di* rho_m0[ia]* interdensity)/(Nloop0[ia]*dloop0[ia]); // Increase in Size of Dislocation Loops due to Climb
    }

    else {
       dNloopdt[ia] = 0;
       ddloopdt[ia] = 0;
    }
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

    Nloop[ia] = Nloop0[ia] + dNloopdt[ia]*dt;
    if(Nloop[ia] < 1e-6) {
      Nloop[ia] = 1e-6;
    }

    dloop[ia] = dloop0[ia] + ddloopdt[ia]*dt;
    if (dloop[ia] < 1e-6) {
      dloop[ia] = 1e-6;
    }
  }

  // Update the Mean Free Paths
  calc_mfp(_num_slip_sys, rho_m, rho_i, Nloop, dloop, H, d_disl);
  d_disl4 = d_disl;
  // for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
  //   Real sum1 = 1e-18;
  //   for(unsigned int ib = 0; ib < num_slip_sys; ib++) {
  //     sum1 = sum1 + H[ia][ib]*(rho_m0[ib] + rho_i0[ib]);
  //   }
  //   d_disl[ia] = 1.0e0/sqrt(sum1);
  // }  


  // Calculate the back stress for each slip system.
  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    dbstressdt[ia] = (k_bs1*k_rho*G*b_mag*sqrt(rho_m[ia] + rho_i[ia])*sgn(tau[ia] - bstress0[ia]) - k_bs2*bstress0[ia])*abs(gdt[ia]);

    bstress[ia] = bstress0[ia] + dbstressdt[ia]*dt;
  }

  // Calculate the athermal slip resistance for each slip system.
  calc_crss(_num_slip_sys, rho_m, rho_i, Nloop, dloop, s_a, A);


  // Calculate thermal slip resistance for each slip system.
  for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
    s_t[ia] = frictional_stress;
  }

  // Calculate new slip rate for each slip system.
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    tau_eff[k] = tau[k] - bstress[k];
  }

  calc_sliprate(_num_slip_sys, tau_eff, s_a, s_t, shr_g, shr);

  // (Vik) Variables for Calculation
  std::vector<Real> zi(num_slip_sys); // Finite Stress Interstitial Capture Efficency
  std::vector<Real> zv(num_slip_sys); // Finite Stress Vacancy Capture Efficency
  std::vector<Real> Kis(num_slip_sys); // Insterstial sink (Dislocation) reaction rate coefficient
  std::vector<Real> Kvs(num_slip_sys); // Vacancy sink (Dislocation) Reaction rate coefficient
  std::vector<Real> vel_climb(num_slip_sys); // Climb Velocity
  std::vector<Real> vacdensity_0(num_slip_sys); // Density of Vacancies on each slip system
  Real dummy1, dummy2; // Dummy Variables used for calculation
  

  // (Vik) Calculate interstial and vacancy capture efficiencies and finite stress thermal equilibrium density of vacancies
  if(_climbmodel){

    calc_climbrate(_num_slip_sys, Di, Dv, vacdensity_th, interdensity, vacdensity,      climbstress, zi, zv, Kis, Kvs, vacdensity_0, vel_climb, rho_m, rho_i, shr_c);
  
  }

  // (Vik) Calculate residual values
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    residual[k] = gdt[k] - shr_g[k];
    if (residual[k] >= 1.0e100) {
      residual[k] = 1.0e100;
    }
    if (residual[k] <= -1.0e100) {
      residual[k] = -1.0e100;
    }
  
  }

  if(_climbmodel){
    for(unsigned int k = 0; k < num_slip_sys; k++){
      int kj = num_slip_sys + k;
      residual[kj] = gdc[k] -shr_c[k];
      if (residual[kj] >= 1.0e100) {
        residual[kj] = 1.0e100;
      }
      if (residual[kj] <= -1.0e100) {
        residual[kj] = -1.0e100;
      }        

    }
  }

  // (Vik) Calculate root mean square error that is to be minimized to zero
  sse = 0.0e0;
  Real gamma_dot_max = abs(gdt[0]);
  if(_climbmodel){
    for(unsigned int k = 0; k < num_slip_sys; k++) {
      int j = num_slip_sys + k;
      sse = sse + 
      abs(gdt[k])*(residual[k]*residual[k]) +
      abs(gdc[k])*(residual[j]*residual[j]);

      if(abs(gdt[k]) > gamma_dot_max) {
        gamma_dot_max = abs(gdt[k]);
      }

      if(abs(gdc[k] > gamma_dot_max)){
        gamma_dot_max = abs(gdc[k]);
      }

    }
    if (gamma_dot_max > 1.e-18) {
      int x = 2*num_slip_sys;
      sse = sqrt(sse/gamma_dot_max);
      sse = sse/x;
    }
  }
  else {
    for(unsigned int k = 0; k < num_slip_sys; k++) {
      sse = sse + 
      abs(gdt[k])*(residual[k]*residual[k]); 

      if(abs(gdt[k]) > gamma_dot_max) {
        gamma_dot_max = abs(gdt[k]);
      }
    }
    if (gamma_dot_max > 1.e-18) {
      int x = num_slip_sys;
      sse = sqrt(sse/gamma_dot_max);
      sse = sse/x;
    }
    
  }
  // std::cout<< "\n _qp:" << _qp << "\t sse:" << sse;

  
} // end NR_residual

void ThermalIrradiationCPUpdate::readPropsFile() {
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

void ThermalIrradiationCPUpdate::assignProperties(){

  // Elastic constants
  C11P = _properties[_qp][0];
  C11_perK = _properties[_qp][1];
  C12P = _properties[_qp][2];
  C12_perK = _properties[_qp][3];
  C44P = _properties[_qp][4];
  C44_perK = _properties[_qp][5];
  GP = _properties[_qp][6];
  G_perK = _properties[_qp][7];

  C11 = C11P - C11_perK*_temp[_qp];
  C12 = C12P - C12_perK*_temp[_qp];
  C44 = C44P - C44_perK*_temp[_qp];
  G = GP - G_perK*_temp[_qp];

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
  Di0 = _properties[_qp][30]; // Interstitial Diffusion Constant
  Dv0 = _properties[_qp][31]; // Vacancy Diffusion Constant
  Em_i = _properties[_qp][32]; // Migration Enthalpy of Vacancy
  Em_v = _properties[_qp][33]; // Migration Enthalpy of Interstitials
  G0 = _properties[_qp][34]; // Gibbs Free Energy of Vacancy Formation
  zi0 = _properties[_qp][35]; //Zero Stress Interstitial Capture Efficency
  zv0 = _properties[_qp][36]; //Zero Stress Vacacny Capture Efficency
  Zsi = _properties[_qp][37]; //ClimbStress Propornality Constant of Interstitial Capture Efficency
  K0 = _properties[_qp][38]; // Point Defect Generation Rates
  Kiv = _properties[_qp][39]; // Point Defect Recombination Rates
  atomvol = _properties[_qp][40]; // Atomic Volume
  lc = _properties[_qp][41]; // Average Climb Distance
  Nv = _properties[_qp][42]; // No of Atoms Per Unit Volume 
  Ecorr = _properties[_qp][43]; // Correction Factor for Energy (Current Unit Energy/1J)
  k_D2_0 = _properties[_qp][44]; // Immobile Dislocations Dynamic Recovery Coefficient
  edot0 = _properties[_qp][45]; //Reference Strain Rate for Dynamic Recovery
  n0 = _properties[_qp][46]; // Reference Strain Rate exponent for Dynamic Recovery
  k_c = _properties[_qp][47]; // Mean Free Path of Precipitates
  Nloop_zero = _properties[_qp][48]; // Initial Density of Dislocation Loops
  dloop_zero = _properties[_qp][49]; // Initial Diameter of Dislocation Loops
  beta_loop = _properties[_qp][50]; //51 Trapping Coefficient of Dislocation Loops
  Ai0_by_i0 = _properties[_qp][51]; //Increase in Area of Dislocation Loop by Addition of a single interstial
  mfpSC =  _properties[_qp][52]; // Hardening of Solute Clusters
  alpha_loop = _properties[_qp][53]; // Used for Calculation of Hardening due to Irradiated Loops
  enthalpy_const2 = _properties[_qp][54]; // Multiplication constant for enthalpy term for glide at Low Temp
  p2 = _properties[_qp][55];  // Shape parameter (dislocation glide) at Low Temp
  q2 = _properties[_qp][56];  // Shape parameter (dislocation glide) at Low Temp
  transition_temp = _properties[_qp][57]; // Temperature at which the change in enthalpy term occurs

  // Initialize Physical Constants
  B_k = 1.3806503e-23;  // Boltzmann constant
  freq = 1e13;  // Debye frequency

  // adjusted to change units of product with stress to SI units
  act_vol = (pow(b_mag,3))*1.0e-3;


  if(_deltaH_eV){
    delF0 = enthalpy_const*1.6e-19;

  }
  else {
    // delF0 scaled to SI units
    delF0 = (enthalpy_const*G*1.0e6)*(pow(b_mag,3))*1.0e-9;
  }

}

// normalize to unit vector
void ThermalIrradiationCPUpdate::normalize_vector(Real *x, Real *y, Real *z) {
  Real length = sqrt(*x * *x + *y * *y + *z * *z);

  *x = *x/length;
  *y = *y/length;
  *z = *z/length;
}

// Rotate 4th order stiffness tensor
void ThermalIrradiationCPUpdate::rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]) {
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
void ThermalIrradiationCPUpdate::forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]) {
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
void ThermalIrradiationCPUpdate::Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]) {
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

void ThermalIrradiationCPUpdate::aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]) {

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

void ThermalIrradiationCPUpdate::aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]) {

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

void ThermalIrradiationCPUpdate::aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product[i][j] = 0.0;
      for (unsigned int k = 0; k < 3; k++) {
        product[i][j] = product[i][j] + a[i][k]*b[k][j];
      }
    }
  }
}

Real ThermalIrradiationCPUpdate::aa_dot_dot_bb(Real a[3][3], Real b[3][3]) {

  Real product = 0.0e0;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product = product + a[i][j]*b[i][j];
    }
  }
  return product;
}

// calculate Euler angles in Bunge notation
void ThermalIrradiationCPUpdate::bunge_angles(Real (&array1)[3][3], Real (&psi0)[3]) {

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
Real ThermalIrradiationCPUpdate::power(Real x, Real y) {

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
Real ThermalIrradiationCPUpdate::sgn(Real x) {

  Real ans;
  ans = 1.0;

  if (x < 0.0e0) {
    ans = -1.0;
  }

  return ans;
}

// return maximum of two scalars
Real ThermalIrradiationCPUpdate::max_val(Real a, Real b)
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

// Read IntMatFile Data Functions
int ThermalIrradiationCPUpdate::numrows(std::string filename){
    int numlines{0};
    
    
    std::ifstream file(filename);
    if (!file.is_open())
    {
       mooseError("Error: Unable to open IntMat file ");
        return -1;
    }

    std::string line;
    while(std::getline(file,line))
    {
        if(!line.empty()){
            numlines++;
        }
        
    }
        
    file.close();    

    return numlines;
}

int ThermalIrradiationCPUpdate::numcols(std::string filename){
    int columncount{0};
    std::ifstream file(filename);
    if (!file.is_open()){
        mooseError("Error to Open IntMat File" );
    }
    std::string line;
    if (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
        while( ss >> token){
            columncount++;
        }

    }
        file.close();
        return columncount;  
}

void ThermalIrradiationCPUpdate::readfile(std::string filename){

    // Read in Euler Angles from IntMat File
    MooseUtils::checkFileReadable(filename);
    std::ifstream fdata;
    fdata.open(filename.c_str());

    // Find Out No of Rows and Cols in the File
    int rows = numrows(filename);
    int cols = numcols(filename);
    if(rows !=cols){
        mooseError("NumRows != NumCols in IntMat File");
    }

    // Resize the IntMatData Matrix
    sintmat.resize(rows);
    for (auto & row :sintmat){
        row.resize(cols);
    }

    int ia{0};
    std::string line;
    while(getline(fdata, line)){
        std::istringstream ss{line};
        for (int i = 0; i < cols; i++){
            ss >> sintmat[ia][i];
        }
        ia++;
    }  

}

void ThermalIrradiationCPUpdate::calc_crss(
    unsigned int num_slip_sys,
    std::vector<Real> &rho_m,
    std::vector<Real> &rho_i,
    std::vector<Real> &Nloop,
    std::vector<Real> &dloop,
    std::vector<Real> &s_a,
    std::vector<std::vector<Real>> A){
        
    for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      Real sum1 = 0.0;
          
      for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        //  sum1 = sum1 + p0*A[ia][ib]*(rho_m[ib] + rho_i[ib]);
        sum1 = sum1 + A[ia][ib]*(rho_m[ib] + rho_i[ib]);
      }
      Real s_a1{0}; // Accounting for Hardening Due to Dislocation Densities
      s_a1 = k_rho * G * b_mag*sqrt(sum1);

      Real s_a2{0}; // Acounting for Hardening Due to Dislocation Loops
      Real s_a3{0}; // Acounting for Hardening Due to Solute Clusters
      if(_irad_defects)
        {
          Real sum2{0};
          for (unsigned int ia = 0; ia < _num_slip_sys ; ia++){
            sum2 = sum2 + Nloop[ia]*dloop[ia];
          }
          s_a2 = alpha_loop* G * b_mag*sqrt(sum2);
          s_a3 = alpha_loop* G*b_mag/mfpSC;;
        }

        Real s_a_res{0}; // Resultant s_a due to s_a1 and s_a2
        s_a_res = sqrt(s_a1*s_a1 + s_a2*s_a2  + s_a3*s_a3 );

        s_a[ia] = tau0 + hp_coeff/sqrt(grain_size) + s_a_res;
      }

}

void ThermalIrradiationCPUpdate::calc_sliprate(
    unsigned int num_slip_sys,
    std::vector<Real> &tau_eff,
    std::vector<Real> &s_a,
    std::vector<Real> &s_t,
    std::vector<Real> &gamma_dot_g,
    std::vector<Real> &gamma_dot){

    if(_sliplaw == 2) { //* ViscoPlastic Slip Law containing terms s_a and s_t

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
    }
    else if (_sliplaw == 1) { // * Viscoplastic Slip Law

      for(unsigned int k = 0; k < _num_slip_sys; k++) {
        Real refgdg = gammadot0g;
        Real x1 = delF0/(B_k*_temp[_qp]);
        Real x2 = abs(tau_eff[k])/s_a[k];
        Real x3 = pow(x2, p);
        Real x4 = pow(1 -x3, q);

        gamma_dot_g[k] = refgdg*exp(-x1*x4)*sgn(tau_eff[k]);
        gamma_dot[k] = gamma_dot_g[k];
      }
    }

    else {
      mooseError("Invalid Parameter : _sliplaw");
    }

}


void ThermalIrradiationCPUpdate::calc_climbrate(
    unsigned int num_slip_sys,
    Real Di,
    Real Dv,
    Real vacdensity_th,
    Real &interdensity,
    Real &vacdensity,
    std::vector<Real> &climbstress,
    std::vector<Real> &zi,
    std::vector<Real> &zv,
    std::vector<Real> &Kis,
    std::vector<Real> &Kvs,
    std::vector<Real> &vacdensity_0,
    std::vector<Real> &vel_climb,
    std::vector<Real> &rho_m,
    std::vector<Real> &rho_i,
    std::vector<Real> &gamma_dot_c){
        
    // Declare all parameters
    Real temp1{0}, temp2{0}, temp3{0};
    Real dummy1, dummy2;
    Real lg;
           
    // Calculate zi, zv, Kis, Kvs, vacdensity_0
    for (unsigned int i = 0; i < _num_slip_sys ; i++) {
      zi[i] = zi0 + (Zsi/G) * climbstress[i];
      zv[i] = zv0;
      Kis[i] = abs(zi[i]) * Di;
      Kvs[i] = abs(zv[i]) * Dv;
      vacdensity_0[i] = vacdensity_th * exp((climbstress[i] * atomvol * Ecorr)/(B_k*_temp[_qp]));
    }

    // Calculate the Saturated Point Defect Density

    for (unsigned int k = 0; k < _num_slip_sys ; k++) {
      temp1 = temp1 + zi[k]*Di*(rho_m[k] + rho_i[k]);
      temp2 = temp2 + zv[k]*Dv*(rho_m[k] + rho_i[k]);
      temp3 = temp3 + zv[k]*Dv*vacdensity_0[k]*(rho_m[k] + rho_i[k]);
    }

    if (K0 == 0) {
      interdensity = 0;
      vacdensity = vacdensity_th;
    }
    else {
      interdensity = (K0 / temp1);
      vacdensity = (K0 + temp3)/temp2;
    }

    // (Vik) Calculate the first estimate of gammadot_climb for each of the slip system
    for (unsigned int k = 0; k < _num_slip_sys ; k++) {
      dummy1 = zi[k] * Di * interdensity;
      dummy2 = zv[k] * Dv * (vacdensity - vacdensity_0[k] );
      vel_climb[k] = (atomvol/b_mag)* (dummy1 - dummy2);

      if(vel_climb[k] < 1e-25){
        vel_climb[k] = 1e-25;
      }
      gamma_dot_c[k] = rho_m[k] * b_mag  *vel_climb[k];

      if (gamma_dot_c[k] <= 1e-18){
        gamma_dot_c[k] = 1.0e-18; }
      }
}

void ThermalIrradiationCPUpdate::calc_mfp(
    unsigned int num_slip_sys,
    std::vector<Real> &rho_m,
    std::vector<Real> &rho_i,
    std::vector<Real> &Nloop,
    std::vector<Real> &dloop,
    std::vector<std::vector<Real>> H,
    std::vector<Real> &d_disl) {

    for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      Real sum1 = 0.0;
      for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        sum1 = sum1 + H[ia][ib]*(rho_m[ib] + rho_i[ib]);
      }
      sum1 = sqrt(sum1);
      if(_irad_defects) {
        Real sum2{1e-18};
        Real sum3{1e-18};
        Real sum4{1e-18};
        for (unsigned int ia = 0; ia < _num_slip_sys ; ia++){
          sum2 = sum2 + Nloop[ia]*dloop[ia];
        }
        sum2 = beta_loop*sqrt(sum2);
        sum4 = beta_loop/mfpSC;
        sum3 = sum1 + sum2 + sum4;
        sum1 = sum3;
      }
      d_disl[ia] = 1.e0/sum1;
    }
}
        

