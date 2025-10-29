#include "ComputeStressBase.h"
#include "EulerAngleReader.h"
#include "EBSDMeshReader.h"
#include "GrainAreaSize.h"
//#include "IntMatFileReader.h"

class ThermalIrradiationCPUpdate : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ThermalIrradiationCPUpdate(const InputParameters & parameters);
  virtual ~ThermalIrradiationCPUpdate();

protected:

  FileName _propsFile;
  FileName _slipSysFile;

  unsigned int _num_props;
  unsigned int _num_slip_sys;
  unsigned int _num_state_vars;

  int _grainid;

  const Real _tol;
  const VariableValue & _temp;

  const EulerAngleReader * _EulerAngFileReader;
  const EBSDMeshReader * _EBSDFileReader;
  const GrainAreaSize * _GrainAreaSize;

  Real _grain_size;
  int _isEulerRadian;
  int _isEulerBunge;

  bool _climbmodel;
  bool _transition_slip_on;
  
  
  //const IntMatFileReader * _IntMatFileReader;
  FileName _intMatfile;
  bool _modify_intmat;
  std::vector<Real> _coefficients;

  unsigned int _sliplaw;
  unsigned int _crsslaw;
  bool _irad_defects;
  bool _deltaH_eV;

  virtual void initQpStatefulProperties();
  virtual void computeQpStress();

  Real max_val(Real a,Real b);

  MaterialProperty<Point> & _euler_ang;

  // The eigenstrains
  std::vector<MaterialPropertyName> _eigenstrain_names;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains_old;

  MaterialProperty<std::vector<Real> > & _state_var;
  const MaterialProperty<std::vector<Real> > & _state_var_old;
  MaterialProperty<std::vector<Real> > & _properties;
  const MaterialProperty<std::vector<Real> > & _properties_old;

  // Miller indices of slip plane normals
  MaterialProperty<std::vector<std::vector<Real> > > & _y;
  const MaterialProperty<std::vector<std::vector<Real> > > & _y_old;
  // Miller indices of slip directions
  MaterialProperty<std::vector<std::vector<Real> > > & _z;
  const MaterialProperty<std::vector<std::vector<Real> > > & _z_old;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  MaterialProperty<RankFourTensor> & _Cel_cp;

  // parameters/variables used for calculations
  static const int max_loops = 20;
  void readPropsFile();
  void assignProperties();
  void normalize_vector(Real*, Real*, Real*);
  void rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]);
  void forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]);
  void Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]);
  void aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]);
  void aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]);
  void aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]);
  Real aa_dot_dot_bb(Real a[3][3], Real b[3][3]);
  void bunge_angles(Real (&array1)[3][3], Real (&psi0)[3]);
  Real power(Real x, Real y);
  Real sgn(Real x);

  void NR_residual (
    unsigned int num_slip_sys, 
    std::vector<std::vector<Real>> &xs0, 
    std::vector<std::vector<Real>> &xm0, 
    Real temp, 
    Real dt, 
    std::vector<Real> gamma_dot, 
    std::vector<Real> gammadot_climb,
    RankTwoTensor F1, 
    RankTwoTensor &F_el, 
    RankTwoTensor &F_p_inv, 
    RankTwoTensor F_p_inv0, 
    Real C[3][3][3][3], 
    std::vector<Real> &rho_m0, 
    std::vector<Real> &rho_m, 
    std::vector<Real> &rho_i0, 
    std::vector<Real> &rho_i, 
    std::vector<Real> &bstress0, 
    std::vector<Real> &bstress, 
    std::vector<Real> &Nloop0, 
    std::vector<Real> &Nloop,
    std::vector<Real> &dloop0, 
    std::vector<Real> &dloop,
    RankTwoTensor &sig, 
    std::vector<Real> &tau, 
    std::vector<Real> &climbstress,
    std::vector<Real> &tau_eff, 
    std::vector<Real> &s_a, 
    std::vector<Real> &s_t, 
    std::vector<std::vector<Real>> A, 
    std::vector<std::vector<Real>> H, 
    Real interdensity,
    Real vacdensity,
    Real vacdensity_th,
    Real Di,
    Real Dv,
    std::vector<Real> &residual, 
    Real &sse);

    // Function for Calculation of CRSS
    void calc_crss(
      unsigned int num_slip_sys,
      std::vector<Real> &rho_m,
      std::vector<Real> &rho_i,
      std::vector<Real> &Nloop,
      std::vector<Real> &dloop,
      std::vector<Real> &s_a,
      std::vector<std::vector<Real>> A);

    // Function for Calculation of Derivative of CRSS
void calc_dsadgb(
  unsigned int num_slip_sys,
  std::vector<Real> &rho_m,
  std::vector<Real> &rho_i,
  std::vector<Real> &Nloop,
  std::vector<Real> &dloop,
  std::vector<Real> &s_a,
  std::vector<std::vector<Real>> &A,
  std::vector<std::vector<Real>> &drhomdgb,
  std::vector<std::vector<Real>> &drhoidgb,
  std::vector<std::vector<Real>> &dNldgb,
  std::vector<std::vector<Real>> &ddldgb,
  std::vector<std::vector<Real>> &dsadgb);
    
    // Function for Calculation of Slip Rate
    void calc_sliprate(
      unsigned int num_slip_sys,
      std::vector<Real> &tau_eff,
      std::vector<Real> &s_a,
      std::vector<Real> &s_t,
      std::vector<Real> &gamma_dot_g,
      std::vector<Real> &gamma_dot);

    // Function for Calculation of Climb Rate
    void calc_climbrate(
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
      std::vector<Real> &gammadot_c);


  // Function for Calculation of Mean Free Path
  void calc_mfp(
    unsigned int num_slip_sys,
    std::vector<Real> &rho_m,
    std::vector<Real> &rho_i,
    std::vector<Real> &Nloop,
    std::vector<Real> &dloop,
    std::vector<std::vector<Real>> H,
    std::vector<Real> &d_disl);      

  Real tolerance;

  Real act_vol;
  Real delF0;

  // Constants
  Real
  B_k,              // Boltzman Constant
  freq;             // Debyde Freq

  // Variables Inside Props File
  Real                  // Line no in Props File 
  C11P,                  //1 Elastic Constant
  C11_perK,             //2 Temp Variation of C11
  C12P,                  //3 Elastic Constant
  C12_perK,             //4 Temperature Variation of C12
  C44P,                  //5 Elastic Constant
  C44_perK,             //6 Temperature Variation of C44
  GP,                    //7 Shear Modulus
  G_perK,               //8 Shear Modulus Temp variation
  b_mag,                //9 Burgers Vector Magnitude
  gammadot0g,           //10 reference strain rate
  enthalpy_const,       //11 Multiplication Const for Enthalpy Term
  p,                    //12 Inner Exponent of Slip Law
  q,                    //13 Outer Exponent of Slip Law
  tau0,                 //14 Athermal Slip Resistance Constant
  hp_coeff,             //15 Hall Petch Coefficient
  grain_size,           //16 Grain Size
  frictional_stress,    //17 Lattice Frictional Resistance
  p0,                   //18 Dislocation Barrier Strength
  k_rho,                //19 Taylor Hardening Parameter
  k_I,                  //20 Mean Free Path Constant
  Alatent,              //21 Latent Hardening Coefficient
  rho_m_zero,           //22 Initial Mobile Dislocation Density
  rho_i_zero,           //23 Initial Immobile Dislocation Density
  d_disl_zero,          //24 Initial Dislocation Line Length
  k_M,                  //25 Dislocation Line Generation Constant
  R_c,                  //26 Critical Capture Radius for Dislocation Annhilation
  k_ann,                //27 Dislocation Evlution annhilation constant
  k_D,                  //28 Immobile Dislocation Evolution Dynamic Recovery Constant
  k_bs1,                //29 Backstress Evolution Constant 1
  k_bs2,                //30 Backstress Evolution Constant 2
  Di0,                  //31 Interstitial Diffusion Constant
  Dv0,                  //32 Vacancy Diffusion Constant
  Em_i,                 //33 Migration Enthalpy of Interstitials
  Em_v,                 //34 Migration Enthalpy of Vacancies
  G0,                   //35 Gibbs Free Energy of Vacancy Formation
  zi0,                  //36 Zero Stress Interstitial Capture Efficency
  zv0,                  //37 Zero Stress Vacacny Capture Efficency
  Zsi,                  //38 ClimbStress Propornality Constant of Interstitial Capture Efficency
  K0,                   //39 Point Defect Generation Rates
  Kiv,                  //40 Point Defect Recombination Rates
  atomvol,              //41 Atomic Volume
  lc,                   //42 Average Climb Distance
  Nv,                   //43 No of Atoms Per Unit Volume in the Alloy
  Ecorr,                //44 Correction Factor to Make Energy in J  
  k_D2_0,               //45 Immobile Dislocation Dynamic Recovery Coefficient Due to Climb
  edot0,                //46 Reference Strain Rate for Dynamic Recovery
  n0,                   //47 Dynamic Recovery strain rate exponent
  k_c,                   //48 Mean Free Path of Precipitates
  Nloop_zero,           //49 Initial Density of Dislocation Loops
  dloop_zero,           //50 Initial Diameter of Dislocation Loops
  beta_loop,            //51 Trapping Coefficient of Dislocation Loops
  Ai0_by_i0,            //52 Increase in Area of Dislocation Loop by Addition of a single interstial
  mfpSC,                //53 Hardening by Solute Clusters  
  alpha_loop,           //54 Alpha Loop used for Calculation of CRSS 
  enthalpy_const2,      //55 Multiplication Const for Enthalpy Term for Glide at Low Temperatures
  p2,                   //56 Inner Exponent of Slip Law for Glide at Low Temperatures
  q2,                   //57 Outer Exponent of Slip Law for Glide at Low Temperatures
  transition_temp;      //58 Transition Temperature between Low and High Temperature Slip Laws

  Real k_D2;
  Real C11, C12, C44, G;

  std::vector<std::vector<Real>> sintmat; // Interaction Matrix

  int numrows(std::string);
  int numcols(std::string);
  void readfile(std::string);


  Real sse;
};
