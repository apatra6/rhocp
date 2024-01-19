#include "ComputeStressBase.h"
#include "EulerAngleReader.h"
#include "EBSDMeshReader.h"
#include "GrainAreaSize.h"

class DDCPStressUpdate : public ComputeStressBase
{
public:
  static InputParameters validParams();

  DDCPStressUpdate(const InputParameters & parameters);
  virtual ~DDCPStressUpdate();

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

  void NR_residual (unsigned int num_slip_sys, std::vector<std::vector<Real>> &xs0, std::vector<std::vector<Real>> &xm0, Real temp, Real dt, std::vector<Real> gamma_dot, RankTwoTensor F1, RankTwoTensor &F_el, RankTwoTensor &F_p_inv, RankTwoTensor F_p_inv0, Real C[3][3][3][3], std::vector<Real> &rho_m0, std::vector<Real> &rho_m, std::vector<Real> &rho_i0, std::vector<Real> &rho_i, std::vector<Real> &bstress0, std::vector<Real> &bstress, RankTwoTensor &sig, std::vector<Real> &tau, std::vector<Real> &tau_eff, std::vector<Real> &s_a, std::vector<Real> &s_t, std::vector<std::vector<Real>> A, std::vector<std::vector<Real>> H, std::vector<Real> &residual, Real &sse);

  Real tolerance;

  Real act_vol;
  Real delF0;

  Real C11, C11_perK, C12, C12_perK, C44, C44_perK, G, G_perK, b_mag, gammadot0g, enthalpy_const, p, q, p0, tau0, hp_coeff, grain_size, frictional_stress, k_rho, k_I, Alatent, rho_m_zero, rho_i_zero, d_disl_zero, k_M, R_c, k_ann, k_D, k_bs1, k_bs2, B_k, freq;

  Real sse;
};
