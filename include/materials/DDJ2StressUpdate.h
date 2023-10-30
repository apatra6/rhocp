#include "ComputeStressBase.h"
#include "GrainAreaSize.h"

class DDJ2StressUpdate : public ComputeStressBase
{
public:
  static InputParameters validParams();

  DDJ2StressUpdate(const InputParameters & parameters);
  virtual ~DDJ2StressUpdate();

protected:
  FileName _propsFile;
  const GrainAreaSize * _GrainAreaSize;
  unsigned int _num_props;
  unsigned int _num_state_vars;
  const Real _tol;
  const VariableValue & _temp;

  int _grainid;

  Real _grain_size;

  virtual void initQpStatefulProperties();
  virtual void computeQpStress();

  Real max_val(Real a,Real b);

  MaterialProperty<std::vector<Real> > & _state_var;
  const MaterialProperty<std::vector<Real> > & _state_var_old;
  MaterialProperty<std::vector<Real> > & _properties;
  const MaterialProperty<std::vector<Real> > & _properties_old;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  MaterialProperty<RankFourTensor> & _Cel_cp;

  // parameters/variables used for calculations
  void readPropsFile();
  void assignProperties();
  void normalize_vector(Real*, Real*, Real*);
  void rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]);
  void forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]);
  void Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]);
  void aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]);
  void aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]);
  void aa_dot_dot_bbbb(Real a[3][3], Real b[3][3][3][3], Real (&product)[3][3]);
  void aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]);
  Real aa_dot_dot_bb(Real a[3][3], Real b[3][3]);
  Real power(Real x, Real y);
  Real sgn(Real x);

  void NR_residual_J2 (Real temp, Real dt, RankTwoTensor F1, RankTwoTensor &F_p_inv, RankTwoTensor F_p_inv_0, Real C[3][3][3][3], RankTwoTensor Np_star, Real &rho_m0, Real &rho_m, Real &rho_i0, Real &rho_i, RankTwoTensor &bstress0, RankTwoTensor &bstress, RankTwoTensor &sig, Real &sig_vm, Real &sig_vm_star, Real &s_a, Real &s_t, Real &gamma_dot_trial, Real &residual);

  Real tolerance;

  Real act_vol;
  Real delF0;

  Real YM, YM_perK, poisson_ratio, poisson_ratio_perK, G, G_perK, b_mag, gammadot0g, enthalpy_const, p, q, p0, tau0, hp_coeff, grain_size, frictional_stress, k_rho, k_I, Alatent, rho_m_zero, rho_i_zero, d_disl_zero, k_M, R_c, k_ann, k_D, k_bs1, k_bs2, B_k, freq;
};
