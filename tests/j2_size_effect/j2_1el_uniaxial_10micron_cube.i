[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0.0
  xmax = 10.0
  ymin = 0.0
  ymax = 10.0
  zmin = 0.0
  zmax = 10.0
  nx = 1
  ny = 1
  nz = 1
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
[]

[AuxVariables]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_m]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_i]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ep_eff]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain = FINITE
    incremental = true
    use_finite_deform_jacobian = true
    volumetric_locking_correction = false
  [../]
[]

[AuxKernels]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    execute_on = timestep_end
    index_i = 0
    index_j = 0
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    execute_on = timestep_end
    index_i = 0
    index_j = 0
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./rho_m]
    type = StateVariable
    variable = rho_m
    sdv_id = 41
    execute_on = timestep_end
  [../]
  [./rho_i]
    type = StateVariable
    variable = rho_i
    sdv_id = 42
    execute_on = timestep_end
  [../]
  [./Ep_eff]
    type = StateVariable
    variable = Ep_eff
    sdv_id = 29
    execute_on = timestep_end
  [../]
[]

[Functions]
  [./top_pull]
    type = ParsedFunction
    expression = '10.0*0.0005'
  [../]

  [./dts]
    type = PiecewiseLinear
    x = '0 0.2'
    y = '0.0005 0.1'
  [../]
[]

[BCs]

  [./x_bot]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./y_bot]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./z_bot]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./y_pull_function]
    type = PresetVelocity # FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = top_pull
  [../]

[]

[UserObjects]
  [./grain_size] # Note: this does not work with restart, comment if restarting
    type = GrainAreaSize
  [../]
[]

[Materials]
  [./J2StressUpdate]
    type = DDJ2StressUpdate
    propsFile = j2_props.in
    num_state_vars = 51
    num_props = 27
    temp = 300 # K
    tol = 5e-7
    GrainAreaSize = grain_size
  [../]
  [./elasticity_tensor]
    type = ComputeCPElasticityTensor
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full=true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'NEWTON'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  l_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 10
  start_time = 0.0
  end_time = 100.0

  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 1e-8
    cutback_factor_at_failure = 0.1
    growth_factor = 2
  [../]

  [./Predictor]
    type = SimplePredictor
    scale = 1
  [../]
[]

[Postprocessors]
  [./strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  [../]
  [./strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  [../]
  [./strain_xx]
    type = ElementAverageValue
    variable = strain_xx
  [../]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  [../]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./rho_m]
    type = ElementAverageValue
    variable = rho_m
  [../]
  [./rho_i]
    type = ElementAverageValue
    variable = rho_i
  [../]
  [./Ep_eff]
    type = ElementAverageValue
    variable = Ep_eff
  [../]
[]

[Outputs]
  file_base = out_10micron
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 10
  [./exodus]
    type = Exodus
    time_step_interval = 10
  [../]
[]
