[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  displacements = 'disp_x disp_y disp_z'
  construct_side_list_from_node_list = true
  # parallel_type = distributed
  [./fmg]
    type = FileMeshGenerator
    file = 810grains_51840elements_z.inp
  [../]
  [./bot_corner]
    type = ExtraNodesetGenerator
    new_boundary = bot_corner
    input = fmg
    coord = '0 0 0'
  [../]
  [./add_side_sets]
    type = SideSetsFromNormalsGenerator
    normals = '1  0  0
               0  1  0
               0  0  1
              -1  0  0
               0 -1  0
               0  0 -1'
    fixed_normal = false
    new_boundary = 'xp_face yp_face zp_face xn_face yn_face zn_face'
    input=bot_corner
  [../]
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
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./phi1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Phi]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./phi2]
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

  [./Fel11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fel12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fel13]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fel22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fel23]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fel33]
    order = CONSTANT
    family = MONOMIAL
  [../]

  #lattice strains
  [./e111x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e111y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e111z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e200x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e200y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e200z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e220x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e220y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e220z]
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
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    execute_on = timestep_end
    index_i = 0
    index_j = 1
  [../]
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    execute_on = timestep_end
    index_i = 0
    index_j = 2
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    execute_on = timestep_end
    index_i = 1
    index_j = 2
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./phi1]
    type = OutputEuler
    variable = phi1
    angle_id = 0
    execute_on = timestep_end
  [../]
  [./Phi]
    type = OutputEuler
    variable = Phi
    angle_id = 1
    execute_on = timestep_end
  [../]
  [./phi2]
    type = OutputEuler
    variable = phi2
    angle_id = 2
    execute_on = timestep_end
  [../]
  [./rho_m]
    type = StateVariable
    variable = rho_m
    sdv_id = 48
    execute_on = timestep_end
  [../]
  [./rho_i]
    type = StateVariable
    variable = rho_i
    sdv_id = 49
    execute_on = timestep_end
  [../]
  [./Ep_eff]
    type = StateVariable
    variable = Ep_eff
    sdv_id = 47
    execute_on = timestep_end
  [../]

  [./Fel11]
    type = StateVariable
    variable = Fel11
    sdv_id = 19
    execute_on = timestep_end
  [../]
  [./Fel12]
    type = StateVariable
    variable = Fel12
    sdv_id = 20
    execute_on = timestep_end
  [../]
  [./Fel13]
    type = StateVariable
    variable = Fel13
    sdv_id = 21
    execute_on = timestep_end
  [../]
  [./Fel22]
    type = StateVariable
    variable = Fel22
    sdv_id = 23
    execute_on = timestep_end
  [../]
  [./Fel23]
    type = StateVariable
    variable = Fel23
    sdv_id = 24
    execute_on = timestep_end
  [../]
  [./Fel33]
    type = StateVariable
    variable = Fel33
    sdv_id = 27
    execute_on = timestep_end
  [../]

  #lattice strain in individual planes
  # (111)
  [./e111x]
    type = LatticeStrain
    variable = e111x
    h = 1
    k = 1
    l = 1
    g0 = 1
    g1 = 0
    g2 = 0
    execute_on = timestep_end
  [../]
  [./e111y]
    type = LatticeStrain
    variable = e111y
    h = 1
    k = 1
    l = 1
    g0 = 0
    g1 = 1
    g2 = 0
    execute_on = timestep_end
  [../]
  [./e111z]
    type = LatticeStrain
    variable = e111z
    h = 1
    k = 1
    l = 1
    g0 = 0
    g1 = 0
    g2 = 1
    execute_on = timestep_end
  [../]

  #lattice strain in individual planes
  # (200)
  [./e200x]
    type = LatticeStrain
    variable = e200x
    h = 2
    k = 0
    l = 0
    g0 = 1
    g1 = 0
    g2 = 0
    execute_on = timestep_end
  [../]
  [./e200y]
    type = LatticeStrain
    variable = e200y
    h = 2
    k = 0
    l = 0
    g0 = 0
    g1 = 1
    g2 = 0
    execute_on = timestep_end
  [../]
  [./e200z]
    type = LatticeStrain
    variable = e200z
    h = 2
    k = 0
    l = 0
    g0 = 0
    g1 = 0
    g2 = 1
    execute_on = timestep_end
  [../]

  #lattice strain in individual planes
  # (220)
  [./e220x]
    type = LatticeStrain
    variable = e220x
    h = 2
    k = 2
    l = 0
    g0 = 1
    g1 = 0
    g2 = 0
    execute_on = timestep_end
  [../]
  [./e220y]
    type = LatticeStrain
    variable = e220y
    h = 2
    k = 2
    l = 0
    g0 = 0
    g1 = 1
    g2 = 0
    execute_on = timestep_end
  [../]
  [./e220z]
    type = LatticeStrain
    variable = e220z
    h = 2
    k = 2
    l = 0
    g0 = 0
    g1 = 0
    g2 = 1
    execute_on = timestep_end
  [../]
[]

[UserObjects]
  [./euler_angle]
    type = EulerAngleReader
    file_name = texture_810_reduced_14.in
    execute_on = 'initial'
  [../]
[]

[Functions]
  [./top_pull]
    type = ParsedFunction
    expression = '0.0008'
  [../]

  [./dts]
    type = PiecewiseLinear
    x = '0 0.5'
    y = '0.0001 0.05'
  [../]
[]

[BCs]

  [./x_bot]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = xn_face
    value = 0.0
  [../]

  [./y_bot]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = yn_face
    value = 0.0
  [../]

  [./z_bot]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = zn_face
    value = 0.0
  [../]

  [./z_pull_function]
    type = PresetVelocity # FunctionDirichletBC
    variable = disp_z
    boundary = zp_face
    function = top_pull
  [../]

[]

[Materials]
  [./CPStressUpdate]
    type = DDCPTSTStressUpdate
    propsFile = fcc_props_v3.in
    slipSysFile = fcc_slip_sys.in
    num_slip_sys = 12
    num_state_vars = 98 # 50 + 4*num_slip_sys
    num_props = 30
    temp = 300 # K
    tol = 5e-7
    EulerAngFileReader = euler_angle
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

  petsc_options = '-snes_ksp_ew'

  # explicit solve (uncomment if explicit solve is performed)
  # scheme = 'explicit-euler'
  # solve_type = 'LINEAR'
  # line_search = 'none'

  #implicit solve
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  l_tol = 1e-10
  nl_abs_tol = 5e-7
  nl_rel_tol = 1e-5
  nl_forced_its = 1
  nl_max_its = 10
  l_max_its = 1000

  start_time = 0.0
  end_time = 450.0

  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 1e-8
    cutback_factor_at_failure = 0.5
    growth_factor = 1.2
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
  [./stress_xz]
    type = ElementAverageValue
    variable = stress_xz
  [../]
  [./stress_xy]
    type = ElementAverageValue
    variable = stress_xy
  [../]
  [./stress_yz]
    type = ElementAverageValue
    variable = stress_yz
  [../]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./phi_1]
    type = ElementAverageValue
    variable = phi1
  [../]
  [./Phi]
    type = ElementAverageValue
    variable = Phi
  [../]
  [./phi_2]
    type = ElementAverageValue
    variable = phi2
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

  [./Fel11]
    type = ElementAverageValue
    variable = Fel11
  [../]
  [./Fel12]
    type = ElementAverageValue
    variable = Fel12
  [../]
  [./Fel13]
    type = ElementAverageValue
    variable = Fel13
  [../]
  [./Fel22]
    type = ElementAverageValue
    variable = Fel22
  [../]
  [./Fel23]
    type = ElementAverageValue
    variable = Fel23
  [../]
  [./Fel33]
    type = ElementAverageValue
    variable = Fel33
  [../]

  #lattice strain
  [./e111x]
    type=ElementNormalizedValue
    variable=e111x
  [../]
  [./e111y]
    type=ElementNormalizedValue
    variable=e111y
  [../]
  [./e111z]
    type=ElementNormalizedValue
    variable=e111z
  [../]
  [./e200x]
    type=ElementNormalizedValue
    variable=e200x
  [../]
  [./e200y]
    type=ElementNormalizedValue
    variable=e200y
  [../]
  [./e200z]
    type=ElementNormalizedValue
    variable=e200z
  [../]
  [./e220x]
    type=ElementNormalizedValue
    variable=e220x
  [../]
  [./e220y]
    type=ElementNormalizedValue
    variable=e220y
  [../]
  [./e220z]
    type=ElementNormalizedValue
    variable=e220z
  [../]
[]

[Outputs]
  file_base = out_51840el_tension
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 10
   [out]
       type = Checkpoint
       num_files = 3
       time_step_interval = 200
   []
  [./exodus]
    type = Exodus
    time_step_interval = 200
  [../]
[]
