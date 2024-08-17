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
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zx]
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

  [./eff_strain]
     order = CONSTANT
     family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ep_eff]
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
  [./phi1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./BPhi]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./phi2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tf1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tf2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tf3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tf4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tf5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tf6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./twinvariant]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g_act_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g_act_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g_act_3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g_act_4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g_act_5]
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
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xy
    execute_on = timestep_end
    index_i = 0
    index_j = 1
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yz
    execute_on = timestep_end
    index_i = 1
    index_j = 2
  [../]
  [./strain_zx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_zx
    execute_on = timestep_end
    index_i = 0
    index_j = 2
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    execute_on = timestep_end
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./Ep_eff]
    type = StateVariable
    variable = Ep_eff
    sdv_id = 47
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
  [./phi1]
    type = OutputEuler
    variable = phi1
    angle_id = 0
    execute_on = timestep_end
  [../]
  [./BPhi]
    type = OutputEuler
    variable = BPhi
    angle_id = 1
    execute_on = timestep_end
  [../]
  [./phi2]
    type = OutputEuler
    variable = phi2
    angle_id = 2
    execute_on = timestep_end
  [../]
  [./tf1]
    type = StateVariable
    variable = tf1
    sdv_id = 105
    execute_on = timestep_end
  [../]
  [./tf2]
    type = StateVariable
    variable = tf2
    sdv_id = 106
    execute_on = timestep_end
  [../]
  [./tf3]
    type = StateVariable
    variable = tf3
    sdv_id = 107
    execute_on = timestep_end
  [../]
  [./tf4]
    type = StateVariable
    variable = tf4
    sdv_id = 108
    execute_on = timestep_end
  [../]
  [./tf5]
    type = StateVariable
    variable = tf5
    sdv_id = 109
    execute_on = timestep_end
  [../]
  [./tf6]
    type = StateVariable
    variable = tf6
    sdv_id = 110
    execute_on = timestep_end
  [../]
  [./twinvariant]
    type = StateVariable
    variable = twinvariant
    sdv_id = 111
    execute_on = timestep_end
  [../]
  [./g_act_1]
    type = StateVariable
    variable = g_act_1
    sdv_id = 113
    execute_on = timestep_end
  [../]
  [./g_act_2]
    type = StateVariable
    variable = g_act_2
    sdv_id = 114
    execute_on = timestep_end
  [../]
  [./g_act_3]
    type = StateVariable
    variable = g_act_3
    sdv_id = 115
    execute_on = timestep_end
  [../]
  [./g_act_4]
    type = StateVariable
    variable = g_act_4
    sdv_id = 116
    execute_on = timestep_end
  [../]
  [./g_act_5]
    type = StateVariable
    variable = g_act_5
    sdv_id = 117
    execute_on = timestep_end
  [../]

[]


[UserObjects]
  [./euler_angle]
    type = EulerAngleReader
    file_name = orientations_SX.in
    execute_on = 'initial'
  [../]
  [./grain_size] # Note: this does not work with restart, comment if restarting
    type = GrainAreaSize
  [../]
[]

[Functions]
  [./top_push]
    type = ParsedFunction
    expression = '-0.0001*10' # 10 is the sample dimension, 1e-4/s is the strain
  [../]
  [./dts]
    type = PiecewiseLinear
    x = '0 0.001'
    y = '0.00005 0.1'
  [../]
[]

[BCs]

##############################################
# X : Left and Right
# Y : Bottom and Top
# Z : Back and Front
##############################################
# plane strain compression
# compression along y-direction on 'top' boundary
# constraint along x-direction on 'right' and 'left' boundaries
# z-direction is unconstrained
##############################################
  [./y_roller]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./x_roller1]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0.0
  [../]
  [./x_roller2]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./y_push_function]
    type = PresetVelocity
    variable = disp_y
    boundary = top
    function = top_push
  [../]

[]

[Materials]
  [./elasticity_tensor]
    type = ComputeCPElasticityTensor
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./CPMaterial]
    type = DDCPHCPStressUpdate
    block = 'ANY_BLOCK_ID 0'
    num_slip_sys = 24
    num_twin_sys = 6
    propsFile = Mg_props.in
    slipSysFile = Mg_slip_sys.in
    num_state_vars = 117
    num_props = 109
    temp = 300 # K
    tol = 5e-7
    EulerAngFileReader = euler_angle
    GrainAreaSize = grain_size
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

  # # explicit solve (uncomment if explicit solve is performed)
  # scheme = 'explicit-euler'
  # solve_type = 'LINEAR'
  # line_search = 'none'

  # implicit solve
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  l_tol = 1e-10
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-5
  nl_forced_its = 1
  nl_max_its = 10
  l_max_its = 1000

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
  [./strain_xy]
    type = ElementAverageValue
    variable = strain_xy
  [../]
  [./strain_yz]
    type = ElementAverageValue
    variable = strain_yz
  [../]
  [./strain_zx]
    type = ElementAverageValue
    variable = strain_zx
  [../]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./Ep_eff]
    type = ElementAverageValue
    variable = Ep_eff
  [../]
  [./rho_m]
    type = ElementAverageValue
    variable = rho_m
  [../]
  [./rho_i]
    type = ElementAverageValue
    variable = rho_i
  [../]
  [./phi1]
    type = ElementAverageValue
    variable = phi1
  [../]
  [./BPhi]
    type = ElementAverageValue
    variable = BPhi
  [../]
  [./phi2]
    type = ElementAverageValue
    variable = phi2
  [../]
  [./tf1]
    type = ElementAverageValue
    variable = tf1
  [../]
  [./tf2]
    type = ElementAverageValue
    variable = tf2
  [../]
  [./tf3]
    type = ElementAverageValue
    variable = tf3
  [../]
  [./tf4]
    type = ElementAverageValue
    variable = tf4
  [../]
  [./tf5]
    type = ElementAverageValue
    variable = tf5
  [../]
  [./tf6]
    type = ElementAverageValue
    variable = tf6
  [../]
  [./twinvariant]
    type = ElementAverageValue
    variable = twinvariant
  [../]
  [./g_act_1]
    type = ElementAverageValue
    variable = g_act_1
  [../]
  [./g_act_2]
    type = ElementAverageValue
    variable = g_act_2
  [../]
  [./g_act_3]
    type = ElementAverageValue
    variable = g_act_3
  [../]
  [./g_act_4]
    type = ElementAverageValue
    variable = g_act_4
  [../]
  [./g_act_5]
    type = ElementAverageValue
    variable = g_act_5
  [../]
[]

[Outputs]
  file_base = out_10micron
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 10
  # [out]
  #     type = Checkpoint
  #     num_files = 2
  #     time_step_interval = 50
  # []
  [./exodus]
    type = Exodus
    time_step_interval = 10
  [../]
[]
