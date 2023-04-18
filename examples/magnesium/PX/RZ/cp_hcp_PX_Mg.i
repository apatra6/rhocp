[Mesh]
  displacements = 'disp_x disp_y disp_z'
  construct_side_list_from_node_list = true
  # parallel_type = distributed
  [./fmg]
    type = FileMeshGenerator
    file = 512grains_512elems.e
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
    scaling = 1e-5
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-5
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-5
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

[Kernels]
  [./TensorMechanics]
    use_displaced_mesh = true
    displacements = 'disp_x disp_y disp_z'
    use_finite_deform_jacobian = true
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
    file_name = orientations512.in
    execute_on = 'initial'
  [../]
[]

[Functions]
  [./top_push]
    type = ParsedFunction
    expression = '-0.00008'
  [../]
  [./traction_function]
    type = PiecewiseLinear
    x = '0  2'
    y = '-1  -1'
  [../]
  [./dts]
    type = PiecewiseLinear
    x = '0 0.005'
    y = '0.00005 0.005'
  [../]
[]

[BCs]

##############################################
#
# X : Left and Right (T)
# Y : Bottom and Top (R)
# Z : Back and Front (Z)
#
##############################################
  # Bottom roller
  [./bottom_roller]
    type = DirichletBC
    variable = disp_y
    boundary = yn_face
    value = 0.0
  [../]

  # Constraint along transverse direction
  [./x_fix1]
    type = DirichletBC
    variable = disp_x
    boundary = xn_face
    value = 0.0
  [../]
  [./x_fix2]
    type = DirichletBC
    variable = disp_x
    boundary = xp_face
    value = 0.0
  [../]

  # fixed BCs
  # corner node fixed in all DOFs
  # [./z_bot]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = bot_corner
  #   value = 0.0
  # [../]
  # [./y_bot]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = bot_corner
  #   value = 0.0
  # [../]
  # [./x_bot]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = bot_corner
  #   value = 0.0
  # [../]

  [./y_push_function]
    type = PresetVelocity
    variable = disp_y
    boundary = yp_face
    function = top_push
  [../]

[]

[Materials]
  [./strain]
    type = ComputeFiniteStrain
    volumetric_locking_correction = false
    displacements = 'disp_x disp_y disp_z'
    block = 'ANY_BLOCK_ID 0'
  [../]
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
    tol = 2e-6
    EulerAngFileReader = euler_angle
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

  #explicit solve (uncomment if explicit solve is performed)
  # scheme = 'explicit-euler'
  # solve_type = 'LINEAR'
  # line_search = 'none'

  #implicit solve
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  l_tol = 1e-10
  nl_abs_tol = 5e-6
  nl_rel_tol = 1e-5
  nl_forced_its = 1
  nl_max_its = 10
  l_max_its = 1000

  start_time = 0.0
  end_time = 1000.0

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
    block = 'ANY_BLOCK_ID 0'
    variable = strain_zz
  [../]
  [./strain_yy]
    type = ElementAverageValue
    block = 'ANY_BLOCK_ID 0'
    variable = strain_yy
  [../]
  [./strain_xx]
    type = ElementAverageValue
    block = 'ANY_BLOCK_ID 0'
    variable = strain_xx
  [../]

  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./Ep_eff]
    type = ElementAverageValue
    variable = Ep_eff
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./rho_m]
    type = ElementAverageValue
    variable = rho_m
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./rho_i]
    type = ElementAverageValue
    variable = rho_i
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./phi1]
    type = ElementAverageValue
    variable = phi1
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./BPhi]
    type = ElementAverageValue
    variable = BPhi
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./phi2]
    type = ElementAverageValue
    variable = phi2
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./tf1]
    type = ElementAverageValue
    variable = tf1
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./tf2]
    type = ElementAverageValue
    variable = tf2
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./tf3]
    type = ElementAverageValue
    variable = tf3
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./tf4]
    type = ElementAverageValue
    variable = tf4
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./tf5]
    type = ElementAverageValue
    variable = tf5
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./tf6]
    type = ElementAverageValue
    variable = tf6
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./twinvariant]
    type = ElementAverageValue
    variable = twinvariant
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./g_act_1]
    type = ElementAverageValue
    variable = g_act_1
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./g_act_2]
    type = ElementAverageValue
    variable = g_act_2
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./g_act_3]
    type = ElementAverageValue
    variable = g_act_3
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./g_act_4]
    type = ElementAverageValue
    variable = g_act_4
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./g_act_5]
    type = ElementAverageValue
    variable = g_act_5
    block = 'ANY_BLOCK_ID 0'
  [../]
[]

[Outputs]
  file_base = out_poly_RZ
  csv = true
  print_linear_residuals = true
  perf_graph = true
  interval = 10
   [out]
       type = Checkpoint
       num_files = 5
       interval = 500
   []
  [./exodus]
    type = Exodus
    interval = 500
  [../]
[]
