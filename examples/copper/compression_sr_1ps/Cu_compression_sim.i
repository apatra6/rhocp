[Mesh]
  displacements = 'disp_x disp_y disp_z'
  construct_side_list_from_node_list = false
  [./fmg]
    type = FileMeshGenerator
    file = 64grains_512elems.e
    # file = 512grains_4096elems.e
  [../]
  [./bot_corner]
    type = ExtraNodesetGenerator
    new_boundary = bot_corner
    input = fmg
    coord = '0 0.0 0.0'
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
  [./stress_xx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = FIRST
    family = MONOMIAL
  [../]
  [./strain_xx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = FIRST
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = FIRST
    family = MONOMIAL
  [../]
  [./vonmises]
    order = FIRST
    family = MONOMIAL
  [../]
  [./Ep_eff]
    order = FIRST
    family = MONOMIAL
  [../]
  [./rho_m]
    order = FIRST
    family = MONOMIAL
  [../]
  [./rho_i]
    order = FIRST
    family = MONOMIAL
  [../]
  [./phi_1]
    order = FIRST
    family = MONOMIAL
  [../]
  [./Phi]
    order = FIRST
    family = MONOMIAL
  [../]
  [./phi_2]
    order = FIRST
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
    strain = FINITE
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
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
  [./phi_1]
    type = OutputEuler
    variable = phi_1
    angle_id = 0
    execute_on = timestep_end
  [../]
  [./Phi]
    type = OutputEuler
    variable = Phi
    angle_id = 1
    execute_on = timestep_end
  [../]
  [./phi_2]
    type = OutputEuler
    variable = phi_2
    angle_id = 2
    execute_on = timestep_end
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
[]

[UserObjects]
  [./euler_angle]
    type = EulerAngleReader
    file_name = orientations64.in # orientations512.in
    execute_on = 'initial'
  [../]
[]

[Functions]
  [./top_push]
    type = ParsedFunction
    value = '-0.4*1' # 0.4 is the sample dimension, 1/s is the strain rate
  [../]

  [./dts]
    type = PiecewiseLinear
    x = '0         0.0001'
    y = '0.0001    0.00001'
  [../]
[]

[BCs]
 [./x_roller]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = yn_face
    value = 0.0
  [../]
  [./y_roller]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = xn_face
    value = 0.0
  [../]
  [./z_roller]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = zn_face
    value = 0.0
  [../]

  [./z_bot]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = bot_corner
    value = 0.0
  [../]

  [./y_bot]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bot_corner
    value = 0.0
  [../]

  [./x_bot]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = bot_corner
    value = 0.0
  [../]
  [./z_push_function]
    type = PresetVelocity
    variable = disp_z
    boundary = zp_face
    function = top_push
  [../]
[]

[Materials]
  [./strain]
    type = ComputeFiniteStrain
    volumetric_locking_correction = false
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./CPStressUpdate]
    type = DDCPStressUpdate
    propsFile = fcc_props.in
    slipSysFile = fcc_slip_sys.in
    num_slip_sys = 12
    num_state_vars = 86 # 50 + 3*num_slip_sys
    num_props = 30
    temp = 298 # K
    tol = 5e-3
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

  solve_type = 'NEWTON'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  l_tol = 1e-8
  nl_abs_tol = 5e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 10
  start_time = 0.0
  end_time = 0.1

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
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]

  [./phi_1]
    type = ElementAverageValue
    variable = phi_1
  [../]
  [./Phi]
    type = ElementAverageValue
    variable = Phi
  [../]
  [./phi_2]
    type = ElementAverageValue
    variable = phi_2
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
[]

[Outputs]
  file_base = Cu_compression_1ps
  csv = true
  print_linear_residuals = true
  perf_graph = true
  interval = 10
  [./exodus]
   type = Exodus
   interval = 50
  [../]
  [./cp]
    type = Checkpoint
    interval = 100
    num_files = 2
  [../]
[]
