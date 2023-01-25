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
  [./backstress01]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress02]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress03]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress04]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress05]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress06]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress07]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress08]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress09]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress10]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress12]
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

  [./backstress01]
    type = StateVariable
    variable = backstress01
    sdv_id = 75
    execute_on = timestep_end
  [../]
  [./backstress02]
    type = StateVariable
    variable = backstress02
    sdv_id = 76
    execute_on = timestep_end
  [../]
  [./backstress03]
    type = StateVariable
    variable = backstress03
    sdv_id = 77
    execute_on = timestep_end
  [../]
  [./backstress04]
    type = StateVariable
    variable = backstress04
    sdv_id = 78
    execute_on = timestep_end
  [../]
  [./backstress05]
    type = StateVariable
    variable = backstress05
    sdv_id = 79
    execute_on = timestep_end
  [../]
  [./backstress06]
    type = StateVariable
    variable = backstress06
    sdv_id = 80
    execute_on = timestep_end
  [../]
  [./backstress07]
    type = StateVariable
    variable = backstress07
    sdv_id = 81
    execute_on = timestep_end
  [../]
  [./backstress08]
    type = StateVariable
    variable = backstress08
    sdv_id = 82
    execute_on = timestep_end
  [../]
  [./backstress09]
    type = StateVariable
    variable = backstress09
    sdv_id = 83
    execute_on = timestep_end
  [../]
  [./backstress10]
    type = StateVariable
    variable = backstress10
    sdv_id = 84
    execute_on = timestep_end
  [../]
  [./backstress11]
    type = StateVariable
    variable = backstress11
    sdv_id = 85
    execute_on = timestep_end
  [../]
  [./backstress12]
    type = StateVariable
    variable = backstress12
    sdv_id = 86
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
  [./top_pull]
    type = PiecewiseConstant
   # x = '0  1.2459 3.61748 5.9881 8.35968 10.73078 13.10236 15.47346'
   # y = '-0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 0.004'
     x = '0  31.15 90.4 149.65 208.9 268.15 327.4 386.65 445.9 505.15 564.4 623.65 682.9 742.15 801.4 860.65 919.9 979.15 1038.4 1097.65
          1156.9 1216.15 1275.4 1334.65 1393.9 1453.15 1512.4 1571.65 1630.9 1690.15 1749.4 1808.65 1867.9 1927.15 1986.4 2045.65 2104.9 2164.15 2223.4 2282.65 2341.9 2445
          2594.2 2743.4 2892.6 3041.8 3191 3340.2 3489.4 3638.6 3787.8 3937'
     y = '-0.00016 0.00016 -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016 -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016  
          -0.00016 0.00016 -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016 -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016  -0.00016 0.00016 -0.00016
           0.00016 -0.00016 0.00016 -0.00016 0.00016 -0.00016 0.00016 -0.00016 0.00016 -0.00016'
  [../]
  [./dts]
    type = PiecewiseLinear
  #  x = '0 0.2 1.245 1.25 1.26 3.61 3.62 3.63 5.98 6.00 6.01 8.35 8.36 8.37 10.73 10.74 10.75 13.10 13.11 13.12 15.47 15.48 15.49'
  #  y = '0.00005 0.005 0.005 0.0005 0.005 0.005 0.0005 0.005 0.005 0.0005 0.005 0.005 0.0005 0.005 0.005 0.0005 0.005 0.005 0.0005 0.005 0.005 0.005 0.0005'
     x = '0 0.2'
     y = '0.00005 0.04'

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
  [./z_pull_function]
    type = PresetVelocity
    variable = disp_z
    boundary = zp_face
    function = top_pull
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
    tol = 5e-6
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
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 10
  start_time = 0.0
  end_time = 3937

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
  [./backstress01]
    type = ElementAverageValue
    variable = backstress01
  [../]
  [./backstress02]
    type = ElementAverageValue
    variable = backstress02
  [../]
  [./backstress03]
    type = ElementAverageValue
    variable = backstress03
  [../]
  [./backstress04]
    type = ElementAverageValue
    variable = backstress04
  [../]
  [./backstress05]
    type = ElementAverageValue
    variable = backstress05
  [../]
  [./backstress06]
    type = ElementAverageValue
    variable = backstress06
  [../]
  [./backstress07]
    type = ElementAverageValue
    variable = backstress07
  [../]
  [./backstress08]
    type = ElementAverageValue
    variable = backstress08
  [../]
  [./backstress09]
    type = ElementAverageValue
    variable = backstress09
  [../]
  [./backstress10]
    type = ElementAverageValue
    variable = backstress10
  [../]
  [./backstress11]
    type = ElementAverageValue
    variable = backstress11
  [../]
  [./backstress12]
    type = ElementAverageValue
    variable = backstress12
  [../]
[]

[Outputs]
  file_base = Cu_cyclic
  csv = true
  print_linear_residuals = true
  perf_graph = true
  interval = 10
  [./exodus]
   type = Exodus
   interval = 100
  [../]
  [./cp]
    type = Checkpoint
    interval = 100
    num_files = 2
  [../]
[]
