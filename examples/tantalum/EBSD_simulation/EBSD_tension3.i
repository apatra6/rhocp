[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  construct_side_list_from_node_list = false
  [./emg]
    type = EBSDMeshGenerator
    filename = 'tantalum_input_original_euler.txt'
  []
  [./bottom_center]
    type = ExtraNodesetGenerator
    new_boundary = bottom_center
    coord = '675.0  0.0  15.0'
    input = emg
  [../]
  [./add_side_sets]
    type = SideSetsFromNormalsGenerator
    normals = '1  0  0
               0  1  0
               0  0  1
              -1  0  0
               0 -1  0
               0  0  -1'
    fixed_normal = false
    new_boundary = 'xp_face yp_face zp_face xn_face yn_face zn_face'
    input= bottom_center
  [../]  
  [./bottom_nodes]
    type = BoundingBoxNodeSetGenerator
    input = add_side_sets
    new_boundary = bottom_nodes
    bottom_left = '675.0 0.0  0.0'
    top_right = '675.0   0.0  15.0'
  [../]  
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-6
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-6
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-6
  [../]
[]

[AuxVariables]
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
  [./strain_xy]
    order = FIRST
    family = MONOMIAL
  [../]
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
  [./stress_xy]
    order = FIRST
    family = MONOMIAL
  [../]
  [./vonmises]
    order = FIRST
    family = MONOMIAL
  [../]
  [./grain_id]
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
  [./rho_m]
    order = FIRST
    family = MONOMIAL
  [../]
  [./rho_i]
    order = FIRST
    family = MONOMIAL
  [../]
  [./Ep_eff]
    order = FIRST
    family = MONOMIAL
  [../]
  [./eeq]
    order = FIRST
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
  [./eeq]
    type = StateVariable
    variable = eeq
    sdv_id = 46
    execute_on = timestep_end
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
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xy
    execute_on = timestep_end
    index_i = 0
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
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./grain_id]
    type = EBSDMeshReaderPointDataAux
    variable = grain_id
    ebsd_reader = ebsd_reader
    data_name = 'feature_id'
    execute_on = timestep_end
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
[]

[UserObjects]
  [./ebsd_reader]
    # Read in the EBSD data. Uses the filename given in the mesh block.
    type = EBSDMeshReader
    execute_on = initial
  [../]
[]

[Functions]
  [./pull]
    type = ParsedFunction
    expression = '5280*1.0e-3' # 5280 \mu m is the sample dimension, 1e-3/s is the strain rate
  [../]

  [./dts]
    type = PiecewiseLinear
    x = '0       0.01'
    y = '0.00001 0.002'
  [../]
[]

[BCs]
  [./y_roller]
    type = DirichletBC
    variable = disp_y
    boundary = yn_face
    value = 0.0
  [../]
  [./x_roller]
    type = DirichletBC
    variable = disp_x
    boundary = yn_face
    value = 0.0
  [../]

  [./z_roller]
    type = DirichletBC
    variable = disp_z
    boundary = zn_face
    value = 0.0
  [../]

  [./y_pull_function]
    type = PresetVelocity
    variable = disp_y
    boundary = yp_face
    function = pull
  [../]
[]

[Materials]
  [./CPStressUpdate]
    type = DDCPStressUpdate
    propsFile = bcc_props.in
    slipSysFile = bcc_slip_sys.in
    num_slip_sys = 24
    num_state_vars = 122 # 50 + 3*num_slip_sys
    num_props = 30
    temp = 298 # K
    tol = 1e-6
    EBSDFileReader = ebsd_reader
    isEulerBunge = 1
  [../]
  [./elasticity_tensor]
    type = ComputeCPElasticityTensor
  [../]
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  ## uncomment for implicit solve
  # solve_type = 'NEWTON'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu superlu_dist'

  ## uncomment for explicit solve
  scheme = 'explicit-euler'
  solve_type = 'LINEAR'
  line_search = 'none'

  petsc_options = '-snes_ksp_ew'

  l_tol = 1e-10
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6
  nl_max_its = 10
  nl_forced_its = 1
  l_max_its = 10000

  start_time = 0.0
  end_time = 100.0
  dt = 0.002
  dtmax = 0.002

  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 1e-8
    cutback_factor_at_failure = 0.5
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
  [./stress_xy]
    type = ElementAverageValue
    variable = stress_xy
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
  [./eeq]
    type = ElementAverageValue
    variable = eeq
  [../]
[]

[Outputs]
  file_base = out_EBSD_3
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 10
  [out]
    type = Checkpoint
    num_files = 5
    time_step_interval = 200
  []
  [./exodus]
    type = Exodus
    time_step_interval = 200
  [../]
[]
