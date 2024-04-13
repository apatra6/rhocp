[Mesh]
  displacements = 'disp_x disp_y disp_z'
  construct_side_list_from_node_list = false
  [./emg]
    # Create a mesh representing the EBSD data
    type = EBSDMeshGenerator
    filename = example_1_file_rhocp.txt
  [../]
  [./assignphase]
    # Assign a phase ID based on EBSD data
    type = AssignSubdomainIDfromPhase
    EBSDFilename = example_1_file_rhocp.txt
    input = emg
  [../]  
  [./bottom_nodes]
    type = BoundingBoxNodeSetGenerator
    input = assignphase
    new_boundary = bottom_nodes
    bottom_left = '0 1 0'
    top_right = '0 1 1'
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
    input= bottom_nodes
  [../] 
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    scaling = 1.0e-4
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    scaling = 1.0e-4
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    scaling = 1.0e-4
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
  [./grain_id]
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
  [./grain_id]
    type = EBSDMeshReaderPointDataAux
    variable = grain_id
    ebsd_reader = ebsd_reader
    data_name = 'feature_id'
    execute_on = timestep_end
  [../]  
[]

[UserObjects]
  [./ebsd_reader]
    # Read in the EBSD data. Uses the filename given in the mesh block.
    type = EBSDMeshReader
    execute_on = 'initial'
  [../]
[]

[Functions]
  [./top_push]
    type = ParsedFunction
    expression = '49.5 * 1.0e-3' # 49.5 is the sample dimension, 1.0e-3/s is the strain rate
  [../]

  [./dts]
    type = PiecewiseLinear
    x = '0         1.0'
    y = '0.0001    0.1'
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

  [./y_push_function]
    type = PresetVelocity
    variable = disp_y
    boundary = yp_face
    function = top_push
  [../]
[]

[Materials]
  [./strain]
    # For both phases
    type = ComputeFiniteStrain
    volumetric_locking_correction = false
    displacements = 'disp_x disp_y disp_z'
    block = 'ANY_BLOCK_ID 0' 
  [../]
  [./CPStressUpdate_p1]
    # For phase 1 - p1
    type = DDCPStressUpdate
    propsFile = fcc_props_p1.in
    slipSysFile = fcc_slip_sys_p1.in
    num_slip_sys = 12
    num_state_vars = 86 # 50 + 3*num_slip_sys
    num_props = 30
    temp = 298 # K
    tol = 1e-6
    EBSDFileReader = ebsd_reader
    block = '1'     # For phase 1 - p1
  [../]
  [./CPStressUpdate_p2]
    # For phase 2 - p2
    type = DDCPStressUpdate
    propsFile = fcc_props_p2.in
    slipSysFile = fcc_slip_sys_p2.in
    num_slip_sys = 12
    num_state_vars = 86 # 50 + 3*num_slip_sys
    num_props = 30
    temp = 298 # K
    tol = 1e-6
    EBSDFileReader = ebsd_reader
    block = '2'     # For phase 2 - p2
  [../]  
  [./elasticity_tensor]
    # For both phases
    type = ComputeCPElasticityTensor
    block = 'ANY_BLOCK_ID 0'
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
  [./stress_xx_p1]
    type = ElementAverageValue
    variable = stress_xx
    block = '1'
  [../]
  [./stress_yy_p1]
    type = ElementAverageValue
    variable = stress_yy
    block = '1'
  [../]
  [./stress_zz_p1]
    type = ElementAverageValue
    variable = stress_zz
    block = '1'
  [../]
  [./strain_zz_p1]
    type = ElementAverageValue
    variable = strain_zz
    block = '1'
  [../]
  [./strain_yy_p1]
    type = ElementAverageValue
    variable = strain_yy
    block = '1'
  [../]
  [./strain_xx_p1]
    type = ElementAverageValue
    variable = strain_xx
    block = '1'
  [../]
  [./phi_1_p1]
    type = ElementAverageValue
    variable = phi_1
    block = '1'
  [../]
  [./Phi_p1]
    type = ElementAverageValue
    variable = Phi
    block = '1'
  [../]
  [./phi_2_p1]
    type = ElementAverageValue
    variable = phi_2
    block = '1'
  [../]
  [./rho_m_p1]
    type = ElementAverageValue
    variable = rho_m
    block = '1'
  [../]
  [./rho_i_p1]
    type = ElementAverageValue
    variable = rho_i
    block = '1'
  [../]
[./stress_xx_p2]
    type = ElementAverageValue
    variable = stress_xx
    block = '2'
  [../]
  [./stress_yy_p2]
    type = ElementAverageValue
    variable = stress_yy
    block = '2'
  [../]
  [./stress_zz_p2]
    type = ElementAverageValue
    variable = stress_zz
    block = '2'
  [../]
  [./strain_zz_p2]
    type = ElementAverageValue
    variable = strain_zz
    block = '2'
  [../]
  [./strain_yy_p2]
    type = ElementAverageValue
    variable = strain_yy
    block = '2'
  [../]
  [./strain_xx_p2]
    type = ElementAverageValue
    variable = strain_xx
    block = '2'
  [../]
  [./phi_1_p2]
    type = ElementAverageValue
    variable = phi_1
    block = '2'
  [../]
  [./Phi_p2]
    type = ElementAverageValue
    variable = Phi
    block = '2'
  [../]
  [./phi_2_p2]
    type = ElementAverageValue
    variable = phi_2
    block = '2'
  [../]
  [./rho_m_p2]
    type = ElementAverageValue
    variable = rho_m
    block = '2'
  [../]
  [./rho_i_p2]
    type = ElementAverageValue
    variable = rho_i
    block = '2'
  [../]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]  
  [./Ep_eff]
    type = ElementAverageValue
    variable = Ep_eff
  [../]
[]

[Outputs]
  file_base = ebsd_test
  csv = true
  print_linear_residuals = true
  perf_graph = true
  interval = 1
  [./exodus]
   type = Exodus
   interval = 1
  [../]
  [./cp]
    type = Checkpoint
    interval = 100
    num_files = 2
  [../]
[]
