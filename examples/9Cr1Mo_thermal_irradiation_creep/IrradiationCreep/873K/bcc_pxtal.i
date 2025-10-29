# Hoop Stress for Simulation
sigmaV = 13
sigmaH = ${fparse (2/1.732)*sigmaV}


[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = 512grains_4096elems.e
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
  [./phi_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Phi]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./phi_2]
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
  [./euler_angle]
    type = EulerAngleReader
    file_name = BCC_roll_texture2.in
    execute_on = 'initial'
  [../]
[]

[Functions]
  [./top_pull]
    type = ParsedFunction
    expression = '0.4*0.0001' # 0.4 is the sample dimension, 1e-4/s is the strain rate
  [../]

  [./dts]
    type = PiecewiseLinear
    x = '0         0.1  150'
    y = '0.0001    1    623297'
  [../]
  [presspull]
      type = PiecewiseLinear
      x = '0 150 1e6'
      y = '0 1  1'
  []     
[]

[BCs]
  # roller BCs on lateral faces
  [./bottom_roller]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = yn_face
    value = 0.0
  [../]
  [./left_roller]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = xn_face
    value = 0.0
  [../]
  [./back_roller]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = zn_face
    value = 0.0
  [../]

  # fixed BCs
  # corner node fixed in all DOFs
  [./bottom_nodes_x]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = bot_corner
    value = 0.0
  [../]
  [./bottom_nodes_y]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bot_corner
    value = 0.0
  [../]
  [./bottom_nodes_z]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = bot_corner
    value = 0.0
  [../]
  [y_pull_press]
    type = Pressure
    variable = disp_y
    boundary = yp_face
    factor = -${sigmaH}
    function = presspull
  []
  [z_pull_press]
    type = Pressure
    variable = disp_z
    boundary = zp_face
    factor = -${fparse 0.5*sigmaH}
    function = presspull
    []  
  # # tensile loading along y-direction
  # [./y_pull_function]
  #   type = PresetVelocity
  #   variable = disp_y
  #   boundary = yp_face
  #   function = top_pull
  # [../]  

[]

[Materials]
  [./CPStressUpdate]
    type = ThermalIrradiationCPUpdate
    propsFile = bcc_props.in
    slipSysFile = bcc_slip_sys1.in
    num_slip_sys = 24
    num_state_vars = 221 # 50 + 7*num_slip_sys +1 
    num_props = 54
    temp = 873 # K
    climbmodel = true
    tol = 5e-7
    EulerAngFileReader = euler_angle
    irad_defects = true
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

!include'slipratef.i'
!include'irad.i'

[Executioner]
  type = Transient

  solve_type = 'NEWTON'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  # line_search = 'none'

  l_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 10

  start_time = 0.0
  end_time = 1e8

  [./TimeStepper]
    # type = FunctionDT
    # function = dts
    # min_dt = 1e-12
    # cutback_factor_at_failure = 0.2
    # growth_factor = 1.2
    type = LogConstantDT
    log_dt = 0.25
    first_dt = 50
    cutback_factor_at_failure = 0.1
    growth_factor = 2.0    
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
  # [./phi_1]
  #   type = ElementAverageValue
  #   variable = phi_1
  # [../]
  # [./Phi]
  #   type = ElementAverageValue
  #   variable = Phi
  # [../]
  # [./phi_2]
  #   type = ElementAverageValue
  #   variable = phi_2
  # [../]
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
  file_base = out_ICreep_873K_13MPa
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 10
  [./exodus]
    type = Exodus
    time_step_interval = 100
  [../]
  [./cp]
    type = Checkpoint
    time_step_interval = 5
    num_files = 2
  [../]
[]
