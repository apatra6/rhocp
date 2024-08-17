[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0.0
  xmax = 0.1
  ymin = 0.0
  ymax = 0.1
  zmin = 0.0
  zmax = 0.1
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
  [./Ep_eff]
    type = StateVariable
    variable = Ep_eff
    sdv_id = 29
    execute_on = timestep_end
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

  [./z_pull_function]
    type = PresetVelocity
    variable = disp_z
    boundary = front
    function = top_pull
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
  file_base = j2_cyclic
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
    time_step_interval = 100
    num_files = 2
  [../]
[]
