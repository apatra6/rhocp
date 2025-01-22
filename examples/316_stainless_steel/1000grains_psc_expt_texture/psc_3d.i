[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  construct_side_list_from_node_list = true
  [./fmg]
    type = FileMeshGenerator
    file = 1000grains_27000elems.e
    # file = 64grains_64elems.e
  [../]
  [./bot_corner]
    type = ExtraNodesetGenerator
    new_boundary = bot_corner
    input = fmg
    coord = '0.0 0.0 0.0'
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
 
  [./temperature]
    #order = CONSTANT
    #family = MONOMIAL
  [../]
  [./tolerance_aux]
    #order = CONSTANT
    #family = MONOMIAL
  [../]  
  [./eff_pls]
    order = FIRST
    family = MONOMIAL
  [../]

#stress components
  [./stress_xx]
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
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]  
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]    

#strain components  
  [./total_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./total_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./total_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./total_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./total_strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./total_strain_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]

#effective components   
  [./eeq]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./vonmises]
    order = FIRST
    family = MONOMIAL
  [../]

#dislocation and backstress evolution  
  [./ssd]
    order=FIRST
    family=MONOMIAL
  [../]
  [./eff_bs]
    order = FIRST
    family = MONOMIAL
  [../] 

#directional backstress
  [./back1]
    order = FIRST
    family = MONOMIAL
  [../]
  [./back2]
    order = FIRST
    family = MONOMIAL
  [../]
  [./back3]
    order = FIRST
    family = MONOMIAL
  [../]    
  [./back4]
    order = FIRST
    family = MONOMIAL
  [../]
  [./back5]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./back6]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./back7]
    order = FIRST
    family = MONOMIAL
  [../]
  [./back8]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./back9]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./back10]
    order = FIRST
    family = MONOMIAL
  [../]
  [./back11]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./back12]
    order = FIRST
    family = MONOMIAL
  [../]    
    
#euler
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

#lattice strain: austenite
  [./g111x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g200x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g220x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g311x]
    order = CONSTANT
    family = MONOMIAL
  [../]

#lattice strain: austentite
  [./g111y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g200y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g220y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g311y]
    order = CONSTANT
    family = MONOMIAL
  [../]

#lattice strain: austentite
  [./g111z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g200z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g220z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./g311z]
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
    decomposition_method = EigenSolution
    eigenstrain_names = gamma_thermal
  [../]
[]

[AuxKernels]
  [./temperature]
    variable = temperature
    type = FunctionAux
    function = temp_func
    execute_on = timestep_begin
  [../]
  [./tolerance_aux]
    variable = tolerance_aux
    type = FunctionAux
    function = tol_func
    execute_on = timestep_begin
  [../] 

#1. (111)
  [./g111x]
    type = LatticeStrain
    variable = g111x
    h = 1
    k = 1
    l = 1
    g0 = 1
    g1 = 0
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g111y]
    type = LatticeStrain
    variable = g111y
    h = 1
    k = 1
    l = 1
    g0 = 0
    g1 = 1
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g111z]
    type = LatticeStrain
    variable = g111z
    h = 1
    k = 1
    l = 1
    g0 = 0
    g1 = 0
    g2 = 1
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]

#lattice strain in individual planes 
#1. (200)
  [./g200x]
    type = LatticeStrain
    variable = g200x
    h = 2
    k = 0
    l = 0
    g0 = 1
    g1 = 0
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g200y]
    type = LatticeStrain
    variable = g200y
    h = 2
    k = 0
    l = 0
    g0 = 0
    g1 = 1
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g200z]
    type = LatticeStrain
    variable = g200z
    h = 2
    k = 0
    l = 0
    g0 = 0
    g1 = 0
    g2 = 1
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]

#lattice strain in individual planes 
#1. (220)
  [./g220x]
    type = LatticeStrain
    variable = g220x
    h = 2
    k = 2
    l = 0
    g0 = 1
    g1 = 0
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g220y]
    type = LatticeStrain
    variable = g220y
    h = 2
    k = 2
    l = 0
    g0 = 0
    g1 = 1
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g220z]
    type = LatticeStrain
    variable = g220z
    h = 2
    k = 2
    l = 0
    g0 = 0
    g1 = 0
    g2 = 1
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]

#lattice strain in individual planes 
#1. (311)
  [./g311x]
    type = LatticeStrain
    variable = g311x
    h = 3
    k = 1
    l = 1
    g0 = 1
    g1 = 0
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g311y]
    type = LatticeStrain
    variable = g311y
    h = 3
    k = 1
    l = 1
    g0 = 0
    g1 = 1
    g2 = 0
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]
  [./g311z]
    type = LatticeStrain
    variable = g311z
    h = 3
    k = 1
    l = 1
    g0 = 0
    g1 = 0
    g2 = 1
    eigenstrain_names = gamma_thermal
    execute_on = timestep_end
  [../]  

#euler angles
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

#stress components
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
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    execute_on = timestep_end
    index_i = 1
    index_j = 2
  [../]  
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../] 

#strain components
  [./total_strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_xx
    execute_on = timestep_end
    index_i = 0
    index_j = 0
  [../]
  [./total_strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]  
  [./total_strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./total_strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_xy
    execute_on = timestep_end
    index_i = 0
    index_j = 1
  [../]
  [./total_strain_xz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_xz
    execute_on = timestep_end
    index_i = 0
    index_j = 2
  [../] 
  [./total_strain_yz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_yz
    execute_on = timestep_end
    index_i = 1
    index_j = 2
  [../]    

#effective components
  [./eeq]
    type = StateVariable
    variable = eeq
    sdv_id = 46
    execute_on = timestep_end
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../] 
  [./eff_pls]
    type = StateVariable
    variable = eff_pls
    sdv_id = 47
    execute_on = timestep_end
  [../]

#SSD and backstress components  
  [./ssd]
    type = StateVariable
    variable = ssd
    sdv_id = 48
    execute_on = timestep_end
  [../]
  [./eff_bs]
    type = StateVariable
    variable = eff_bs
    sdv_id = 49
    execute_on = timestep_end
  [../]

  [./back1]
    type = StateVariable
    variable = back1
    sdv_id = 62
    execute_on = timestep_end
  [../]
  [./back2]
    type = StateVariable
    variable = back2
    sdv_id = 63
    execute_on = timestep_end
  [../]
  [./back3]
    type = StateVariable
    variable = back3
    sdv_id = 64
    execute_on = timestep_end
  [../]
  [./back4]
    type = StateVariable
    variable = back4
    sdv_id = 65
    execute_on = timestep_end
  [../]
  [./back5]
    type = StateVariable
    variable = back5
    sdv_id = 66
    execute_on = timestep_end
  [../]
  [./back6]
    type = StateVariable
    variable = back6
    sdv_id = 67
    execute_on = timestep_end
  [../]
  [./back7]
    type = StateVariable
    variable = back7
    sdv_id = 68
    execute_on = timestep_end
  [../]
  [./back8]
    type = StateVariable
    variable = back8
    sdv_id = 69
    execute_on = timestep_end
  [../]              
  [./back9]
    type = StateVariable
    variable = back9
    sdv_id = 70
    execute_on = timestep_end
  [../] 
  [./back10]
    type = StateVariable
    variable = back10
    sdv_id = 71
    execute_on = timestep_end
  [../]
  [./back11]
    type = StateVariable
    variable = back11
    sdv_id = 72
    execute_on = timestep_end
  [../]              
  [./back12]
    type = StateVariable
    variable = back12
    sdv_id = 73
    execute_on = timestep_end
  [../]                  

[]

[UserObjects]
  [./euler_angle]
    type = EulerAngleReader
    file_name = orientations.in
    # file_name = orientations64.in
    execute_on = 'initial'
  [../]
  [./grain_size]
    type = GrainAreaSize
  [../]
[]

[Functions]
  [temp_func]
     type = ParsedFunction
     value = 'if(t>1 , 300.0 , 345.0-45.0*t)'
  []
  [tol_func]
    type = PiecewiseLinear
    x = '0        0.9999     1      1.0001'
    y = '1.0e-6   1.0e-6   1.0e-5   1.0e-6'
  []
  [./dts]
    type = PiecewiseLinear
    x = '0       0.005   0.99    1      1.0000001    1.0001  285.0   287.0'
    y = '1.0e-4   0.05   0.05  1.0e-8   1.0e-5      0.05     0.05   0.5'
    # 285 seconds corresponds to ~0.25 tensile strain
  [../]
  [./pull_velocity]
    type = ParsedFunction
    value = 'if(t<285.0, -1.0e-3 * 1.0, 0.0)'
    # -1.0e-3 is the strain rate, 1.0 mm is the size
  [../]
  [./cte_func_gamma]
    # deA_dt = (0.0009472 + 2*0.000001031*(T_range)-3*0.0000000002978*(T_range).^2)*1e-2
    type = ParsedFunction
    symbol_names = 'a1 a21 a22 a31 a32 tref'
    symbol_values = '0.0009472e-2 2 0.000001031e-2 3 0.0000000002978e-2 0'
    expression = 'a1 + a21*a22*(t-tref) - a31*a32*(t-tref)*(t-tref)'
  [../]  
[]

[BCs]
 [./zn_roller]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = 'zn_face'
    value = 0.0
 [../]
 [./zp_roller]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = 'zp_face'
    value = 0.0
 [../] 
 [./yn_roller]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = 'yn_face'
    value = 0.0
 [../]
 [./yp_roller]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = 'yp_face'
    value = 0.0
  [../]  
 [./xn_roller]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = 'xn_face'
    value = 0.0
  [../]
 [./xp_roller]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = 'xp_face'
    value = 0.0
  [../]     

# preset velocity
  [./x_push_function]
    type = PresetVelocity
    variable = disp_x
    boundary = xp_face
    function = pull_velocity
  [../]
[]


[Controls]
  [./disp_control]
    type = TimePeriod
    disable_objects = 'BCs::x_push_function'
    start_time = 0
    end_time = 1.0
  [../]

  # xp is loading
  [./disp_control_1]
    type = TimePeriod
    disable_objects = 'BCs::xp_roller'
    start_time = 1.0
    end_time = 1587.0
  [../]  

  # z face is free
  [./disp_control_2]
    type = TimePeriod
    disable_objects = 'BCs::zp_roller'
    start_time = 1.0
    end_time = 1587.0
  [../] 
  [./disp_control_3]
    type = TimePeriod
    disable_objects = 'BCs::zn_roller'
    start_time = 1.0
    end_time = 1587.0
  [../]   
[]


[Materials]
  [./CPStressUpdate]
    type = DDCP_SSD_StressUpdate
    propsFile = fcc_props.in
    slipSysFile = fcc_slip_sys.in
    num_slip_sys = 12
    num_state_vars = 73
    num_props = 25
    temp = temperature # K
    tol = tolerance_aux  #1.0e-6
    EulerAngFileReader = euler_angle
    eigenstrain_name = gamma_thermal
  [../]  
  [./thermal_expansion_austenite]
    type = ComputeInstantaneousThermalExpansionFunctionEigenstrain
    thermal_expansion_function = cte_func_gamma
    stress_free_temperature = 345.0
    temperature = temperature
    eigenstrain_name = gamma_thermal
  [../]    
  [./elasticity_tensor]
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
  # line_search = 'none'

  l_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 12
  nl_forced_its = 1
  l_max_its = 10

  start_time = 0.0
  end_time = 1587.0

  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 1e-12
    cutback_factor_at_failure = 0.2
    growth_factor = 1.2
  [../]

  [./Predictor]
    type = SimplePredictor
    scale = 1
  [../]
[]

[Postprocessors]

#auxkernel variables
  [./temperature]
    type=ElementAverageValue
    variable=temperature
  [../]

#strain
  [./strain_zz]
    type = ElementAverageValue
    variable=total_strain_zz
  [../]
  [./strain_yy]
    type = ElementAverageValue
    variable=total_strain_yy
  [../]
  [./strain_xx]
    type = ElementAverageValue
    variable=total_strain_xx
  [../]
         
#stress
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

#lattice strain austenite
  [./g111x]
    type=ElementNormalizedValue
    variable=g111x
  [../] 
  [./g111y]
    type=ElementNormalizedValue
    variable=g111y
  [../]  
  [./g111z]
    type=ElementNormalizedValue
    variable=g111z
  [../]

  [./g200x]
    type=ElementNormalizedValue
    variable=g200x
  [../] 
  [./g200y]
    type=ElementNormalizedValue
    variable=g200y
  [../]  
  [./g200z]
    type=ElementNormalizedValue
    variable=g200z
  [../]

  [./g220x]
    type=ElementNormalizedValue
    variable=g220x
  [../]
  [./g220y]
    type=ElementNormalizedValue
    variable=g220y
  [../]  
  [./g220z]
    type=ElementNormalizedValue
    variable=g220z
  [../]

  [./g311x]
    type=ElementNormalizedValue
    variable=g311x
  [../]
  [./g311y]
    type=ElementNormalizedValue
    variable=g311y
  [../] 
  [./g311z]
    type=ElementNormalizedValue
    variable=g311z
  [../]
  
  [./ssd]
    type = ElementAverageValue
    variable = ssd
  [../]
[]

[Outputs]
  file_base = out_1000_psc
  # file_base = out_64_tension
  csv = true
  print_linear_residuals = true
  perf_graph = true
  interval = 5
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
