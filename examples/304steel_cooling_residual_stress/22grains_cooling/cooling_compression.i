[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
 construct_side_list_from_node_list=true
 [./fmg]
   type = FileMeshGenerator
   file = 22grain_mesh.e
 []

 [./new_nodeset1]
   type = ExtraNodesetGenerator
   coord = '0 0 0'
   new_boundary = 101
   input = fmg
 [../]
 [./new_nodeset2]
   type = ExtraNodesetGenerator
   coord = '0 15 0'
   new_boundary = 102
   input = new_nodeset1
 [../]
 [./new_nodeset3]
   type = ExtraNodesetGenerator
   coord = '15 0 0'
   new_boundary = 103
   input = new_nodeset2
 [../]
 [./new_nodeset4]
   type = ExtraNodesetGenerator
   coord = '15 15 0'
   new_boundary = 104
   input = new_nodeset3
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
  [./resid_x]
  [../]
  [./resid_y]
  [../]
  [./resid_z]
  [../]
  [./mechanical_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mechanical_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mechanical_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mechanical_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mechanical_strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mechanical_strain_zx]
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
  [./temperature]
    # order = CONSTANT
    # family = MONOMIAL
  [../]
  [./a110x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a200x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a211x]
    order = CONSTANT
    family = MONOMIAL
  [../]
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
  [./a110y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a200y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a211y]
    order = CONSTANT
    family = MONOMIAL
  [../]
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
  [./sdv47]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sdv48]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sdv49]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Physics/SolidMechanics/QuasiStatic]
  [./alpha]
    strain = FINITE
    incremental = true
    use_finite_deform_jacobian = true
    volumetric_locking_correction = false
    decomposition_method = EigenSolution
    save_in = 'resid_x resid_y resid_z'
    block = 10
    eigenstrain_names = alpha_thermal
  [../]
  [./gamma]
    strain = FINITE
    incremental = true
    use_finite_deform_jacobian = true
    volumetric_locking_correction = false
    decomposition_method = EigenSolution
    save_in = 'resid_x resid_y resid_z'
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    eigenstrain_names = gamma_thermal
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
  [./mechanical_strain_zz]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = mechanical_strain_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./mechanical_strain_yy]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = mechanical_strain_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./mechanical_strain_xx]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = mechanical_strain_xx
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
  [./temperature]
    variable = temperature
    type = FunctionAux
    function = temp_func
    execute_on = timestep_begin
  [../]
  [./a110x]
    type = LatticeStrain
    variable = a110x
    h = 1
    k = 1
    l = 0
    g0 = 1
    g1 = 0
    g2 = 0
    eigenstrain_names = alpha_thermal
    execute_on = timestep_end
    block = 10
  [../]
  [./a110y]
    type = LatticeStrain
    variable = a110y
    h = 1
    k = 1
    l = 0
    g0 = 0
    g1 = 1
    g2 = 0
    eigenstrain_names = alpha_thermal
    execute_on = timestep_end
    block = 10
  [../]
  [./a200x]
    type = LatticeStrain
    variable = a200x
    h = 2
    k = 0
    l = 0
    g0 = 1
    g1 = 0
    g2 = 0
    eigenstrain_names = alpha_thermal
    execute_on = timestep_end
    block = 10
  [../]
  [./a200y]
    type = LatticeStrain
    variable = a200y
    h = 2
    k = 0
    l = 0
    g0 = 0
    g1 = 1
    g2 = 0
    eigenstrain_names = alpha_thermal
    execute_on = timestep_end
    block = 10
  [../]
  [./a211x]
    type = LatticeStrain
    variable = a211x
    h = 2
    k = 1
    l = 1
    g0 = 1
    g1 = 0
    g2 = 0
    eigenstrain_names = alpha_thermal
    execute_on = timestep_end
    block = 10
  [../]
  [./a211y]
    type = LatticeStrain
    variable = a211y
    h = 2
    k = 1
    l = 1
    g0 = 0
    g1 = 1
    g2 = 0
    eigenstrain_names = alpha_thermal
    execute_on = timestep_end
    block = 10
  [../]
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
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
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./sdv47]
    type = StateVariable
    variable = sdv47
    sdv_id = 47
    execute_on = timestep_end
  [../]
  [./sdv48]
    type = StateVariable
    variable = sdv48
    sdv_id = 48
    execute_on = timestep_end
  [../]
  [./sdv49]
    type = StateVariable
    variable = sdv49
    sdv_id = 49
    execute_on = timestep_end
  [../]
[]

[UserObjects]
  [./euler_angle]
    type = EulerAngleReader
    file_name = orientations.in
    execute_on = 'initial'
  [../]
[]

[Functions]
  # [./disp_func]
  #   type = ParsedFunction
  #   vars = 'strain_rate t_start u0 delu'
  #   vals = '-0.0002 816 15 -1.882985e-1' # -1.882985e-1 is the disp_y of top face at the end of cooling
  #   value = 'strain_rate*(t-t_start)*(u0+delu)+delu'
  # [../]
  [./temp_func]
    type = PiecewiseLinear
    x = '0 880 5000'
    # x = '0 1 5000'
    y = '1348 293 293'
  [../]
  [./cte_func_alpha]
    # deF_dt = (0.0003865 + 2*0.000001152*(T_range)-3*0.0000000004379*(T_range).^2)*1e-2
    type = ParsedFunction
    symbol_names = 'a1 a21 a22 a31 a32 tref'
    symbol_values = '0.0003865e-2 2 0.000001152e-2 3 0.0000000004379e-2 0'
    expression = 'a1 + a21*a22*(t-tref) - a31*a32*(t-tref)*(t-tref)'
  [../]
  [./cte_func_gamma]
    # deA_dt = (0.0009472 + 2*0.000001031*(T_range)-3*0.0000000002978*(T_range).^2)*1e-2
    type = ParsedFunction
    symbol_names = 'a1 a21 a22 a31 a32 tref'
    symbol_values = '0.0009472e-2 2 0.000001031e-2 3 0.0000000002978e-2 0'
    expression = 'a1 + a21*a22*(t-tref) - a31*a32*(t-tref)*(t-tref)'
  [../]
  [./velocity_y]
    type = ParsedFunction
    expression = '-0.0002*15'
  [../]
  [./dts]
    type = PiecewiseLinear
    # x = '0 20 816 816.001 850 1066'
    # y = '0.02 0.1 2.4 0.001 0.2 0.2'
    x = '0 20 880 880.001 930 1130'
    y = '0.02 0.5 2.4 0.001 0.2 0.2'
  [../]
[]

[BCs]

  [./bc3]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = 3
    value = 0.0
  [../]
  [./bc4]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = 4
    value = 0.0
  [../]
  [./bc5]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = 5
    value = 0.0
  [../]
  [./bc61]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = 6
    value = 0.0
  [../]
  [./bc62]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = 6
    value = 0.0
  [../]
  [./bc63]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = 6
    value = 0.0
  [../]

  [./y_pull]
    # type = PostprocessorDirichletBC
    type = PresetVelocity
    variable = disp_y
    boundary = 2
    # postprocessor = applied_disp
    function = velocity_y
  [../]

[]

[Controls]
  [./disp_control]
    type = TimePeriod
    disable_objects = 'BCs::y_pull'
    start_time = 0
    end_time = 880 # 816
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeCPElasticityTensor
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./thermal_expansion_ferrite]
    type = ComputeInstantaneousThermalExpansionFunctionEigenstrain
    block = 10
    thermal_expansion_function = cte_func_alpha
    stress_free_temperature = 1348
    temperature = temperature
    eigenstrain_name = alpha_thermal
  [../]
  [./thermal_expansion_austenite]
    type = ComputeInstantaneousThermalExpansionFunctionEigenstrain
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    thermal_expansion_function = cte_func_gamma
    stress_free_temperature = 1348
    temperature = temperature
    eigenstrain_name = gamma_thermal
  [../]
  [./CP_austenite]
    type = DDCPStressUpdate
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    propsFile = fcc_props.in
    slipSysFile = fcc_slip_sys.in
    num_slip_sys = 12
    num_state_vars = 86 # 50 + 3*num_slip_sys
    num_props = 30
    temp = temperature
    tol = 1e-7
    eigenstrain_names = gamma_thermal
    EulerAngFileReader = euler_angle
  [../]
  [./CP_ferrite]
    type = DDCPStressUpdate
    block = 10
    propsFile = bcc_props.in
    slipSysFile = bcc_slip_sys.in
    num_slip_sys = 24
    num_state_vars = 122 # 50 + 3*num_slip_sys
    num_props = 30
    temp = temperature
    tol = 1e-7
    eigenstrain_names = alpha_thermal
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
  end_time = 1130
  # dt = 0.2
  num_steps = 10000

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

  [./avg_strain_xx]
    type = ElementAverageValue
    variable = mechanical_strain_xx
  [../]
  [./avg_strain_yy]
    type = ElementAverageValue
    variable = mechanical_strain_yy
  [../]
  [./avg_strain_zz]
    type = ElementAverageValue
    variable = mechanical_strain_zz
  [../]
  [./avg_vm_stress]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./top_disp]
    type = SideAverageValue
    variable = disp_y
    boundary = 2
  [../]
  # [./applied_disp]
  #   type = FunctionValuePostprocessor
  #   function = disp_func
  # [../]
  [./a_strain_zz]
    type = ElementAverageValue
    block = 10
    variable = mechanical_strain_zz
  [../]
  [./a_strain_yy]
    type = ElementAverageValue
    block = 10
    variable = mechanical_strain_yy
  [../]
  [./a_strain_xx]
    type = ElementAverageValue
    block = 10
    variable = mechanical_strain_xx
  [../]
  [./a_strain_xy]
    type = ElementAverageValue
    block = 10
    variable = mechanical_strain_xy
  [../]
  [./a_strain_yz]
    type = ElementAverageValue
    block = 10
    variable = mechanical_strain_yz
  [../]
  [./a_strain_zx]
    type = ElementAverageValue
    block = 10
    variable = mechanical_strain_zx
  [../]
  [./g_strain_zz]
    type = ElementAverageValue
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    variable = mechanical_strain_zz
  [../]
  [./g_strain_yy]
    type = ElementAverageValue
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    variable = mechanical_strain_yy
  [../]
  [./g_strain_xx]
    type = ElementAverageValue
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    variable = mechanical_strain_xx
  [../]
  [./g_strain_xy]
    type = ElementAverageValue
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    variable = mechanical_strain_xy
  [../]
  [./g_strain_yz]
    type = ElementAverageValue
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    variable = mechanical_strain_yz
  [../]
  [./g_strain_zx]
    type = ElementAverageValue
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
    variable = mechanical_strain_zx
  [../]
  [./a_vm_stress]
    type = ElementAverageValue
    variable = vonmises
    block = 10
  [../]
  [./g_vm_stress]
    type = ElementAverageValue
    variable = vonmises
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./a_stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    block = 10
  [../]
  [./g_stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./a_stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    block = 10
  [../]
  [./g_stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./a_stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    block = 10
  [../]
  [./g_stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./temperature]
    type = ElementAverageValue
    variable = temperature
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./a110x]
    type = ElementNormalizedValue
    variable = a110x
    block = 10
  [../]
  [./a200x]
    type = ElementNormalizedValue
    variable = a200x
    block = 10
  [../]
  [./a211x]
    type = ElementNormalizedValue
    variable = a211x
    block = 10
  [../]
  [./g111x]
    type = ElementNormalizedValue
    variable = g111x
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./g200x]
    type = ElementNormalizedValue
    variable = g200x
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./g220x]
    type = ElementNormalizedValue
    variable = g220x
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./g311x]
    type = ElementNormalizedValue
    variable = g311x
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./a110y]
    type = ElementNormalizedValue
    variable = a110y
    block = 10
  [../]
  [./a200y]
    type = ElementNormalizedValue
    variable = a200y
    block = 10
  [../]
  [./a211y]
    type = ElementNormalizedValue
    variable = a211y
    block = 10
  [../]
  [./g111y]
    type = ElementNormalizedValue
    variable = g111y
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./g200y]
    type = ElementNormalizedValue
    variable = g200y
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./g220y]
    type = ElementNormalizedValue
    variable = g220y
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./g311y]
    type = ElementNormalizedValue
    variable = g311y
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./a_rho_m]
    type = ElementAverageValue
    variable = sdv48
    block = 10
  [../]
  [./a_rho_i]
    type = ElementAverageValue
    variable = sdv49
    block = 10
  [../]
  [./g_rho_m]
    type = ElementAverageValue
    variable = sdv48
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./g_rho_i]
    type = ElementAverageValue
    variable = sdv49
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./a_Epeff]
    type = ElementAverageValue
    variable = sdv47
    block = 10
  [../]
  [./g_Epeff]
    type = ElementAverageValue
    variable = sdv47
    block = '1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22'
  [../]
  [./force_x]
    type = NodalSum
    variable = resid_x
    boundary = 2
  [../]
  [./force_y]
    type = NodalSum
    variable = resid_y
    boundary = 2
  [../]
  [./force_z]
    type = NodalSum
    variable = resid_z
    boundary = 2
  [../]
[]

[Outputs]
  file_base = cooling_22gr
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 10
  [./exodus]
   type = Exodus
   time_step_interval = 100
  [../]
[]
