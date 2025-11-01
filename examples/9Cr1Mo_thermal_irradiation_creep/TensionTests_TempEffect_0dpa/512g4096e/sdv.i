[AuxVariables]
	 [gstrain_1]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_2]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_3]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_4]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_5]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_6]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_7]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_8]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_9]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_10]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_11]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_12]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_13]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_14]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_15]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_16]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_17]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_18]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_19]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_20]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_21]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_22]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_23]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [gstrain_24]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []


[]


[AuxKernels]
	 [gstrain_1]
	 	  type = StateVariable 
	 	 variable = gstrain_1
	 	 sdv_id = 124
	 	 execute_on = timestep_end
	 []
	 [gstrain_2]
	 	  type = StateVariable 
	 	 variable = gstrain_2
	 	 sdv_id = 125
	 	 execute_on = timestep_end
	 []
	 [gstrain_3]
	 	  type = StateVariable 
	 	 variable = gstrain_3
	 	 sdv_id = 126
	 	 execute_on = timestep_end
	 []
	 [gstrain_4]
	 	  type = StateVariable 
	 	 variable = gstrain_4
	 	 sdv_id = 127
	 	 execute_on = timestep_end
	 []
	 [gstrain_5]
	 	  type = StateVariable 
	 	 variable = gstrain_5
	 	 sdv_id = 128
	 	 execute_on = timestep_end
	 []
	 [gstrain_6]
	 	  type = StateVariable 
	 	 variable = gstrain_6
	 	 sdv_id = 129
	 	 execute_on = timestep_end
	 []
	 [gstrain_7]
	 	  type = StateVariable 
	 	 variable = gstrain_7
	 	 sdv_id = 130
	 	 execute_on = timestep_end
	 []
	 [gstrain_8]
	 	  type = StateVariable 
	 	 variable = gstrain_8
	 	 sdv_id = 131
	 	 execute_on = timestep_end
	 []
	 [gstrain_9]
	 	  type = StateVariable 
	 	 variable = gstrain_9
	 	 sdv_id = 132
	 	 execute_on = timestep_end
	 []
	 [gstrain_10]
	 	  type = StateVariable 
	 	 variable = gstrain_10
	 	 sdv_id = 133
	 	 execute_on = timestep_end
	 []
	 [gstrain_11]
	 	  type = StateVariable 
	 	 variable = gstrain_11
	 	 sdv_id = 134
	 	 execute_on = timestep_end
	 []
	 [gstrain_12]
	 	  type = StateVariable 
	 	 variable = gstrain_12
	 	 sdv_id = 135
	 	 execute_on = timestep_end
	 []
	 [gstrain_13]
	 	  type = StateVariable 
	 	 variable = gstrain_13
	 	 sdv_id = 136
	 	 execute_on = timestep_end
	 []
	 [gstrain_14]
	 	  type = StateVariable 
	 	 variable = gstrain_14
	 	 sdv_id = 137
	 	 execute_on = timestep_end
	 []
	 [gstrain_15]
	 	  type = StateVariable 
	 	 variable = gstrain_15
	 	 sdv_id = 138
	 	 execute_on = timestep_end
	 []
	 [gstrain_16]
	 	  type = StateVariable 
	 	 variable = gstrain_16
	 	 sdv_id = 139
	 	 execute_on = timestep_end
	 []
	 [gstrain_17]
	 	  type = StateVariable 
	 	 variable = gstrain_17
	 	 sdv_id = 140
	 	 execute_on = timestep_end
	 []
	 [gstrain_18]
	 	  type = StateVariable 
	 	 variable = gstrain_18
	 	 sdv_id = 141
	 	 execute_on = timestep_end
	 []
	 [gstrain_19]
	 	  type = StateVariable 
	 	 variable = gstrain_19
	 	 sdv_id = 142
	 	 execute_on = timestep_end
	 []
	 [gstrain_20]
	 	  type = StateVariable 
	 	 variable = gstrain_20
	 	 sdv_id = 143
	 	 execute_on = timestep_end
	 []
	 [gstrain_21]
	 	  type = StateVariable 
	 	 variable = gstrain_21
	 	 sdv_id = 144
	 	 execute_on = timestep_end
	 []
	 [gstrain_22]
	 	  type = StateVariable 
	 	 variable = gstrain_22
	 	 sdv_id = 145
	 	 execute_on = timestep_end
	 []
	 [gstrain_23]
	 	  type = StateVariable 
	 	 variable = gstrain_23
	 	 sdv_id = 146
	 	 execute_on = timestep_end
	 []
	 [gstrain_24]
	 	  type = StateVariable 
	 	 variable = gstrain_24
	 	 sdv_id = 147
	 	 execute_on = timestep_end
	 []
	
[]


 [Postprocessors]
	 [gstrain_1]
	 	 type = ElementAverageValue
	 	 variable = gstrain_1
	 []
	 [gstrain_2]
	 	 type = ElementAverageValue
	 	 variable = gstrain_2
	 []
	 [gstrain_3]
	 	 type = ElementAverageValue
	 	 variable = gstrain_3
	 []
	 [gstrain_4]
	 	 type = ElementAverageValue
	 	 variable = gstrain_4
	 []
	 [gstrain_5]
	 	 type = ElementAverageValue
	 	 variable = gstrain_5
	 []
	 [gstrain_6]
	 	 type = ElementAverageValue
	 	 variable = gstrain_6
	 []
	 [gstrain_7]
	 	 type = ElementAverageValue
	 	 variable = gstrain_7
	 []
	 [gstrain_8]
	 	 type = ElementAverageValue
	 	 variable = gstrain_8
	 []
	 [gstrain_9]
	 	 type = ElementAverageValue
	 	 variable = gstrain_9
	 []
	 [gstrain_10]
	 	 type = ElementAverageValue
	 	 variable = gstrain_10
	 []
	 [gstrain_11]
	 	 type = ElementAverageValue
	 	 variable = gstrain_11
	 []
	 [gstrain_12]
	 	 type = ElementAverageValue
	 	 variable = gstrain_12
	 []
	 [gstrain_13]
	 	 type = ElementAverageValue
	 	 variable = gstrain_13
	 []
	 [gstrain_14]
	 	 type = ElementAverageValue
	 	 variable = gstrain_14
	 []
	 [gstrain_15]
	 	 type = ElementAverageValue
	 	 variable = gstrain_15
	 []
	 [gstrain_16]
	 	 type = ElementAverageValue
	 	 variable = gstrain_16
	 []
	 [gstrain_17]
	 	 type = ElementAverageValue
	 	 variable = gstrain_17
	 []
	 [gstrain_18]
	 	 type = ElementAverageValue
	 	 variable = gstrain_18
	 []
	 [gstrain_19]
	 	 type = ElementAverageValue
	 	 variable = gstrain_19
	 []
	 [gstrain_20]
	 	 type = ElementAverageValue
	 	 variable = gstrain_20
	 []
	 [gstrain_21]
	 	 type = ElementAverageValue
	 	 variable = gstrain_21
	 []
	 [gstrain_22]
	 	 type = ElementAverageValue
	 	 variable = gstrain_22
	 []
	 [gstrain_23]
	 	 type = ElementAverageValue
	 	 variable = gstrain_23
	 []
	 [gstrain_24]
	 	 type = ElementAverageValue
	 	 variable = gstrain_24
	 []

[]
