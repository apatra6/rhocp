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

	 [cstrain_1]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_2]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_3]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_4]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_5]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_6]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_7]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_8]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_9]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_10]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_11]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_12]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_13]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_14]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_15]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_16]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_17]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_18]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_19]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_20]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_21]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_22]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_23]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrain_24]
	 	 order = FIRST
	 	 family = MONOMIAL
	 []

	 [cstrainfrac1]
		order = FIRST
		family = MONOMIAL
   	 []

	 [cstrainfrac2]
		order = FIRST
		family = MONOMIAL
	 []	 

[]


[AuxKernels]
	 [gstrain_1]
	 	  type = StateVariable 
	 	 variable = gstrain_1
	 	 sdv_id = 123
	 	 execute_on = timestep_end
	 []
	 [gstrain_2]
	 	  type = StateVariable 
	 	 variable = gstrain_2
	 	 sdv_id = 124
	 	 execute_on = timestep_end
	 []
	 [gstrain_3]
	 	  type = StateVariable 
	 	 variable = gstrain_3
	 	 sdv_id = 125
	 	 execute_on = timestep_end
	 []
	 [gstrain_4]
	 	  type = StateVariable 
	 	 variable = gstrain_4
	 	 sdv_id = 126
	 	 execute_on = timestep_end
	 []
	 [gstrain_5]
	 	  type = StateVariable 
	 	 variable = gstrain_5
	 	 sdv_id = 127
	 	 execute_on = timestep_end
	 []
	 [gstrain_6]
	 	  type = StateVariable 
	 	 variable = gstrain_6
	 	 sdv_id = 128
	 	 execute_on = timestep_end
	 []
	 [gstrain_7]
	 	  type = StateVariable 
	 	 variable = gstrain_7
	 	 sdv_id = 129
	 	 execute_on = timestep_end
	 []
	 [gstrain_8]
	 	  type = StateVariable 
	 	 variable = gstrain_8
	 	 sdv_id = 130
	 	 execute_on = timestep_end
	 []
	 [gstrain_9]
	 	  type = StateVariable 
	 	 variable = gstrain_9
	 	 sdv_id = 131
	 	 execute_on = timestep_end
	 []
	 [gstrain_10]
	 	  type = StateVariable 
	 	 variable = gstrain_10
	 	 sdv_id = 132
	 	 execute_on = timestep_end
	 []
	 [gstrain_11]
	 	  type = StateVariable 
	 	 variable = gstrain_11
	 	 sdv_id = 133
	 	 execute_on = timestep_end
	 []
	 [gstrain_12]
	 	  type = StateVariable 
	 	 variable = gstrain_12
	 	 sdv_id = 134
	 	 execute_on = timestep_end
	 []
	 [gstrain_13]
	 	  type = StateVariable 
	 	 variable = gstrain_13
	 	 sdv_id = 135
	 	 execute_on = timestep_end
	 []
	 [gstrain_14]
	 	  type = StateVariable 
	 	 variable = gstrain_14
	 	 sdv_id = 136
	 	 execute_on = timestep_end
	 []
	 [gstrain_15]
	 	  type = StateVariable 
	 	 variable = gstrain_15
	 	 sdv_id = 137
	 	 execute_on = timestep_end
	 []
	 [gstrain_16]
	 	  type = StateVariable 
	 	 variable = gstrain_16
	 	 sdv_id = 138
	 	 execute_on = timestep_end
	 []
	 [gstrain_17]
	 	  type = StateVariable 
	 	 variable = gstrain_17
	 	 sdv_id = 139
	 	 execute_on = timestep_end
	 []
	 [gstrain_18]
	 	  type = StateVariable 
	 	 variable = gstrain_18
	 	 sdv_id = 140
	 	 execute_on = timestep_end
	 []
	 [gstrain_19]
	 	  type = StateVariable 
	 	 variable = gstrain_19
	 	 sdv_id = 141
	 	 execute_on = timestep_end
	 []
	 [gstrain_20]
	 	  type = StateVariable 
	 	 variable = gstrain_20
	 	 sdv_id = 142
	 	 execute_on = timestep_end
	 []
	 [gstrain_21]
	 	  type = StateVariable 
	 	 variable = gstrain_21
	 	 sdv_id = 143
	 	 execute_on = timestep_end
	 []
	 [gstrain_22]
	 	  type = StateVariable 
	 	 variable = gstrain_22
	 	 sdv_id = 144
	 	 execute_on = timestep_end
	 []
	 [gstrain_23]
	 	  type = StateVariable 
	 	 variable = gstrain_23
	 	 sdv_id = 145
	 	 execute_on = timestep_end
	 []
	 [gstrain_24]
	 	  type = StateVariable 
	 	 variable = gstrain_24
	 	 sdv_id = 146
	 	 execute_on = timestep_end
	 []
	 [cstrain_1]
	 	  type = StateVariable 
	 	 variable = cstrain_1
	 	 sdv_id = 147
	 	 execute_on = timestep_end
	 []
	 [cstrain_2]
	 	  type = StateVariable 
	 	 variable = cstrain_2
	 	 sdv_id = 148
	 	 execute_on = timestep_end
	 []
	 [cstrain_3]
	 	  type = StateVariable 
	 	 variable = cstrain_3
	 	 sdv_id = 149
	 	 execute_on = timestep_end
	 []
	 [cstrain_4]
	 	  type = StateVariable 
	 	 variable = cstrain_4
	 	 sdv_id = 150
	 	 execute_on = timestep_end
	 []
	 [cstrain_5]
	 	  type = StateVariable 
	 	 variable = cstrain_5
	 	 sdv_id = 151
	 	 execute_on = timestep_end
	 []
	 [cstrain_6]
	 	  type = StateVariable 
	 	 variable = cstrain_6
	 	 sdv_id = 152
	 	 execute_on = timestep_end
	 []
	 [cstrain_7]
	 	  type = StateVariable 
	 	 variable = cstrain_7
	 	 sdv_id = 153
	 	 execute_on = timestep_end
	 []
	 [cstrain_8]
	 	  type = StateVariable 
	 	 variable = cstrain_8
	 	 sdv_id = 154
	 	 execute_on = timestep_end
	 []
	 [cstrain_9]
	 	  type = StateVariable 
	 	 variable = cstrain_9
	 	 sdv_id = 155
	 	 execute_on = timestep_end
	 []
	 [cstrain_10]
	 	  type = StateVariable 
	 	 variable = cstrain_10
	 	 sdv_id = 156
	 	 execute_on = timestep_end
	 []
	 [cstrain_11]
	 	  type = StateVariable 
	 	 variable = cstrain_11
	 	 sdv_id = 157
	 	 execute_on = timestep_end
	 []
	 [cstrain_12]
	 	  type = StateVariable 
	 	 variable = cstrain_12
	 	 sdv_id = 158
	 	 execute_on = timestep_end
	 []
	 [cstrain_13]
	 	  type = StateVariable 
	 	 variable = cstrain_13
	 	 sdv_id = 159
	 	 execute_on = timestep_end
	 []
	 [cstrain_14]
	 	  type = StateVariable 
	 	 variable = cstrain_14
	 	 sdv_id = 160
	 	 execute_on = timestep_end
	 []
	 [cstrain_15]
	 	  type = StateVariable 
	 	 variable = cstrain_15
	 	 sdv_id = 161
	 	 execute_on = timestep_end
	 []
	 [cstrain_16]
	 	  type = StateVariable 
	 	 variable = cstrain_16
	 	 sdv_id = 162
	 	 execute_on = timestep_end
	 []
	 [cstrain_17]
	 	  type = StateVariable 
	 	 variable = cstrain_17
	 	 sdv_id = 163
	 	 execute_on = timestep_end
	 []
	 [cstrain_18]
	 	  type = StateVariable 
	 	 variable = cstrain_18
	 	 sdv_id = 164
	 	 execute_on = timestep_end
	 []
	 [cstrain_19]
	 	  type = StateVariable 
	 	 variable = cstrain_19
	 	 sdv_id = 165
	 	 execute_on = timestep_end
	 []
	 [cstrain_20]
	 	  type = StateVariable 
	 	 variable = cstrain_20
	 	 sdv_id = 166
	 	 execute_on = timestep_end
	 []
	 [cstrain_21]
	 	  type = StateVariable 
	 	 variable = cstrain_21
	 	 sdv_id = 167
	 	 execute_on = timestep_end
	 []
	 [cstrain_22]
	 	  type = StateVariable 
	 	 variable = cstrain_22
	 	 sdv_id = 168
	 	 execute_on = timestep_end
	 []
	 [cstrain_23]
	 	  type = StateVariable 
	 	 variable = cstrain_23
	 	 sdv_id = 169
	 	 execute_on = timestep_end
	 []
	 [cstrain_24]
	 	  type = StateVariable 
	 	 variable = cstrain_24
	 	 sdv_id = 170
	 	 execute_on = timestep_end
	 []
	 [cstrainfrac1]
		type = StateVariable 
		variable = cstrainfrac1
		sdv_id = 172
		execute_on = timestep_end
   	 []	
	 [cstrainfrac2]
		type = StateVariable 
		variable = cstrainfrac2
		sdv_id = 173
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
	 [cstrain_1]
	 	 type = ElementAverageValue
	 	 variable = cstrain_1
	 []
	 [cstrain_2]
	 	 type = ElementAverageValue
	 	 variable = cstrain_2
	 []
	 [cstrain_3]
	 	 type = ElementAverageValue
	 	 variable = cstrain_3
	 []
	 [cstrain_4]
	 	 type = ElementAverageValue
	 	 variable = cstrain_4
	 []
	 [cstrain_5]
	 	 type = ElementAverageValue
	 	 variable = cstrain_5
	 []
	 [cstrain_6]
	 	 type = ElementAverageValue
	 	 variable = cstrain_6
	 []
	 [cstrain_7]
	 	 type = ElementAverageValue
	 	 variable = cstrain_7
	 []
	 [cstrain_8]
	 	 type = ElementAverageValue
	 	 variable = cstrain_8
	 []
	 [cstrain_9]
	 	 type = ElementAverageValue
	 	 variable = cstrain_9
	 []
	 [cstrain_10]
	 	 type = ElementAverageValue
	 	 variable = cstrain_10
	 []
	 [cstrain_11]
	 	 type = ElementAverageValue
	 	 variable = cstrain_11
	 []
	 [cstrain_12]
	 	 type = ElementAverageValue
	 	 variable = cstrain_12
	 []
	 [cstrain_13]
	 	 type = ElementAverageValue
	 	 variable = cstrain_13
	 []
	 [cstrain_14]
	 	 type = ElementAverageValue
	 	 variable = cstrain_14
	 []
	 [cstrain_15]
	 	 type = ElementAverageValue
	 	 variable = cstrain_15
	 []
	 [cstrain_16]
	 	 type = ElementAverageValue
	 	 variable = cstrain_16
	 []
	 [cstrain_17]
	 	 type = ElementAverageValue
	 	 variable = cstrain_17
	 []
	 [cstrain_18]
	 	 type = ElementAverageValue
	 	 variable = cstrain_18
	 []
	 [cstrain_19]
	 	 type = ElementAverageValue
	 	 variable = cstrain_19
	 []
	 [cstrain_20]
	 	 type = ElementAverageValue
	 	 variable = cstrain_20
	 []
	 [cstrain_21]
	 	 type = ElementAverageValue
	 	 variable = cstrain_21
	 []
	 [cstrain_22]
	 	 type = ElementAverageValue
	 	 variable = cstrain_22
	 []
	 [cstrain_23]
	 	 type = ElementAverageValue
	 	 variable = cstrain_23
	 []
	 [cstrain_24]
	 	 type = ElementAverageValue
	 	 variable = cstrain_24
	 []
	 [cstrainfrac1]
		type = ElementAverageValue
		variable = cstrainfrac1
     []
	 [cstrainfrac2]
		type = ElementAverageValue
		variable = cstrainfrac2
     []		 	 
[]
