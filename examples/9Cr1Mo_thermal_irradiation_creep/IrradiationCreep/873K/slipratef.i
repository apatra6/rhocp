[AuxVariables]
	 [gslipr_1]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_2]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_3]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_4]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_5]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_6]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_7]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_8]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_9]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_10]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_11]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_12]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_13]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_14]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_15]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_16]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_17]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_18]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_19]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_20]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_21]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_22]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_23]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [gslipr_24]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_1]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_2]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_3]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_4]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_5]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_6]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_7]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_8]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_9]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_10]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_11]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_12]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_13]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_14]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_15]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_16]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_17]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_18]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_19]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_20]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_21]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_22]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_23]
	 	 order = CONSTANT
	 	 family = MONOMIAL
	 []

	 [cslipr_24]
	 	 order = CONSTANT
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
	 [gslipr_1]
	 	  type = StateVariable 
	 	 variable = gslipr_1
	 	 sdv_id = 171
	 	 execute_on = timestep_end
	 []
	 [gslipr_2]
	 	  type = StateVariable 
	 	 variable = gslipr_2
	 	 sdv_id = 172
	 	 execute_on = timestep_end
	 []
	 [gslipr_3]
	 	  type = StateVariable 
	 	 variable = gslipr_3
	 	 sdv_id = 173
	 	 execute_on = timestep_end
	 []
	 [gslipr_4]
	 	  type = StateVariable 
	 	 variable = gslipr_4
	 	 sdv_id = 174
	 	 execute_on = timestep_end
	 []
	 [gslipr_5]
	 	  type = StateVariable 
	 	 variable = gslipr_5
	 	 sdv_id = 175
	 	 execute_on = timestep_end
	 []
	 [gslipr_6]
	 	  type = StateVariable 
	 	 variable = gslipr_6
	 	 sdv_id = 176
	 	 execute_on = timestep_end
	 []
	 [gslipr_7]
	 	  type = StateVariable 
	 	 variable = gslipr_7
	 	 sdv_id = 177
	 	 execute_on = timestep_end
	 []
	 [gslipr_8]
	 	  type = StateVariable 
	 	 variable = gslipr_8
	 	 sdv_id = 178
	 	 execute_on = timestep_end
	 []
	 [gslipr_9]
	 	  type = StateVariable 
	 	 variable = gslipr_9
	 	 sdv_id = 179
	 	 execute_on = timestep_end
	 []
	 [gslipr_10]
	 	  type = StateVariable 
	 	 variable = gslipr_10
	 	 sdv_id = 180
	 	 execute_on = timestep_end
	 []
	 [gslipr_11]
	 	  type = StateVariable 
	 	 variable = gslipr_11
	 	 sdv_id = 181
	 	 execute_on = timestep_end
	 []
	 [gslipr_12]
	 	  type = StateVariable 
	 	 variable = gslipr_12
	 	 sdv_id = 182
	 	 execute_on = timestep_end
	 []
	 [gslipr_13]
	 	  type = StateVariable 
	 	 variable = gslipr_13
	 	 sdv_id = 183
	 	 execute_on = timestep_end
	 []
	 [gslipr_14]
	 	  type = StateVariable 
	 	 variable = gslipr_14
	 	 sdv_id = 184
	 	 execute_on = timestep_end
	 []
	 [gslipr_15]
	 	  type = StateVariable 
	 	 variable = gslipr_15
	 	 sdv_id = 185
	 	 execute_on = timestep_end
	 []
	 [gslipr_16]
	 	  type = StateVariable 
	 	 variable = gslipr_16
	 	 sdv_id = 186
	 	 execute_on = timestep_end
	 []
	 [gslipr_17]
	 	  type = StateVariable 
	 	 variable = gslipr_17
	 	 sdv_id = 187
	 	 execute_on = timestep_end
	 []
	 [gslipr_18]
	 	  type = StateVariable 
	 	 variable = gslipr_18
	 	 sdv_id = 188
	 	 execute_on = timestep_end
	 []
	 [gslipr_19]
	 	  type = StateVariable 
	 	 variable = gslipr_19
	 	 sdv_id = 189
	 	 execute_on = timestep_end
	 []
	 [gslipr_20]
	 	  type = StateVariable 
	 	 variable = gslipr_20
	 	 sdv_id = 190
	 	 execute_on = timestep_end
	 []
	 [gslipr_21]
	 	  type = StateVariable 
	 	 variable = gslipr_21
	 	 sdv_id = 191
	 	 execute_on = timestep_end
	 []
	 [gslipr_22]
	 	  type = StateVariable 
	 	 variable = gslipr_22
	 	 sdv_id = 192
	 	 execute_on = timestep_end
	 []
	 [gslipr_23]
	 	  type = StateVariable 
	 	 variable = gslipr_23
	 	 sdv_id = 193
	 	 execute_on = timestep_end
	 []
	 [gslipr_24]
	 	  type = StateVariable 
	 	 variable = gslipr_24
	 	 sdv_id = 194
	 	 execute_on = timestep_end
	 []
	 [cslipr_1]
	 	  type = StateVariable 
	 	 variable = cslipr_1
	 	 sdv_id = 195
	 	 execute_on = timestep_end
	 []
	 [cslipr_2]
	 	  type = StateVariable 
	 	 variable = cslipr_2
	 	 sdv_id = 196
	 	 execute_on = timestep_end
	 []
	 [cslipr_3]
	 	  type = StateVariable 
	 	 variable = cslipr_3
	 	 sdv_id = 197
	 	 execute_on = timestep_end
	 []
	 [cslipr_4]
	 	  type = StateVariable 
	 	 variable = cslipr_4
	 	 sdv_id = 198
	 	 execute_on = timestep_end
	 []
	 [cslipr_5]
	 	  type = StateVariable 
	 	 variable = cslipr_5
	 	 sdv_id = 199
	 	 execute_on = timestep_end
	 []
	 [cslipr_6]
	 	  type = StateVariable 
	 	 variable = cslipr_6
	 	 sdv_id = 200
	 	 execute_on = timestep_end
	 []
	 [cslipr_7]
	 	  type = StateVariable 
	 	 variable = cslipr_7
	 	 sdv_id = 201
	 	 execute_on = timestep_end
	 []
	 [cslipr_8]
	 	  type = StateVariable 
	 	 variable = cslipr_8
	 	 sdv_id = 202
	 	 execute_on = timestep_end
	 []
	 [cslipr_9]
	 	  type = StateVariable 
	 	 variable = cslipr_9
	 	 sdv_id = 203
	 	 execute_on = timestep_end
	 []
	 [cslipr_10]
	 	  type = StateVariable 
	 	 variable = cslipr_10
	 	 sdv_id = 204
	 	 execute_on = timestep_end
	 []
	 [cslipr_11]
	 	  type = StateVariable 
	 	 variable = cslipr_11
	 	 sdv_id = 205
	 	 execute_on = timestep_end
	 []
	 [cslipr_12]
	 	  type = StateVariable 
	 	 variable = cslipr_12
	 	 sdv_id = 206
	 	 execute_on = timestep_end
	 []
	 [cslipr_13]
	 	  type = StateVariable 
	 	 variable = cslipr_13
	 	 sdv_id = 207
	 	 execute_on = timestep_end
	 []
	 [cslipr_14]
	 	  type = StateVariable 
	 	 variable = cslipr_14
	 	 sdv_id = 208
	 	 execute_on = timestep_end
	 []
	 [cslipr_15]
	 	  type = StateVariable 
	 	 variable = cslipr_15
	 	 sdv_id = 209
	 	 execute_on = timestep_end
	 []
	 [cslipr_16]
	 	  type = StateVariable 
	 	 variable = cslipr_16
	 	 sdv_id = 210
	 	 execute_on = timestep_end
	 []
	 [cslipr_17]
	 	  type = StateVariable 
	 	 variable = cslipr_17
	 	 sdv_id = 211
	 	 execute_on = timestep_end
	 []
	 [cslipr_18]
	 	  type = StateVariable 
	 	 variable = cslipr_18
	 	 sdv_id = 212
	 	 execute_on = timestep_end
	 []
	 [cslipr_19]
	 	  type = StateVariable 
	 	 variable = cslipr_19
	 	 sdv_id = 213
	 	 execute_on = timestep_end
	 []
	 [cslipr_20]
	 	  type = StateVariable 
	 	 variable = cslipr_20
	 	 sdv_id = 214
	 	 execute_on = timestep_end
	 []
	 [cslipr_21]
	 	  type = StateVariable 
	 	 variable = cslipr_21
	 	 sdv_id = 215
	 	 execute_on = timestep_end
	 []
	 [cslipr_22]
	 	  type = StateVariable 
	 	 variable = cslipr_22
	 	 sdv_id = 216
	 	 execute_on = timestep_end
	 []
	 [cslipr_23]
	 	  type = StateVariable 
	 	 variable = cslipr_23
	 	 sdv_id = 217
	 	 execute_on = timestep_end
	 []
	 [cslipr_24]
	 	  type = StateVariable 
	 	 variable = cslipr_24
	 	 sdv_id = 218
	 	 execute_on = timestep_end
	 []
	 [cstrainfrac1]
		type = StateVariable 
		variable = cstrainfrac1
		sdv_id = 220
		execute_on = timestep_end
   	 []	
	 [cstrainfrac2]
		type = StateVariable 
		variable = cstrainfrac2
		sdv_id = 221
		execute_on = timestep_end
	[]		 
[]


 [Postprocessors]
	 [gslipr_1]
	 	 type = ElementAverageValue
	 	 variable = gslipr_1
	 []
	 [gslipr_2]
	 	 type = ElementAverageValue
	 	 variable = gslipr_2
	 []
	 [gslipr_3]
	 	 type = ElementAverageValue
	 	 variable = gslipr_3
	 []
	 [gslipr_4]
	 	 type = ElementAverageValue
	 	 variable = gslipr_4
	 []
	 [gslipr_5]
	 	 type = ElementAverageValue
	 	 variable = gslipr_5
	 []
	 [gslipr_6]
	 	 type = ElementAverageValue
	 	 variable = gslipr_6
	 []
	 [gslipr_7]
	 	 type = ElementAverageValue
	 	 variable = gslipr_7
	 []
	 [gslipr_8]
	 	 type = ElementAverageValue
	 	 variable = gslipr_8
	 []
	 [gslipr_9]
	 	 type = ElementAverageValue
	 	 variable = gslipr_9
	 []
	 [gslipr_10]
	 	 type = ElementAverageValue
	 	 variable = gslipr_10
	 []
	 [gslipr_11]
	 	 type = ElementAverageValue
	 	 variable = gslipr_11
	 []
	 [gslipr_12]
	 	 type = ElementAverageValue
	 	 variable = gslipr_12
	 []
	 [gslipr_13]
	 	 type = ElementAverageValue
	 	 variable = gslipr_13
	 []
	 [gslipr_14]
	 	 type = ElementAverageValue
	 	 variable = gslipr_14
	 []
	 [gslipr_15]
	 	 type = ElementAverageValue
	 	 variable = gslipr_15
	 []
	 [gslipr_16]
	 	 type = ElementAverageValue
	 	 variable = gslipr_16
	 []
	 [gslipr_17]
	 	 type = ElementAverageValue
	 	 variable = gslipr_17
	 []
	 [gslipr_18]
	 	 type = ElementAverageValue
	 	 variable = gslipr_18
	 []
	 [gslipr_19]
	 	 type = ElementAverageValue
	 	 variable = gslipr_19
	 []
	 [gslipr_20]
	 	 type = ElementAverageValue
	 	 variable = gslipr_20
	 []
	 [gslipr_21]
	 	 type = ElementAverageValue
	 	 variable = gslipr_21
	 []
	 [gslipr_22]
	 	 type = ElementAverageValue
	 	 variable = gslipr_22
	 []
	 [gslipr_23]
	 	 type = ElementAverageValue
	 	 variable = gslipr_23
	 []
	 [gslipr_24]
	 	 type = ElementAverageValue
	 	 variable = gslipr_24
	 []
	 [cslipr_1]
	 	 type = ElementAverageValue
	 	 variable = cslipr_1
	 []
	 [cslipr_2]
	 	 type = ElementAverageValue
	 	 variable = cslipr_2
	 []
	 [cslipr_3]
	 	 type = ElementAverageValue
	 	 variable = cslipr_3
	 []
	 [cslipr_4]
	 	 type = ElementAverageValue
	 	 variable = cslipr_4
	 []
	 [cslipr_5]
	 	 type = ElementAverageValue
	 	 variable = cslipr_5
	 []
	 [cslipr_6]
	 	 type = ElementAverageValue
	 	 variable = cslipr_6
	 []
	 [cslipr_7]
	 	 type = ElementAverageValue
	 	 variable = cslipr_7
	 []
	 [cslipr_8]
	 	 type = ElementAverageValue
	 	 variable = cslipr_8
	 []
	 [cslipr_9]
	 	 type = ElementAverageValue
	 	 variable = cslipr_9
	 []
	 [cslipr_10]
	 	 type = ElementAverageValue
	 	 variable = cslipr_10
	 []
	 [cslipr_11]
	 	 type = ElementAverageValue
	 	 variable = cslipr_11
	 []
	 [cslipr_12]
	 	 type = ElementAverageValue
	 	 variable = cslipr_12
	 []
	 [cslipr_13]
	 	 type = ElementAverageValue
	 	 variable = cslipr_13
	 []
	 [cslipr_14]
	 	 type = ElementAverageValue
	 	 variable = cslipr_14
	 []
	 [cslipr_15]
	 	 type = ElementAverageValue
	 	 variable = cslipr_15
	 []
	 [cslipr_16]
	 	 type = ElementAverageValue
	 	 variable = cslipr_16
	 []
	 [cslipr_17]
	 	 type = ElementAverageValue
	 	 variable = cslipr_17
	 []
	 [cslipr_18]
	 	 type = ElementAverageValue
	 	 variable = cslipr_18
	 []
	 [cslipr_19]
	 	 type = ElementAverageValue
	 	 variable = cslipr_19
	 []
	 [cslipr_20]
	 	 type = ElementAverageValue
	 	 variable = cslipr_20
	 []
	 [cslipr_21]
	 	 type = ElementAverageValue
	 	 variable = cslipr_21
	 []
	 [cslipr_22]
	 	 type = ElementAverageValue
	 	 variable = cslipr_22
	 []
	 [cslipr_23]
	 	 type = ElementAverageValue
	 	 variable = cslipr_23
	 []
	 [cslipr_24]
	 	 type = ElementAverageValue
	 	 variable = cslipr_24
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
