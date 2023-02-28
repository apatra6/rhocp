# This Python 3.x script uses the SEACAS libraries for extracting spatio-temporal data for all variables from the exodus output file to a csv file
# You need to install SEACAS from https://github.com/sandialabs/seacas and additional Python libraries before running this script

# Replace file paths appropriately

import sys
sys.path.insert(0, '/home/anirban/projects/seacas/lib')
from exodus3 import *
import numpy as np
from itertools import zip_longest
import pandas as pd
import csv
import time

EXODUS_FILE = "/home/anirban/projects/rhocp/examples/tantalum/temperature_effect/compression_512/298K_sr_5000_512grains/out_298K.e"

exodusHandle = exodus(EXODUS_FILE, 'r', array_type='numpy')
timesteps = exodusHandle.num_times()

def grouper(iterable, n, fillvalue=None):
	args = [iter(iterable)] * n
	return zip_longest(fillvalue=fillvalue, *args)

def find_elem_id(arr):
	result = set(arr[0])
	for curr_set in arr[1:]:
		result.intersection_update(curr_set)
	return list(result)[0]

def find_elem_centroid(nodes):
	lmap = lambda func, *iterable: list(map(func, *iterable))
	coord_map = lmap(exodusHandle.get_coord, nodes)

	centroid = [sum(i[0]/8 for i in coord_map), sum(i[1]/8 for i in coord_map), sum(i[2]/8 for i in coord_map)]
	return centroid

def get_stress_strain_values(grain):
	elem_var_names = exodusHandle.get_element_variable_names()

	df2 = {"time":[]}

	# replace with desired time step number
	time = timesteps

	for elem_name in elem_var_names:
		new_val = exodusHandle.get_element_variable_values(grain, elem_name,time)
		if elem_name in df2.keys():
			df2[elem_name] = np.concatenate((df2[elem_name], new_val))
		else:
			df2[elem_name]= new_val
		num_elem = len(new_val)
	df2["time"] = df2["time"] + [time]*num_elem
	return pd.DataFrame.from_dict(df2)

def main():
	elem_block_ids = exodusHandle.get_elem_blk_ids()
	elem_variable_names = exodusHandle.get_element_variable_names

	connectivity = [] # Nodes which constitute a given element
	localNodeToLocalElems = [] # Elem ids of which the node is a part
	localElemToLocalElems = [] # Lists the neighbors of a given element
	collectLocalElemToLocalElems(exodusHandle, connectivity, localNodeToLocalElems, localElemToLocalElems)

	node_ids = exodusHandle.get_node_id_map()
	elem_ids = exodusHandle.get_elem_id_map()
	blk_ids = exodusHandle.get_elem_blk_ids()

	lmap = lambda func, *iterable: list(map(func, *iterable))

	stress_strain_values_map = lmap(get_stress_strain_values, blk_ids)
	print ("Done creating map.")

	for blk_id in blk_ids:
		elem_conn, num_blk_elems, num_elem_nodes = exodusHandle.get_elem_connectivity(blk_id)

		elem_id = []
		elem_x = []
		elem_y = []
		elem_z = []
		phases = []
		for node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8 in grouper(elem_conn, 8):
			node_to_elems = [localNodeToLocalElems[node_1], localNodeToLocalElems[node_2]
			, localNodeToLocalElems[node_3], localNodeToLocalElems[node_4], localNodeToLocalElems[node_5], localNodeToLocalElems[node_6], localNodeToLocalElems[node_7], localNodeToLocalElems[node_8]]
			elem_id.append(find_elem_id(node_to_elems))
			temp = find_elem_centroid(connectivity[elem_id[-1]])
			elem_x.append(temp[0])
			elem_y.append(temp[1])
			elem_z.append(temp[2])
			phases.append(0)
		stress_strain_values_map[blk_id-1]['elem_id'] = elem_id
		stress_strain_values_map[blk_id-1]['elem_x'] = elem_x
		stress_strain_values_map[blk_id-1]['elem_y'] = elem_y
		stress_strain_values_map[blk_id-1]['elem_z'] = elem_z
		stress_strain_values_map[blk_id-1]['phases'] = phases
		stress_strain_values_map[blk_id-1]['blk_id'] = [blk_id]*len(elem_id)

	print (stress_strain_values_map[0].keys())
	for key in stress_strain_values_map[0].keys():
		print (key, len(stress_strain_values_map[0][key]))

	all_data = pd.concat(stress_strain_values_map)
	print(all_data.dtypes)
	all_data.to_csv("/home/anirban/projects/rhocp/examples/tantalum/temperature_effect/compression_512/298K_sr_5000_512grains/exodus_data_output.csv", index=False)


if __name__ == '__main__':
	start_time = time.time()
	main()
	print("--- %s seconds ---" % (time.time() - start_time))






