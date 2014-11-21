import numpy as np

from scipy.interpolate import interp1d

import geode

import jaccard_plot

def interp1d_wrapper(xs,ys):
	"""This method is needed because the normal interp1d method expects the xs to be monotone# increasing
t
	UPDATE: I realize that I was probably drunk when I wrote this... the xs ARE monotone increasing, so this is probably useless.
	xs: 1D array of x values
	ys: 1D array of y values
	"""
	if xs[1] < xs[0]:
		xs=[xs[1],xs[0]]
		ys=[ys[1],ys[0]]
	return interp1d(xs, ys)

def get_expected_IBD_per_pair_at_dist(xs,expected_IBD_per_pair_x,expected_IBD_per_pair_y):
	"""This method takes a list of distances and returns a list of expected IBD"""
	return_array = np.zeros(len(xs))
	temp_xs = np.copy(xs)
	temp_xs[(temp_xs<expected_IBD_per_pair_x[0]) | (temp_xs >= expected_IBD_per_pair_x[-1])] = expected_IBD_per_pair_x[0] #This is also temporary, will be fixed for extreme cases after.
	index_insert = np.searchsorted(expected_IBD_per_pair_x, temp_xs)
	
	index_insert[index_insert == len(expected_IBD_per_pair_x) ] -= 1 #this is a temporary fix, all these values will be dealt with in the extreme cases section of this method.
	#This line creates a list of functions that interpolate between two points and apply them on the appropriate distance to that funtion
	#return_array = np.array([ f(x) for f,x in zip(map(lambda x_s,y_s: interp1d_wrapper(x_s,y_s),zip(expected_IBD_per_pair_x[index_insert-1],expected_IBD_per_pair_x[index_insert]),zip(expected_IBD_per_pair_y[index_insert-1],expected_IBD_per_pair_y[index_insert])),temp_xs)])
	
	f = interp1d(expected_IBD_per_pair_x,expected_IBD_per_pair_y)
	return_array = f(temp_xs) 
	#Here we fix the extreme cases
	return_array[ xs < expected_IBD_per_pair_x[0] ] = expected_IBD_per_pair_y[0]
	return_array[ xs >= expected_IBD_per_pair_x[-1] ] = expected_IBD_per_pair_y[-1]
	
	return return_array
	
def load_and_filter_grid(grid_file_path,bounds):
	"""Loads and filter a list of coordinates. Returns a list of point within bounds in degrees.  
	file_path: the path to the file to load
	bounds: a list of two points in radiants longitude before latitudes 
	"""
	euclidean_grid = np.loadtxt(grid_file_path,delimiter=' ')
	grid=geode.get_spherical_coor(euclidean_grid)
	new_grid=[]
	for coordinates in grid:
		if jaccard_plot.is_in_bounds(geode.toDegrees(coordinates),bounds):
			new_grid.append(coordinates)
	return np.array(new_grid)

def gen_indiv(grid,EXP_INDIV_PER_POINT):
	"""This function takes a grid (list of points (lon,lat) ) and returns a list of points in the same format. One for each individual that is generated """
	nb_indiv_per_point=np.random.poisson(EXP_INDIV_PER_POINT,len(grid))
	indiv_positions = []
	for pt_index in range(len(grid)):
		for i in range(nb_indiv_per_point[pt_index]):
			indiv_positions.append(grid[pt_index])
			
	#Now we return the position for every individuals
	return np.array(indiv_positions) 

def get_IBD_Matrix(best_dist_matrix,expected_IBD_per_pair_x,expected_IBD_per_pair_y,pertube=False):
	""" """
	triu = np.triu_indices(len(best_dist_matrix),1)
	expected_IBD_per_pair = get_expected_IBD_per_pair_at_dist(best_dist_matrix[triu],expected_IBD_per_pair_x,expected_IBD_per_pair_y)
	IBD_per_pair = expected_IBD_per_pair #[np.random.poisson(pair,1)[0] for pair in expected_IBD_per_pair]
	return_matrix = np.zeros([len(best_dist_matrix)]*2)
	return_matrix[triu] = IBD_per_pair
	if pertube:
		pert=np.random.uniform(0.,1.2,size=len(return_matrix[triu]))
		return_matrix[triu] = return_matrix[triu]*pert*np.random.normal(1, 0.12, size=len(return_matrix[triu])) #np.random.normal(0,0.15,len(IBD_per_pair))
		#return_matrix[triu] = return_matrix[triu]*pert #np.random.normal(0,0.15,len(IBD_per_pair))
	if False: #randomize:
		return_matrix[triu] = np.random.permutation(IBD_per_pair)	
	return return_matrix + return_matrix.T
