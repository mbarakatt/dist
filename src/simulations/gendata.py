import geode

from simulations import *

from helperplot import *

import jaccard_plot

import os

import sys

import numpy as np

import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

args = sys.argv

def usage():
	print "python gendata.py outputfolder -t threshold [-s seed -p indiv_position_file.txt]"

if len(args) == 1:
	print usage()
	exit(1)

def next_arg(s):
	"""returns the arguments that follows the string s."""
	return args[ args.index(s) + 1 ]
try:
	out_folder = args[1] + '/'
except:
	out_folder = ""

try:
	os.mkdir(out_folder)
except:
	pass

EXP_INDIV_PER_POINT = 4 
try:
	THRESHOLD = float(next_arg('-t')) * (10**7)
except:
	print "No Threshold specified."
	print usage()
	exit(1)
try:
	np.random.seed(int(next_arg('-s')))
except:
	pass

try: 
	HAS_INDIV_POSITION = True 
	indiv_position_path = next_arg('-p') 
except: 
	HAS_INDIV_POSITION = False


indiv_positions = np.loadtxt(indiv_position_path,delimiter='\t')

#Positon of the fake train given in degrees in the format [[start_lon,start_lat],[end_lon,end_lat]]
fake_train = geode.line(-90.2,32.3, -86.3,32.3, inputformat = 'degrees')
#fake_train = geode.line(-87.14, 30.7, -86.47, 35.17, inputformat = 'degrees')
print "This is the position of the invented train: ", fake_train.p1, fake_train.p2
#fake_train = geode.line(0.,0.,0.,0., inputformat = 'degrees')


expected_IBD_per_pair_path = ""

expected_IBD_per_pair_x = np.array([0.05, 51.0, 151.0, 251.0, 351.0, 451.0, 551.0, 651.0, 751.0, 851.0, 951.0, 1051.0, 1151.0, 1251.0, 1351.0, 1451.0, 1551.0, 1651.0])
expected_IBD_per_pair_y = 1000000000 * np.array([0.50018956238238499, 0.10480660594102226, 0.056378384421043225, 0.039495200486426803, 0.031637567838294432, 0.030829774551511858, 0.027024850589105798, 0.024467258780431496, 0.023300694698354666, 0.021073417945136664, 0.020524226100485841, 0.019149564282697996, 0.019069512655836798, 0.018569512655836798, 0.018069512655836798, 0.017569512655836798, 0.017069512655836798, 0.016569512655836798])
#expected_IBD_per_pair_y = np.linspace(0.5,0.2,len(expected_IBD_per_pair_x))

#Format will be 2 columns. First; lons. Second: lats
indiv_position_out_path = "indiv_position.txt"

relatedness_matrix_out_path = "relatedness_matrix.txt"

grid_file_path = "../searchspheres/searchspace17.txt"


bounds = jaccard_plot.BOUNDS["SCCS"]

if not HAS_INDIV_POSITION:
	grid = load_and_filter_grid(grid_file_path,bounds)
	indiv_positions=gen_indiv(grid)
	print "Number of point in grid: %s Number of individual: %s" % (len(grid),len(indiv_positions))

print "Number of individuals in grid:", len(indiv_positions)


indiv_positions_euclidean = geode.get_euclidean_coor(*zip(*indiv_positions))

geodesic_dist_matrix = geode.get_dist(*zip(*indiv_positions))

best_dist_matrix, took_train_mask = geode.computeDistances(fake_train,geodesic_dist_matrix,indiv_positions_euclidean,fullOutput=True)


print "Computing distances done, now generating IBD"	

temp_dist = np.copy(best_dist_matrix)
temp_dist[np.tril_indices(len(temp_dist),0)] = 0
IBD_matrix = get_IBD_Matrix(best_dist_matrix,expected_IBD_per_pair_x,expected_IBD_per_pair_y)



draw_point_IBD(geode.toDegrees(indiv_positions),IBD_matrix, out_folder + 'point_IBD.png' )
np.savetxt(out_folder  + indiv_position_out_path, indiv_positions,delimiter='\t')



show_scatterxy(geodesic_dist_matrix,IBD_matrix,out_folder)

print "test3"
np.savetxt(out_folder + relatedness_matrix_out_path, IBD_matrix ,delimiter='\t')

#Test if real train offers better correlation than no trains
import weirdTest
def test_soundness():
	myweirdTest = weirdTest.weirdTest(geodesic_dist_matrix, IBD_matrix,threshold=THRESHOLD,negateyaxis = True)
	no_trains = myweirdTest([0], geodesic_dist_matrix, [0]) 
	real_train = myweirdTest([0], best_dist_matrix, [0])
	f_soundness = open(out_folder + "soundness.txt",'w')
	f_soundness.write("W/O train: " + str(no_trains) + "W real train:" + str(real_train))
	f_soundness.write("spear wo: " + str(myweirdTest.spearman(geodesic_dist_matrix)) + " spear w:" + str(myweirdTest.spearman(best_dist_matrix)) )

test_soundness()


print "Done gendata.py"
