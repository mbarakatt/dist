#!/usr/bin/env python2
print("Starting train_geo.py")
print "Importing modules...",
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.optimize
import scipy.special
#from mpl_toolkits.mplot3d import axes3d
from  pylab import *
import random,math,scipy,sys,os
import numpy as np
np.set_printoptions(threshold=np.nan)
np.seterr(divide='ignore')
from itertools import product
from mpl_toolkits.basemap import Basemap, shiftgrid
from jaccard_plot import *
import time
from subprocess import call
#---IMPORTING TESTS-------
#from spearmanFlat import spearmanFlat
#from mantelTest import mantelTest
#from spearmanMantelTest import spearmanMantelTest
from weirdTest import weirdTest
from geode import *
#---END IMPORTING TESTS---
print(" Done")

#importing bounds
#bounds_pt1, bounds_pt2= [(-100, 25), (-70, 42)] #These are the bounds I used for the SCCS dataset
bounds_pt1, bounds_pt2= [(-180, -90), (180, 90)]
#Input should be in degrees
def is_in_bounds(lon,lat): 
	result = bounds_pt1[0] < lon and lon < bounds_pt2[0] and bounds_pt1[1] < lat and lat < bounds_pt2[1]
	#if result:
	#	print lon,lat 
	return bounds_pt1[0] < lon and lon < bounds_pt2[0] and bounds_pt1[1] < lat and lat < bounds_pt2[1]

try:
	mintdist=float(sys.argv[1])
	maxtdist=float(sys.argv[2])
except:
	print "No arguments given, using 0 and 3000"
	mindist = 0
	maxdist = 1

try:
	outfolder = sys.argv[3] + '/'
except:
	outfolder = ""

print "MinDist:",mintdist,"MaxDist:",maxtdist
#---CONSTANT DEFINITION-------
DISTANCE_MASK=20000 #Km
THRESHOLD = float(sys.argv[sys.argv.index("-t") + 1])
print "THRESHOLD:",THRESHOLD
JACCARD_MASK = 1.0 #below this number is included
OUTPUT_FOLDER=outfolder #os.path.join("../results" , time.strftime("%Y.%m.%d.%Hh%mm"))
EARTH_RADIUS=6371
NEGLIGEABLE_RADIAN=0.005
#jaccard_file_path="../data/language/Ruhlen2014jaccard.txt"
#longlatfile = "../data/language/longitudelatitude_americasnegative.txt"
#jaccard_file_path="../SCCS/cM_Matrix_18.txt"
#longlatfile = "../SCCS/SCCS_lonlat.txt"
jaccard_file_path = outfolder + "relatedness_matrix.txt"
pos = np.loadtxt( outfolder + "/indiv_position.txt", delimiter='\t')
#---END CONSTANT DEFINITION---

#---PARSE RELATEDNESS-------
call(['mkdir','-p',OUTPUT_FOLDER])
parse_jaccard(open(jaccard_file_path,'r'))
tempDistJaccard = np.array(parse_jaccard(open(jaccard_file_path,'r')))
maskJaccard = (tempDistJaccard<=JACCARD_MASK) 
distJaccard = tempDistJaccard
applyMaskJaccard = lambda distM : distM
#print "distJaccard", distJaccard
#distJaccard=np.array([[0.,0.5,0.2],[0.5,0.,0.1],[0.2,0.1,0.]])
#---END PARSE RELATEDNESS---

#---PARSING LONG LAT FILE-------
#f_longlatfile = open(longlatfile, 'r')
#lines_longlatfile = f_longlatfile.read().split('\n')
#if lines_longlatfile[-1]=='':
#	lines_longlatfile=lines_longlatfile[:-1]
#pos = [np.array([0.,10.,10.]),np.array([0.,10.,20.])]
#pos =  map(np.array,np.array([map(float,line.split('\t')) for line in lines_longlatfile ]).T)
#print "POS", pos
lons = pos[:,0]#(360+pos[0])*(pos[0]<=0) + pos[0]*(pos[0]>0)
lats = pos[:,1]#(90-pos[1])*(pos[1]>0) + ( -pos[1] + 90 ) * ( pos[1] <= 0 )
#print "lons:",np.round(lons,1), "lats:", np.round(lats,1)

#Include this line if the lons and lats or saved in degree format
#lons, lats = map(lambda x: x*2*np.pi/360,[lons,lats]) #correct that
#---END PARSING LONG LAT FILE---

#distIndiv = parse_file_geo_dist(file("geographicGeoDistLanguages.txt"))
#distIndiv = parse_file_geo_dist(file("../SCCS/geodesicDist.txt"))
distIndiv = get_dist(lons,lats) 

#Destription get_spherical_coor: Return points in lons, lats

#list of points
searchspace_euclidean=np.array(map(lambda x : x.split(" "), open("searchspheres/searchspace3.txt",'r').read().split("\n"))[:-1])#*EARTH_RADIUS
searchspace_euclidean=[map(lambda x : float(x)*EARTH_RADIUS, item) for item in searchspace_euclidean]
#print len(searchspace_euclidean)
#lons, lats
searchspace_radians=np.array(get_spherical_coor(searchspace_euclidean))
searchspace_degree=searchspace_radians * 360.0/(2*np.pi)
searchspace_valid=[]
#print "what" ,searchspace_degree
for i in range(len(searchspace_degree)):
	if is_in_bounds(*searchspace_degree[i]):
		searchspace_valid.append(True)
	else:
		searchspace_valid.append(False)

print 'nbptin search space', np.sum(searchspace_valid)

pos_euclidean = np.array(get_euclidean_coor(lons,lats))
#xs,ys,zs = pos_euclidean

NB_BINS=10

def dist(c,p2):
    return math.hypot(p2[0]-p1[0],p2[1]-p1[1])

def get_dist_2pt_single(p1,p2):
	retval= get_dist_2pt(p1,np.array([p2]))[0]
	return retval

maxdist = DISTANCE_MASK #np.max(b)
bins = np.linspace(0.000001, maxdist, NB_BINS+1) #create NB_BINS bins
dc = (bins + (bins[1]-bins[0])/2.0)[0:-1]
best_min = "nan"
best_sigmoid_params = []


size_dist_L=np.size(distJaccard)
	

def find_index_in_array(array, c):
	for i in range(np.size(array)):
		if c <= array[i]:
			return i
	return len(array)

numberPlotCounter=0
def simplePlot(ax,xs,ys,Title,xLabel,yLabel,savepath,Label,save=True):
	ax.plot(xs, ys,label=Label)
	ax.set_xlabel(xLabel)
	ax.set_ylabel(yLabel)
	ax.set_title(Title)
	ax.set_xlim(0,DISTANCE_MASK)
	if save:
		handles, labels = ax.get_legend_handles_labels()
		ax.legend(handles, labels,loc=4)
		plt.savefig(''.join([c for c in savepath if c not in ["'","[","]"]]) + '.jpg')
		plt.close()

def geoDistVsJaccardDist(lon1, lat1, lon2, lat2, bestDistances,score,tookTrainMask):
	fig = plt.figure()
	ax=fig.add_subplot(111)
	applyMaskTookTrain= lambda M : M #*tookTrainMask
	distIndivFlatten=(applyMaskTookTrain(applyMaskJaccard(distIndiv))).flatten()
	countPerBin = np.histogram(distIndivFlatten, bins)[0]
	distJaccardFlatten=(applyMaskTookTrain(distJaccard)).flatten()
	countJaccardPerBin = np.histogram(distIndivFlatten, bins, weights=distJaccardFlatten)[0]
	simplePlot(ax,dc,countJaccardPerBin/countPerBin, str(map(lambda x : str(x*360/(np.pi*2)),[lon1,lat1,lon2,lat2])),"Geodesic Distance(Km)", "JaccardDistance",os.path.join(OUTPUT_FOLDER,str(map(lambda x : str(x*360/(np.pi*2)),[lon1,lat1,lon2,lat2]))),"Without Train",save=False)
	ax.plot(distIndivFlatten,distJaccardFlatten,'.')
	#print "asdf'asdf" , countJaccardPerBin/countPerBin, countPerBin

	bestDistancesFlatten=(applyMaskTookTrain(applyMaskJaccard(bestDistances))).flatten()
	countPerBin = np.histogram(bestDistancesFlatten, bins)[0]
	#countPerBin[0]-=len(bestDistances)
	countJaccardPerBin = np.histogram(bestDistancesFlatten, bins, weights=distJaccardFlatten)[0]
	ax.plot(bestDistancesFlatten*(bestDistancesFlatten != distIndivFlatten),distJaccardFlatten*(bestDistancesFlatten != distIndivFlatten),'.')
	simplePlot(ax,dc, countJaccardPerBin/countPerBin, str(map(lambda x : str(x*360/(np.pi*2)),[lon1,lat1,lon2,lat2])),"Geodesic Distance(Km)", "JaccardDistance",os.path.join(OUTPUT_FOLDER, str(score)+"," + str(map(lambda x : str(x*360/(np.pi*2)) ,[lon1,lat1,lon2,lat2]))),"With Train")
#	numberPlotCounter+=1
	
def write_train(fout,lon1,lat1,lon2,lat2,score):
	fout.write("%f,%f,%f,%f,%f\n" % (lon1,lat1,lon2,lat2,score) )
	fout.flush()


testList=[weirdTest]
testClassList=[thistest(distIndiv,distJaccard,threshold = THRESHOLD, negateyaxis=True) for thistest in testList]


def computeScores(fout,lon1,lat1,lon2,lat2):
	t0 =time.time()
	lon1, lat1, lon2, lat2 = [x*2*np.pi/360.0 for x in [lon1,lat1,lon2,lat2]]
	length = get_dist_2pt_single([lon1,lat1] , [lon2,lat2])
	myLine = line(lon1,lat1,lon2,lat2)
	bestDistances, tookTrainMask = computeDistances(myLine, distIndiv, pos_euclidean, fullOutput = True)
	testResults = []
	for testClass in testClassList:
		result=testClass(myLine,bestDistances,tookTrainMask)
		#result = testClass.spearman(bestDistances)
		testResults.append(result)

	#print "Time for a whole iteration: ", time.time() - t0, "Results", testResults
	write_train(fout, lon1, lat1, lon2, lat2, testResults[0])
	return testResults

np.seterr(divide = 'ignore', invalid = 'ignore')

def get_vect(lon,lat,lons,lats):
	x = np.array([(lons - lon), ( lats - lat )]).T
	norm_x = get_norm(x)
	result = x.T/norm_x
	n = norm_x==0
	result[0][n]=0
	result[1][n]=0
	return result 


#fout=open("temp/out"+ sys.argv[1] + ".txt" ,'w')
#
#JUMP_DIST=[200]*lons.size
#ls_pts=np.array([lons,lats]).T
#print "NBPTS:", ls_pts.shape[0]
#for i in range(int(sys.argv[1]), min(int(sys.argv[2]), ls_pts.shape[0]) ):
#	lon, lat = ls_pts[i] 
#	cur_dists = get_dist_one_many(lon,lat,lons,lats)
#	cur_vects = get_vect(lon, lat, lons, lats)
#	#print cur_vects
#	#print "test", np.unique(testClassList[0].sortdistIndiv[0:3000])
#	#print "qwr", len(np.argsort(testClassList[0].distargsort)), len(testClassList[0].get_index_pts(i))
#	#print testClassList[0].get_index_pts(i)
#	index_dists = testClassList[0].get_index_pts(i) #[np.argsort(testClassList[0].distargsort)] #np.searchsorted( testClassList[0].sortdistIndiv, cur_dists ) - 1
#	#print np.unique(index_dists)
#	new_dists_forward = cur_dists - np.amin([cur_dists, JUMP_DIST], axis=0) 
#	new_dists_backward = cur_dists + np.amin([cur_dists, JUMP_DIST], axis=0) 
#	#print  np.unique(testClassList[0].sortdistIndiv)
#	index_forward = np.searchsorted(testClassList[0].sortdistIndiv, new_dists_forward)
#	index_backward = np.searchsorted(testClassList[0].sortdistIndiv, new_dists_backward)
#	#print "sortdistIndiv", testClassList[0].sortdistIndiv[0::10000]
#	bug=(index_backward < index_dists) | (index_forward > index_dists)
#	#print "fb_index",zip(index_forward[bug],index_dists[bug],index_backward[bug])
#	#print "fb_values", zip(new_dists_forward[bug],cur_dists[bug] ,new_dists_backward[bug])
#	#print index_forward.shape, np.unique(index_forward)
#	#print testClassList[0].tempname(400,800) 
#	#print np.array([index_dists,index_forward ]).T 	
#	#print testClassList[0].tempname(1000,10000)
#	#print zip(index_dists,index_forward,index_backward)
#	scores_forward = np.array( [ testClassList[0].tempname(a1, a2) for a1,a2 in np.array([index_dists,index_forward ]).T ])
#	scores_backward = np.array( [ testClassList[0].tempname(a1, a2) for a1,a2 in np.array([index_dists,index_backward ]).T ])
#	isBelow = [testClassList[0].sortMaskDistIndiv[a] for a in index_dists ]
#	scores_forward[cur_dists-new_dists_forward==0]=0
#	scores_backward[cur_dists-new_dists_backward==0]=0
#	temp_vect = np.sum((scores_forward - scores_backward) * cur_vects, axis=1)
#	vect=temp_vect/float(len(scores_forward)) #/np.sqrt(np.power(temp_vect[0],2)+ np.power(temp_vect[1],2))
#	#print vect
#	scores_fb = zip(scores_forward, scores_backward)
#	print zip(testClassList[0].sortdistJaccard[index_dists],scores_forward,scores_backward)#print scores_fb 
#	fout.write("\t".join([str(i),str(lon), str(lat), str(vect[0]), str(vect[1])]) +'\t' + "\t".join(map(lambda (x,y): str(x) + ',' + str(y),scores_fb)) + "\n")
#	fout.flush()
#exit()

torad=lambda x : np.array(x)*2*np.pi/360.0
#testParams=[(-88.98,31.0,-84.63,33.83)]
testParams=[(0.0,0.0,0.0,0.0)#empty train
			] 
			#,(148.710938,-5.615986,-97.031250,13.923404)#pacific
			#,(115.3,-27.37,139.2,-26.11) #Australia
			#,(118.13,13.92,178.24,-17.64) #Oceania
			#,(-74.88,-1.75,-4.57,40.7)#spain - america
			#,(21.4453,10.83,109.68,0.0)#africa-oceania
			#]


#testParams=[(-104.765625,15.284185,76.992188,1.054628)]#stupidtrain1

#Initialize all the tests that we are going to do.

#fout = open(os.path.join(OUTPUT_FOLDER, "train_geo_out_" + str(maxtdist)) + ".txt",'w')
fout = open(os.path.join(OUTPUT_FOLDER, "train_geo_out.txt"),'w')

for testParam in testParams:
	computeScores(fout,*testParam)

#exit(1)
#Testing one in particular
#computeScores( fout, *np.append(search_point1, search_point2) )

#searchspace_degree=[]
value=[]
value_param=[]
torad=lambda x : x*2*np.pi/360.0
for search_point1 in [searchspace_degree[i] for i in range(len(searchspace_degree)) if searchspace_valid[i]]:
	for search_point2 in [searchspace_degree[i] for i in range(len(searchspace_degree)) if searchspace_valid[i]]:
		#if is_in_bounds(*search_point1) and is_in_bounds(*search_point2):
		#	continue
		rpt1, rpt2 = torad(search_point1), torad(search_point2)
		length = get_dist_2pt_single(rpt1,rpt2)
		if length <mintdist or length >maxtdist:
			continue
		if search_point1[0] > search_point2[1] or search_point1[1] > search_point2[1]:
			continue
		#print "Length:", length
		computeScores(fout,*np.append(search_point1, search_point2))
		#temp=wrapper_compute_llh(np.append(search_point1, search_point2))
		#value.append(temp)
		#value_param.append(np.append(search_point1, search_point2))
	
fout.close()
print "DONE"
print OUTPUT_FOLDER #This line will is for the reproducible script.

"""
print "Count",count
argtosort=np.argsort(value)
min_value=[value[i] for i in argtosort]
min_value_param=[value_param[i] for i in argtosort]

fout=open(os.path.join(OUTPUT_FOLDER, "train_geo_out") + ".txt",'w')
for item in zip(min_value,min_value_param):
	fout.write(str(item[0])+ ',' + str(item[1].tolist()) + '\n')

fout.close()
"""

















