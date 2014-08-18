#!/usr/bin/env python2
print("Starting train_geo.py")
print "Importing modules...",
import matplotlib
matplotlib.use('Agg')
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
from statsmodels.graphics.boxplots import violinplot
import time
from subprocess import call
#---IMPORTING TESTS-------
from spearmanFlat import spearmanFlat
from mantelTest import mantelTest
from spearmanMantelTest import spearmanMantelTest
from weirdTest import weirdTest
#---END IMPORTING TESTS---
print(" Done")

#---CONSTANT DEFINITION-------
DISTANCE_MASK=20000 #Km
JACCARD_MASK=1.0 #below this number is included
OUTPUT_FOLDER=os.path.join("../results" , time.strftime("%Y.%m.%d.%Hh%mm"))
EARTH_RADIUS=6371
NEGLIGEABLE_RADIAN=0.005
jaccard_file_path="../data/language/Ruhlen2014jaccard.txt"
longlatfile = "../data/language/longitudelatitude_americasnegative.txt"
#---END CONSTANT DEFINITION---

#---PARSE RELATEDNESS-------
call(['mkdir','-p',OUTPUT_FOLDER])
parse_jaccard(open(jaccard_file_path,'r'))
tempDistJaccard = np.array(parse_jaccard(open(jaccard_file_path,'r')))
maskJaccard = (tempDistJaccard<=JACCARD_MASK) 
distJaccard = maskJaccard*tempDistJaccard
applyMaskJaccard= lambda distM : (maskJaccard*distM)

#distJaccard=np.array([[0.,0.5,0.2],[0.5,0.,0.1],[0.2,0.1,0.]])
#---END PARSE RELATEDNESS---

#---PARSING LONG LAT FILE-------
f_longlatfile = open(longlatfile, 'r')
lines_longlatfile = f_longlatfile.read().split('\n')
#pos = [np.array([0.,10.,10.]),np.array([0.,10.,20.])]
pos=map(lambda x : np.array(map(float,x)), zip(*[ line.split('\t') for line in lines_longlatfile ]))
lons = pos[0]#(360+pos[0])*(pos[0]<=0) + pos[0]*(pos[0]>0)
lats = pos[1]#(90-pos[1])*(pos[1]>0) + ( -pos[1] + 90 ) * ( pos[1] <= 0 )
#print "lons:",np.round(lons,1), "lats:", np.round(lats,1)
lats, lons = map(lambda x: x*2*np.pi/360,[lats,lons]) #correct that
#---END PARSING LONG LAT FILE---

#Destription get_spherical_coor: Return points in lons, lats
def get_spherical_coor(points):
	xs, ys, zs = map(np.array,zip(*points))
	lons=np.mod(np.arctan(ys/xs) + (xs<0)*np.pi, 2*np.pi )
	lats=np.arccos(zs/np.sqrt(xs*xs+ys*ys+zs*zs))
	result = np.array(zip(*[ lons*(lons <= np.pi) + (lons-2*np.pi)*(lons>np.pi),(np.pi/2-lats) ]))
	return result

def get_euclidean_coor(lons,lats):
	adjusted_lats = np.copy(lats)
	return EARTH_RADIUS * np.cos(lons) * np.cos(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.cos(adjusted_lats), EARTH_RADIUS*np.sin(adjusted_lats)

#list of points
searchspace_euclidean=np.array(map(lambda x : x.split(" "), open("searchspace3.txt",'r').read().split("\n"))[:-1])#*EARTH_RADIUS
searchspace_euclidean=[map(lambda x : float(x)*EARTH_RADIUS, item) for item in searchspace_euclidean]
#print len(searchspace_euclidean)
#lons, lats
searchspace_radians=np.array(get_spherical_coor(searchspace_euclidean))
searchspace_degree=searchspace_radians * 360.0/(2*np.pi)

#print "eucledean", np.round(searchspace_euclidean,decimals=3)
#print "radians",np.round(searchspace_radians,decimals=3)
#print "degree", np.round(searchspace_degree,decimals=3)

#plot_map_distJaccard(map(list,list(distJaccard)), parse_locations(open("../data/language/latitudelongitude_americasnegative.txt",'r'), flip=True))
#savefig("test.jpg")
def get_dist(lons,lats): # great circle distance.
	matrix_lats, matrix_lons = map(lambda v : np.array([v for i in range(len(v))]) , [lats,lons])
	arg = np.power(np.sin((matrix_lats-np.array([lats]).transpose())/2.0), 2) + np.cos(np.array([lats]).transpose()) * np.cos(matrix_lats) * np.power(np.sin((matrix_lons-np.array([lons]).transpose())/2),2)
	arg = 2*np.arcsin(np.sqrt(arg)) #np.where( np.fabs(arg) < 1., arg, 0.999999)
	return EARTH_RADIUS * arg

def get_dist_to_pt(lon,lat,lons,lats):
	arg = np.power(np.sin((lats-lat)/2.0), 2) + np.cos(lat) * np.cos(lats) * np.power(np.sin((lons-lon)/2), 2)
	arg = 2*np.arcsin(np.sqrt(arg)) #np.where( np.fabs(arg) < 1., arg, 0.999999)
	return EARTH_RADIUS * arg

#expect radian lons lats
def get_dist_2pt(p1,p2):
	#print 'p1',p1,'p2',p2
	return get_dist(*map(np.array,zip(*[p1,p2])))[0,1]


#The arc is from A to B and C is the point. Here A, B are point on sphere(or vector starting at origin of sphere) in euclidean 3D space. C is a vector of the same nature
def get_dist_point_arc(A,B,C):
	#print "cross",A, B, C
	t0= time.time()
	D=np.cross(A,B)
	E=np.cross(D,C)
	U=np.cross(E,D)
	get_norm  = lambda x : np.array(map(lambda y : np.sqrt(np.sum(np.power(y,2))),x))
	get_angle = lambda p1, p2 : np.abs(np.arccos(np.dot(p1,np.array(p2).transpose())/(get_norm(p1)*get_norm(p2)))) #always positive
	U=U/np.array([get_norm(U)]).transpose() * EARTH_RADIUS
	is_on_arc =  np.abs(get_angle([A],[B]) - ((get_angle([A],U)) + (get_angle([B],U)))) < NEGLIGEABLE_RADIAN
	t1=time.time()
	sphere_U, sphere_C, sphere_A, sphere_B = map(lambda x: np.array(get_spherical_coor(x)), [U,C,[A],[B]] )
	#print "part1 time:", t1-t0
	sphere_C_unzipped = np.array(zip(*sphere_C))
	sphere_U_unzipped = np.array(zip(*sphere_U))
	dist_A = get_dist_to_pt(sphere_A[0][0], sphere_A[0][1], sphere_C_unzipped[0], sphere_C_unzipped[1])
	dist_B = get_dist_to_pt(sphere_B[0][0], sphere_B[0][1], sphere_C_unzipped[0], sphere_C_unzipped[1])
	dist_U = get_dist_to_pt(sphere_U_unzipped[0], sphere_U_unzipped[1], sphere_C_unzipped[0], sphere_C_unzipped[1])
	
	temp  = np.where(is_on_arc, dist_U, np.minimum(dist_A, dist_B))
	#temp2 = np.where(is_on_arc, map(lambda x : get_dist_2pt(*x), zip(*[sphere_U,sphere_C])), np.minimum(map(lambda c: get_dist_2pt(sphere_A[0],c),sphere_C), map(lambda c: get_dist_2pt(sphere_B[0],c),sphere_C)))
	#print "equal?", temp[0][0:10], temp2[0][0:10]
	#print "part2 time:", time.time() - t1
	return temp


distIndiv =  get_dist(lons,lats)
pos_euclidean = np.array(get_euclidean_coor(lons,lats))
#xs,ys,zs = pos_euclidean

NB_BINS=10
#print distIndiv
#dist_points_arc=get_dist_point_arc(map(lambda x: x[2], [xs,ys,zs]),map(lambda x: x[0], [xs,ys,zs]), map(lambda x: x[1], [xs,ys,zs]))
#print dist_points_arc

def dist(c,p2):
    return math.hypot(p2[0]-p1[0],p2[1]-p1[1])

def sigmoidf(d,params):
	d=d/10000.0# to prevent math overflow
	#print d,(params[1])*max(d, 0.01) + params[2]
	return params[0]/(1.0 + math.exp((params[1])*max(d, 100/10000.0) + params[2])) #smaller number means more IBD

class line:
	def __init__(self,lon1,lat1,lon2,lat2):
		self.p1=[lon1,lat1]
		self.p2=[lon2,lat2]

mylines=[]
def myrandomrange(a,b):
	return random.random()*(abs(a-b)) + min(a,b)

sigmoid_params=map(float,[1,20.0,0.1])
#sigmoid_params_search=[np.linspace(0.5,2,6),np.linspace(10,25,6),np.linspace(0, 0.5, 6)]
sigmoid_params_search=[[1.0],[20],[0.2]]

#mylines.append(line(0.65,0.17, 0.85, 0.17))
#mylines.append(line(myrandomrange(0.55,1),myrandomrange(0,0.45),myrandomrange(0.55,1),myrandomrange(0,0.45)))

#plot(pos[0], pos[1],".")
#for cline in mylines:
#    plot([cline.p1[0],cline.p2[0]],[cline.p1[1],cline.p2[1]])




#def llh_case(dc,IBDc, nb_pairs_c, sigmoid_params):
#	#if nb_pairs_c==0:
#	#	return 0
#	s_dc=[]
#	for dc_i in dc:
#		s_dc.append(sigmoidf(dc_i,sigmoid_params))
#	s_dc = np.array(s_dc)
#	#return -nb_pairs_c*s_dc + IBDc*np.log(s_dc*nb_pairs_c) - logfact[IBDc])
#	print "test", (nb_pairs_c>0)
#	return np.where(nb_pairs_c>0,-nb_pairs_c*s_dc + IBDc*np.log(s_dc*nb_pairs_c) - scipy.special.gammaln(IBDc+1),0)


#The option fullOutput gives the binary matrix of only the pairs of points that are using the train (ie that are closer to each other when they take the train
def computeDistances(myline,distIndiv,distJaccard,fullOutput=False):
	t,v = map(lambda x: list(get_euclidean_coor(x[0],x[1])), [myline.p1, myline.p2])
	a=get_dist_point_arc(t,v,zip(*pos_euclidean))
	b=np.array( a for i in range(len(distIndiv)))
	trans= a.transpose()
	added= a + trans
	tookTrainMask= (added < distIndiv)
	dist_L = np.minimum(distIndiv , added)
	np.fill_diagonal(dist_L,0)
	if fullOutput:
		return dist_L, tookTrainMask
	else:
		return dist_L


maxdist = DISTANCE_MASK #np.max(b)
bins = np.linspace(0.000001, maxdist, NB_BINS+1) #create NB_BINS bins
dc = (bins + (bins[1]-bins[0])/2.0)[0:-1]
best_min = "nan"
best_sigmoid_params = []


#def llh_bin(dist):
#	bins = np.linspace(0, 1, 10+1)
#	count_per_bin = np.histogram(dist,bins)[0]/2.0
#	log_prob_per_bin = np.where((count_per_bin>0), np.log(count_per_bin/(len(dist)/2.0)), [0,0,0,0,0,0,0,0,0,0] )
#print "count_per_bin", count_per_bin, log_prob_per_bin
#	return np.sum(log_prob_per_bin*count_per_bin)


size_dist_L=np.size(distJaccard)
distJaccard_mean=np.mean(distJaccard)
distJaccard_std=np.std(distJaccard)
np.random.seed(seed=int(time.time()))#random enough for this purpuse.
	
#list_corr=np.sort(spearmanMantelTest(distIndiv,permutated_distJaccard))
#corr=spearmanMantelTest(distIndiv,np.array([distJaccard]))

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
	


testList=[weirdTest]
#testList=[mantelTest,spearmanMantelTest,spearmanFlat]
testClassList=[thistest(applyMaskJaccard(distIndiv),applyMaskJaccard(distJaccard)) for thistest in testList]
def computeScores(fout,lon1,lat1,lon2,lat2):
	t0=time.time()
	fout.write("Parameters:%f,%f,%f,%f\n" % (lon1,lat1,lon2,lat2))
	lon1,lat1,lon2,lat2 = [x*2*np.pi/360.0 for x in [lon1,lat1,lon2,lat2]]
	length= get_dist_2pt([lon1,lat1] , [lon2,lat2])
	fout.write("Length of the line: " + str(length)+ '\n')
	myLine=line(lon1,lat1,lon2,lat2)
	bestDistances,tookTrainMask=computeDistances(myLine,distIndiv,distJaccard,fullOutput=True)
	testResults=[]
	

	for testClass in testClassList:
		print "Begin: " , str(testClass).split(' ')[0].split('.')[1] , "..." ,  
		testResults.append(testClass(myLine,applyMaskJaccard(bestDistances),tookTrainMask))
		print "Done"
	for testResultIndex in range(len(testResults)):
		fout.write("Test "+str(testClassList[testResultIndex]).split(' ')[0].split('.')[1] + ": " + str(testResults[testResultIndex]) + '\n')
	geoDistVsJaccardDist(lon1,lat1,lon2,lat2,bestDistances,testResults[0],tookTrainMask)
	print "Time for a whole iteration: " , time.time()-t0
	return testResults

testParams=[(148.710938,-5.615986,-97.031250,13.923404)#pacific #(0.0,0.0,0.0,0.0)#empty train
			,(115.3,-27.37,139.2,-26.11) #Australia
			,(118.13,13.92,178.24,-17.64) #Oceania
			,(-74.88,-1.75,-4.57,40.7)#spain - america
			,(21.4453,10.83,109.68,0.0)#africa-oceania
			]


#testParams=[(-104.765625,15.284185,76.992188,1.054628)]#stupidtrain1

#Initialize all the tests that we are going to do.

fout = open(os.path.join(OUTPUT_FOLDER, "train_geo_out") + ".txt",'w')

for testParam in testParams:
	computeScores(fout,*testParam)


#searchspace_degree=[]
value=[]
value_param=[]
for search_point1 in searchspace_degree:
	for search_point2 in searchspace_degree:
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

















