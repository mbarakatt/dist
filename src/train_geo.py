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
import random,math,scipy,sys
import numpy as np
np.set_printoptions(threshold=np.nan)
from itertools import product
from mpl_toolkits.basemap import Basemap, shiftgrid
from jaccard_plot import *
from statsmodels.graphics.boxplots import violinplot
import time
print(" Done")


OUTPUT_FOLDER="../results" + time.strftime("%Y.%m.%d.%Hh%mm")


EARTH_RADIUS=6371
NEGLIGEABLE_RADIAN=0.005
jaccard_file_path="../data/language/Ruhlen2014jaccard.txt"
parse_jaccard(open(jaccard_file_path,'r'))

relatedness = np.array(parse_jaccard(open(jaccard_file_path,'r')))
#relatedness=np.array([[0.,0.5,0.2],[0.5,0.,0.1],[0.2,0.1,0.]])

longlatfile = "../data/language/longitudelatitude_americasnegative.txt"
f_longlatfile = open(longlatfile, 'r')
lines_longlatfile = f_longlatfile.read().split('\n')
#pos = [np.array([0.,10.,10.]),np.array([0.,10.,20.])]
pos=map(lambda x : np.array(map(float,x)), zip(*[ line.split('\t') for line in lines_longlatfile ]))

lons = pos[0]#(360+pos[0])*(pos[0]<=0) + pos[0]*(pos[0]>0)
lats = pos[1]#(90-pos[1])*(pos[1]>0) + ( -pos[1] + 90 ) * ( pos[1] <= 0 )
#print "lons:",np.round(lons,1), "lats:", np.round(lats,1)
lats, lons = map(lambda x: x*2*np.pi/360,[lats,lons]) #correct that


def get_spherical_coor(points):#return points in lons, lats
	xs, ys, zs = map(np.array,zip(*points))
	lons=np.mod(np.arctan(ys/xs) + (xs<0)*np.pi, 2*np.pi )
	lats=np.arccos(zs/np.sqrt(xs*xs+ys*ys+zs*zs))
	result = np.array(zip(*[ lons*(lons <= np.pi) + (lons-2*np.pi)*(lons>np.pi),(np.pi/2-lats) ]))
	return result

def get_euclidean_coor(lons,lats):
	adjusted_lats = np.copy(lats)
	#return EARTH_RADIUS * np.cos(lons) * np.sin(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.sin(adjusted_lats), EARTH_RADIUS*np.cos(adjusted_lats)
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

#plot_map_relatedness(map(list,list(relatedness)), parse_locations(open("../data/language/latitudelongitude_americasnegative.txt",'r'), flip=True))
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

"""
torad=lambda x : x*2*np.pi/360.0
p1=map(torad,[142  ,  31.8 ])
p2=map(torad,[18.45, -16.29])
print "THIS" , get_dist_2pt(p1,p2)
print "ANDTHIS",get_dist_to_pt(p1[0], p1[1], np.array([p1[0],p2[0]]), np.array([p1[1],p2[1]]))
exit(1)
"""




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


dist_indiv    = get_dist(lons,lats)
pos_euclidean = np.array(get_euclidean_coor(lons,lats))
#xs,ys,zs = pos_euclidean

NB_BINS=10
#print dist_indiv
#dist_points_arc=get_dist_point_arc(map(lambda x: x[2], [xs,ys,zs]),map(lambda x: x[0], [xs,ys,zs]), map(lambda x: x[1], [xs,ys,zs]))
#print dist_points_arc

def dist(p1,p2):
    return math.hypot(p2[0]-p1[0],p2[1]-p1[1])

def sigmoidf(d,params):
	d=d/10000.0# to prevent math overflow
	#print d,(params[1])*max(d, 0.01) + params[2]
	return params[0]/(1.0 + math.exp((params[1])*max(d, 100/10000.0) + params[2])) #smaller number means more IBD

class line:
	def __init__(self,lon1,lat1,lon2,lat2):
		self.p1=[lon1,lat1]
		self.p2=[lon2,lat2]
		#print lon1,lat1
		#toradians= lambda x : 2*np.pi/360
		#self.p1=[toradians(lon1),toradians(lat1)]
		#self.p2=[toradians(lon2),toradians(lat2)]
		#self.p1, self.p2=map(lambda x : map(lambda y : y*2*np.pi/360.0,x) ,[[float(lon1),float(lat1)],[float(lon2),float(lat2)]])

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




"""
def llh_case(dc,IBDc, nb_pairs_c, sigmoid_params):
	#if nb_pairs_c==0:
	#	return 0
	s_dc=[]
	for dc_i in dc:
		s_dc.append(sigmoidf(dc_i,sigmoid_params))
	s_dc = np.array(s_dc)
	#return -nb_pairs_c*s_dc + IBDc*np.log(s_dc*nb_pairs_c) - logfact[IBDc])
	print "test", (nb_pairs_c>0)
	return np.where(nb_pairs_c>0,-nb_pairs_c*s_dc + IBDc*np.log(s_dc*nb_pairs_c) - scipy.special.gammaln(IBDc+1),0)
"""

def compute_distances(mylines):
	#print 'myline',mylines[0].p1, mylines[0].p2
	t,v = map(lambda x: list(get_euclidean_coor(x[0],x[1])), [mylines[0].p1, mylines[0].p2])
	#print "tv", t,v
	#print "dist_indiv", dist_indiv
	a=get_dist_point_arc(t,v,zip(*pos_euclidean))
	#print "pt_arc_dist",a
	#print "arc_time:", time.time() - t0
	b=np.array( a for i in range(len(dist_indiv)))
	trans= a.transpose()
	added= a + trans
	#dist_L = np.minimum(trans, dist_indiv)
	dist_L = np.minimum(dist_indiv , added)
	np.fill_diagonal(dist_L,0)
	#print "min_dist_L", dist_L
	return dist_L

maxdist = 20000.0 #np.max(b)
bins = np.linspace(0, maxdist, NB_BINS+1) #create NB_BINS bins
dc = (bins + (bins[1]-bins[0])/2.0)[0:-1]
best_min = "nan"
best_sigmoid_params = []


def llh_bin(dist):
	bins = np.linspace(0, 1, 10+1)
	count_per_bin = np.histogram(dist,bins)[0]/2.0
	log_prob_per_bin = np.where((count_per_bin>0), np.log(count_per_bin/(len(dist)/2.0)), [0,0,0,0,0,0,0,0,0,0] )
	#print "count_per_bin", count_per_bin, log_prob_per_bin
	return np.sum(log_prob_per_bin*count_per_bin)


size_dist_L=np.size(relatedness)
relatedness_mean=np.mean(relatedness)
relatedness_std=np.std(relatedness)
np.random.seed(seed=int(time.time()))

print "Creating permutations...",
permutated_relatedness= np.array([np.random.permutation(relatedness)for _ in range(49)])
permutated_relatedness_mean=np.mean(permutated_relatedness)
permutated_relatedness_std=np.std(permutated_relatedness)
print "Done"

#M1 is a single matrix, M2 can be a array of matrix
def mantel_test(M1,M2):
	t1=time.time()
	results=1.0/(size_dist_L - 1) * np.sum((M1 - np.mean(M1))/np.std(M1)*(M2 - permutated_relatedness_mean)/permutated_relatedness_std ,axis=(1,2))
	print "Mantel test time:", time.time() - t1
	return -results
	
list_corr=np.sort(mantel_test(dist_indiv,permutated_relatedness))
corr=mantel_test(dist_indiv,np.array([relatedness]))

def find_index_in_array(array, c):
	for i in range(np.size(array)):
		if c <= array[i]:
			return i
	return len(array)

print (50-find_index_in_array(list_corr,corr[0]))/50

print list_corr
print corr

def compute_llh(dist_L, relatedness, plotviolin=True):
	#results=1.0/(size_dist_L-1) * np.sum((dist_L - np.mean(dist_L))/np.std(dist_L)*(relatedness - relatedness_mean)/relatedness_std)
	list_corr1=np.sort(mantel_test(dist_L,permutated_relatedness))
	corr1=mantel_test(dist_L,np.array([relatedness]))
	#print list
	#print corr1
	#plot([0]*50+[0.1]*50,list_corr.tolist()+[corr]+list_corr1.tolist()+[corr1],'.')
	#xlim(-0.05,0.15)
	#title(str(corr-np.mean(list_corr))+ ',' + str( corr1 - np.mean(list_corr1) ) )
	#plt.show()
	#print (50-find_index_in_array(list_corr,corr[0]))/50
	return (corr1 - np.mean(list_corr1))[0]


"""
	#t0=time.time()
	a = np.ndarray.flatten(relatedness)
	b = np.ndarray.flatten(dist_L)
	#t1=time.time()
	#print "flattening time: ", t1 -t0
	c=[]
	
	for i in range(len(bins)-1):
		c.append(a[(b > bins[i]) * (b<=bins[i+1])])
		
	#map(lambda (x,y): c[int(y/maxdist*NB_BINS)].append(x), zip(a,b))
	#t2=time.time()
	#print "get relatedness in bin time:", t2-t1
	if plotviolin:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		violinplot(c, ax=ax)
		plt.show()
		
		plot(dc, np.array(map(np.sum,c))/np.array(map(len,c)))
		plt.show()

	temp=map(llh_bin, c)
	#print "bin_llh",temp
	sum_llh= -np.sum(temp)

	return sum_llh
"""
count=0

def wrapper_compute_llh(params):
	t0=time.time()
	if params[0] > params[2] or params[1] > params[3]:
		return np.inf
		
	print "Computing llh for", params
	params=map(lambda x: x*2*np.pi/360.0, params)
	length= get_dist_2pt( [params[0], params[1]],[ params[2] , params[3]])
	print "LENGTH:", length
	if  length < 300 or length > 6000:
		return np.inf
	a=compute_llh(compute_distances(fittedlines + [line(params[0], params[1], params[2], params[3])]), relatedness)
	print "llh: ", a
	#print "full_iter_time", time.time() - t0
	return a

dist_Lines = [[0]] #distance_between_lines(mylines)

fittedlines=[]
for i in [0]:
	test_params=[(115.3,-27.37,139.2,-26.11) #Australia
				,(118.13,13.92,178.24,-17.64) #Oceania
				,(-74.88,-1.75,-4.57,40.7)#spain - america
				,(21.4453,10.83,109.68,0.0)#africa-oceania
				]
	
	
	for these_params in test_params:
		print these_params
		llh=compute_llh(compute_distances(fittedlines + [line(*map(lambda x : x*2*np.pi/360,these_params))]), relatedness,plotviolin=False)
		print llh
	value=[]
	value_param=[]
	for search_point1 in searchspace_degree:
		for search_point2 in searchspace_degree:
			temp=wrapper_compute_llh(np.append(search_point1, search_point2))
			value.append(temp)
			value_param.append(np.append(search_point1, search_point2))
		
	print "Count",count
	argtosort=np.argsort(value)
	min_value=[value[i] for i in argtosort]
	min_value_param=[value_param[i] for i in argtosort]

	fout=open(os.path.join(OUTPUT_FOLDER, "train_geo_out") +time.strftime("%Hh%Mm")+".txt",'w')
	for item in zip(min_value,min_value_param):
		fout.write(str(item[0])+ ',' + str(item[1].tolist()) + '\n')

	fout.close()


#This line will is for the reproducible program.
print OUTPUT_FOLDER

	"""
	#res, fval, grid, Jout=scipy.optimize.brute(wrapper_compute_llh, [(mybox[0],mybox[2]),(mybox[1],mybox[3]),(mybox[0],mybox[2]),(mybox[1],mybox[3])], Ns=i , full_output=True, finish=None)
	print "result", res
	print "fval", fval
	bestline=line(res[0],res[1],res[2],res[3])
	fittedlines.append(bestline)
	for cline in fittedlines:
		plot([cline.p1[0],cline.p2[0]],[cline.p1[1],cline.p2[1]], 'o')
		plot([cline.p1[0],cline.p2[0]], [cline.p1[1],cline.p2[1]])
		#plot(pos[0], pos[1],".")
		
	sys.exit(1)


	title("Fitted line:" + str(compute_corr(compute_distances(pos, dist_bet_indiv, fittedlines ,np.zeros([len(fittedlines)+1,len(fittedlines)+ 1])),IBD)) + ", " +"Real line:" + str(compute_corr(compute_distances(pos, dist_bet_indiv, mylines ,np.zeros([len(fittedlines)+1,len(fittedlines)+ 1])),IBD))+ "Params:" +str(best_sigmoid_params))
	savefig("pic/test" + str(args[1]) + ".jpg")

	fig=figure()
	ax = fig.add_subplot(1,1,1)
	drawrefsigm(ax)
	drawsigm(compute_distances(pos, dist_bet_indiv, fittedlines, np.zeros([len(fittedlines)+1,len(fittedlines)+ 1])),IBD,"Fitted",ax)
	drawsigm(compute_distances(pos, dist_bet_indiv, mylines, np.zeros([len(fittedlines)+1,len(fittedlines)+ 1])),IBD,"Real",ax)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles[::-1], labels[::-1]) #reverse the order
	title("Fitted line:" + str(compute_corr(compute_distances(pos, dist_bet_indiv, fittedlines ,np.zeros([len(fittedlines)+1,len(fittedlines)+ 1])),IBD)) + ", " +"Real line:" + str(compute_corr(compute_distances(pos, dist_bet_indiv, mylines ,np.zeros([len(fittedlines)+1,len(fittedlines)+ 1])),IBD))+"Params:"+str(best_sigmoid_params))
	savefig("pic/" + "test" + str(args[1]) + "prob"  + ".jpg")


	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for x in np.linspace(mybox[0], mybox[2],num=40):
		for y in np.linspace(mybox[1],mybox[3],num=40):
			ax.scatter(x,y,wrapper_compute_corr([mylines[0].p1[0], mylines[0].p1[1], x, y]))
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	for ii in xrange(0,360,30):
		ax.view_init(elev=30., azim=ii)
		plt.savefig("landscape" + str(ii) + ".png")
	"""

print "DONE"
#sys.exit(1)

















