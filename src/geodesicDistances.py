import numpy as np
from geode import *
from time import *
import itertools
import sys
import getopt

NEGLIGEABLE_RADIAN=0.005
#longlatfile = "../data/language/longitudelatitude_americasnegative.txt"
#outFile="geographicGeoDistLanguages.txt"
def usage():
	print "Usage: \n \t -l longlatfile.txt \n \t -o outfile.txt"
	print "You entered:" ,sys.argv
optlist,args=getopt.getopt(sys.argv[1:],"o:l:")
outFile=""
longlatfile=""
for flag,val in optlist:
	if flag =="-o":
		outFile = val
		outFileWOTrain=outFile.split('.')[0] + "WO.txt"
	if flag =="-l":
		longlatfile = val 
if outFile=="":
	outFile="geographicGeoDistLanguagesTest.txt"
	outFileWOTrain=outFile.split('.')[0] +"WO.txt" 
if longlatfile=="":
	usage()
	exit(2)

fOutFile=open(outFile,'w')
fOutFileWO=open(outFileWOTrain,'w')

class line:
	def __init__(self,lon1,lat1,lon2,lat2):
		self.p1=[lon1,lat1]
		self.p2=[lon2,lat2]

#Input are numpy arrAY OF PT IN RADIANS(LONS,LATS)
def get_dist_2pt_array(ps1,ps2):
	ps1_lons, ps1_lats = ps1.T
	ps2_lons, ps2_lats = ps2.T #map(np.array, zip(*ps2))
	return get_dist_one_many(ps1_lons, ps1_lats, ps2_lons, ps2_lats)

def get_dist_point_arc(A,B,C):
	D=np.cross(A,B)
	E=np.cross(D,C)
	U=np.cross(E,D)
	U=U/np.array([get_norm(U)]).transpose() * EARTH_RADIUS
	is_on_arc =  np.abs(get_angle([A],[B]) - ((get_angle([A],U)) + (get_angle([B],U)))) < NEGLIGEABLE_RADIAN
	#t1=time.time()
	sphere_U, sphere_C, sphere_A, sphere_B = map(lambda x: np.array(get_spherical_coor(x)), [U,C,[A],[B]] )
	sphere_C_unzipped = np.array(zip(*sphere_C))
	sphere_U_unzipped = np.array(zip(*sphere_U))
	dist_A = get_dist_2pt(sphere_A[0], sphere_C)
	dist_B = get_dist_2pt(sphere_B[0], sphere_C)
	dist_U = get_dist_2pt_array(sphere_U, sphere_C)
	temp  = np.where(is_on_arc, dist_U, np.minimum(dist_A, dist_B))
	return temp

def computeDistances(myline,distIndiv):
	t,v = map(lambda x: list(get_euclidean_coor(x[0],x[1])), [myline.p1, myline.p2])
	a=get_dist_point_arc(t,v,pos_euclidean)
	b=np.array( [a] * len(distIndiv))
	added= b.T + b
	tookTrainMask= (added < distIndiv)
	dist_L = np.minimum(distIndiv , added)
	np.fill_diagonal(dist_L,0)
	return dist_L 

def get_dist_one_many(lon1,lat1,lons,lats): # great circle distance.
	#haversin = lambda theta : np.power(np.sin(theta/2.0),2) #(1-np.cos(theta))/2.0
	#print lon1, lat1, lons, lats, np.cos(lats), np.cos(lats), haversin(lons-lon1), lons-lon1,haversin(lats-lat1), lats-lat1
	lons=np.array(lons)
	arg = haversin(lats-lat1) + np.cos(lats) * np.cos(lat1) * haversin( lons - lon1 )
	arg = 2 * EARTH_RADIUS * np.arcsin(np.sqrt(arg)) 
	return arg

def get_dist_to_pt(lon,lat,lons,lats):
	arg = np.power(np.sin((lats-lat)/2.0), 2) + np.cos(lat) * np.cos(lats) * np.power(np.sin((lons-lon)/2), 2)
	arg = 2*np.arcsin(np.sqrt(arg)) #np.where( np.fabs(arg) < 1., arg, 0.999999)
	return EARTH_RADIUS * arg

def pointIsOnArc(ps,p1,ps2):
#	print "pointIsOnArc:DATA:", ps,p1,ps2
	angle2pt=get_angle([p1],ps2)
	angleothertwo = (get_angle([p1],ps)) + (get_angle(ps,ps2))
	#print "angles",angle2pt,angleothertwo
	#print "pointIsOnArc:", map(lambda x : x*360/(2*np.pi) , [angle2pt,angleothertwo]),np.abs(angle2pt + angleothertwo) 
	#print "pointIsOnArc;",  np.abs(angle2pt - angleothertwo) < NEGLIGEABLE_RADIAN , np.abs(angle2pt + angleothertwo - 2*np.pi) < NEGLIGEABLE_RADIAN

	a= (np.abs(angle2pt - angleothertwo) < NEGLIGEABLE_RADIAN) , (np.abs(angle2pt + angleothertwo-2*np.pi) < NEGLIGEABLE_RADIAN)
	#print "a",a
	return a

def arcIntersect(lon1, lat1, lons2, lats2, lon3, lat3, lon4, lat4):
	p1, p3, p4 = map(lambda (x,y): get_euclidean_coor(x,y), [[lon1,lat1],[lon3,lat3],[lon4,lat4]])
	ps2 = get_euclidean_coor(lons2,lats2)
	#print "arcIntersect", p1,ps2,p3,p4
	#print 'euclidean ps2', ps2
	perps1, perp2=[np.cross(p1,ps2),np.cross(p3,p4)]
	#print "arcIntersect:",perps1,perp2
	#print p1, ps2
	intersections = np.cross(perps1, perp2)
	#print "arcIntersect:", intersections
	#print "arcIntersect:",  pointIsOnArc(intersections,p1,ps2) , pointIsOnArc(intersections,p3,[p4])

	arc1 = pointIsOnArc(intersections,p1,ps2)
	arc2 = pointIsOnArc(intersections,p3,[p4])
	return (arc1[0] & arc2[0] ) | (arc1[1] & arc2[1])
	#return  pointIsOnArc(intersections,p3,[p4])

#radian lon lat
def arcIntersectProhibited(lon1,lat1,lons2,lats2):
#	print 'arcIntersectProhibited', lon1,lat1,lons2,lats2
	temp=np.array([False]*len(lons2))
	for blockLine in blockLines:
		temp = temp | arcIntersect(lon1,lat1,lons2,lats2,blockLine[0],blockLine[1],blockLine[2],blockLine[3]) #np.concatenate([p1,ps2,blockLine]))
		#print "temp",temp
	return temp

def computeDistancesLanguages(lon1,lat1,lons2,lats2):
	p1=[lon1,lat1]
	ps2=np.array([lons2,lats2]).T
	temp=[get_dist_2pt(p1,ps2)]
	for wayPoint1, wayPoint2 in itertools.product(range(len(wayPoints)),range(len(wayPoints))):
		if True: #wayPoint1 >= wayPoint2:
			#print "computeDistancesLanguages",p1,ps2
			temp.append( get_dist_2pt_single(p1,wayPoints[wayPoint1]) + distanceWayPoints[wayPoint1,wayPoint2] + get_dist_2pt(wayPoints[wayPoint2],ps2))
	#print "Chosing for amin", np.array(temp)
	return np.amin(temp,axis=0)


def get_dist_2pt_single(p1,p2):
	retval= get_dist_2pt(p1,np.array([p2]))[0]
	return retval

#Input are numpy array of pt in radians(lons,lats)
def get_dist_2pt(p1,ps2):
	ps2_lons, ps2_lats = map(np.array, zip(*ps2))
	return get_dist_one_many(p1[0], p1[1], ps2_lons, ps2_lats)  + 1000000*arcIntersectProhibited(p1[0],p1[1],ps2_lons,ps2_lats)


#---PARSING LONG LAT FILE-------
f_longlatfile = open(longlatfile, 'r')
lines_longlatfile = f_longlatfile.read().split('\n')
#pos = [np.array([0.,10.,10.]),np.array([0.,10.,20.])]print "lines",Lines_longlatfile
if lines_longlatfile[-1]=='':
	lines_longlatfile=lines_longlatfile[:-1]
pos= map(np.array,np.array([map(float,this_line.split('\t')) for this_line in lines_longlatfile ]).T)
lons = pos[0]#(360+pos[0])*(pos[0]<=0) + pos[0]*(pos[0]>0)
lats = pos[1]#(90-pos[1])*(pos[1]>0) + ( -pos[1] + 90 ) * ( pos[1] <= 0 )
#print "lons:",np.round(lons,1), "lats:", np.round(lats,1)
lats, lons = map(lambda x: x*2*np.pi/360,[lats,lons]) #correct that
#---END PARSING LONG LAT FILE---

pos_euclidean=get_euclidean_coor(lons,lats)
#---DEFINING BLOCKING LINES---
#lon1,lat1,lon2,lat2
toRadians = lambda x : x*2*np.pi/360.0

if False:
	blockLines=toRadians(np.array([[-173.320313,54.534318,-174.023438,-57.136239],#Pacific
				[-28.476563,62.267923,-26.367188,-52.268157],#Atlantic
				[62.929688,-60.759160,62.929688,11.141932],#Indian
				[26.0,33.5,11.6,36.1],#Mediterranean 
				[-9.5,34.3,11.3,39.1],#Strait of Gibraltar
				[-126.38, 77.31, 123.4, 77.88]#Above Russia/Canada
				]))
else:
	blockLines=np.array([])
#---END DEFINING BLOCKING LINES-


#---WAYPOINTS---
wayPoints=map(lambda (x,y): (y*2*np.pi/360.0,x*2*np.pi/360.0) ,[(64,177),
(30,31),
(41,28),
(11,104),
(54,-130)])
triangularToFull = lambda x : x+x.T 
distanceWayPoints = np.array([[       0.   ,9147.0055493 ,7862.376 , 8084.30063025,   3123.47727646],
 [    9147.0055493      ,   0.      ,       1252.69961064  ,   7770.52999011,  12270.48282577],
 [ 7862.376   ,  1252.69961064   ,     0.       ,      8036.88417597,  11146.414],
 [    8084.30063025  ,   7770.52999011   ,  8036.88417597     ,   0.,   11207.77790671],
 [    3123.47727646  ,  12270.48282577,  11146.414,    11207.77790671,         0.        ]])
#distanceWayPoints = triangularToFull(np.array([[0,get_dist_2pt_single(wayPoints[0],wayPoints[1]),get_dist_2pt_single(wayPoints[0],wayPoints[2]),get_dist_2pt_single(wayPoints[0],wayPoints[3]),get_dist_2pt_single(wayPoints[0],wayPoints[4])],
#			[0,0, get_dist_2pt_single(wayPoints[1],wayPoints[2]) , get_dist_2pt_single(wayPoints[1],wayPoints[3]), get_dist_2pt_single(wayPoints[1],wayPoints[0]) + get_dist_2pt_single(wayPoints[0],wayPoints[4])],
#			[0,0,0,get_dist_2pt_single(wayPoints[2],wayPoints[3]),get_dist_2pt_single(wayPoints[0],wayPoints[2]) + get_dist_2pt_single(wayPoints[0],wayPoints[4])],
#			[0,0,0,0,get_dist_2pt_single(wayPoints[3],wayPoints[0]) + get_dist_2pt_single(wayPoints[0],wayPoints[4])],
#			[0,0,0,0,0]]))

print distanceWayPoints
#---END WAYPOINTS

#print arcIntersectProhibited(toRadians(-30),toRadians(0),np.array([toRadians(10)]),np.array([toRadians(0)]))



dist=[]#np.zeros(len(lons))
posZipped=zip(lons,lats)
print posZipped[0][1]
for pos1 in range(len(posZipped)):
	#print "Current:", lons[pos1],lats[pos1],lons,lats
	dist.append(np.concatenate([[0]*pos1,computeDistancesLanguages(posZipped[pos1][0],posZipped[pos1][1],lons[pos1:],lats[pos1:])]))
#	print "This iter:", computeDistancesLanguages(posZipped[pos1][0],posZipped[pos1][1],lons,lats)
	print "pos1 ", pos1, "Time:" #, #time()-t0



dist = np.triu(np.array(dist),k=1)
dist = dist + dist.T

bestDistances=computeDistances(line(0,0,40,40), dist)
np.savetxt(fOutFile,bestDistances,delimiter="\t")	
np.savetxt(fOutFileWO,dist,delimiter='\t')
exit(1)
for i in range(len(lons)):
	for j in range(len(lons)):
		fOutFile.write(str(dist[i][j] ))
		if j != len(lons) - 1 :
			fOutFile.write("\t")
		else:
			fOutFile.write("\n")


