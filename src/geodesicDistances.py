import numpy as np
from geode import *
from time import *
import itertools

NEGLIGEABLE_RADIAN=0.005
longlatfile = "../data/language/longitudelatitude_americasnegative.txt"
outFile="geographicGeoDistLanguages.txt"
fOutFile=open(outFile,'w')


def get_dist_one_many(lon1,lat1,lons,lats): # great circle distance.
	#haversin = lambda theta : np.power(np.sin(theta/2.0),2) #(1-np.cos(theta))/2.0
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
#pos = [np.array([0.,10.,10.]),np.array([0.,10.,20.])]
pos=map(lambda x : np.array(map(float,x)), zip(*[ line.split('\t') for line in lines_longlatfile ]))
lons = pos[0]#(360+pos[0])*(pos[0]<=0) + pos[0]*(pos[0]>0)
lats = pos[1]#(90-pos[1])*(pos[1]>0) + ( -pos[1] + 90 ) * ( pos[1] <= 0 )
#print "lons:",np.round(lons,1), "lats:", np.round(lats,1)
lats, lons = map(lambda x: x*2*np.pi/360,[lats,lons]) #correct that
#---END PARSING LONG LAT FILE---

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


print "Test"

dist=[]#np.zeros(len(lons))
posZipped=zip(lons,lats)
print posZipped[0][1]
for pos1 in range(len(posZipped)):
	t0=time()
	#print "Current:", lons[pos1],lats[pos1],lons,lats
	dist.append(np.concatenate([[0]*pos1,computeDistancesLanguages(posZipped[pos1][0],posZipped[pos1][1],lons[pos1:],lats[pos1:])]))
#	print "This iter:", computeDistancesLanguages(posZipped[pos1][0],posZipped[pos1][1],lons,lats)
	print "pos1 ", pos1, "Time:" , time()-t0

print np.array(dist)
	
for i in range(len(lons)):
	for j in range(len(lons)):
		fOutFile.write(str(dist[i][j] ))
		if j != len(lons) - 1 :
			fOutFile.write("\t")
		else:
			fOutFile.write("\n")


