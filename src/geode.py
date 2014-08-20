import numpy as np

NEGLIGEABLE_RADIAN=0.005
EARTH_RADIUS=6371



def get_spherical_coor(points):
	xs, ys, zs = map(np.array,zip(*points))
	lons=np.mod(np.arctan(ys/xs) + (xs<0)*np.pi, 2*np.pi )
	lats=np.arccos(zs/np.sqrt(xs*xs+ys*ys+zs*zs))
	result = np.array(zip(*[ lons*(lons <= np.pi) + (lons-2*np.pi)*(lons>np.pi),(np.pi/2-lats) ]))
	return result

def get_euclidean_coor(lons,lats):
	adjusted_lats = np.copy(lats)
	return EARTH_RADIUS * np.cos(lons) * np.cos(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.cos(adjusted_lats), EARTH_RADIUS*np.sin(adjusted_lats)

def get_dist(lons,lats): # great circle distance.
	matrix_lats, matrix_lons = map(lambda v : np.array([v for i in range(len(v))]) , [lats,lons])
	arg = np.power(np.sin((matrix_lats-np.array([lats]).transpose())/2.0), 2) + np.cos(np.array([lats]).transpose()) * np.cos(matrix_lats) * np.power(np.sin((matrix_lons-np.array([lons]).transpose())/2),2)
	arg = 2*np.arcsin(np.sqrt(arg)) #np.where( np.fabs(arg) < 1., arg, 0.999999)
	return EARTH_RADIUS * arg

def get_dist_to_pt(lon,lat,lons,lats):
	arg = np.power(np.sin((lats-lat)/2.0), 2) + np.cos(lat) * np.cos(lats) * np.power(np.sin((lons-lon)/2), 2)
	arg = 2*np.arcsin(np.sqrt(arg)) #np.where( np.fabs(arg) < 1., arg, 0.999999)
	return EARTH_RADIUS * arg

def pointIsOnArc(p,p1,p2):
	return np.abs(get_angle([p1],[p2]) - np.mod((get_angle([p1],[p])) + (get_angle([p],[p2])), np.pi/2)) < NEGLIGEABLE_RADIAN

def arcIntersect(lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4):
	p1,p2,p3,p4= map(lambda (x,y): get_euclidean_coor(x,y), [[lon1,lat1],[lon2,lat2],[lon3,lat3],[lon4,lat4]])
	perp1,perp2=[np.cross(p1,p2),np.cross(p3,p4)]
	intersection=np.cross(perp1,perp2)
	return pointIsOnArc(intersection,p1,p2) and pointIsOnArc(intersection,p3,p4)

#radian lon lat
def arcIntersectProhibited(p1,p2):
	temp=False
	for blockLine in blockLines:
		temp = temp or arcIntersect(*np.concatenate([p1,p2,blockLine]))
	return temp

def computeDistancesLanguages(lon1,lat1,lon2,lat2):
	p1=[lon1,lat1]
	p2=[lon2,lat2]
	temp=np.array(get_dist_2pt(p1,p2))
	for wayPoint1 in range(len(wayPoints)):
		for wayPoint2 in range(len(wayPoints)):
			if wayPoint1 <= wayPoint2:
				#print temp,(get_dist_2pt(p1,wayPoints[wayPoint1]) +distanceWayPoints[wayPoint1,wayPoint2]+ get_dist_2pt(wayPoints[wayPoint2],p2))[0]
				temp=np.concatenate([temp.flatten(),( get_dist_2pt(p1,wayPoints[wayPoint1]) +distanceWayPoints[wayPoint1,wayPoint2]+ get_dist_2pt(wayPoints[wayPoint2],p2))[0]])
	return min(temp)

#expect radian lons lats
def get_dist_2pt(p1,p2):
	#print 'p1',p1,'p2',p2
	return get_dist(*map(np.array,zip(*[p1,p2])))[0,1] + np.nan_to_num(np.inf*arcIntersectProhibited(p1,p2))

def get_norm(x):
	return  np.array(map(lambda y : np.sqrt(np.sum(np.power(y,2))),x))
def get_angle(p1,p2):
	return  np.abs(np.arccos(np.dot(p1,np.array(p2).transpose())/(get_norm(p1)*get_norm(p2)))) #always positive

#---DEFINING BLOCKING LINES---
#lon1,lat1,lon2,lat2
toRadians = lambda x : x*2*np.pi/360
blockLines=toRadians(np.array([[-173.320313,54.534318,-174.023438,-57.136239],#Pacific
		[-28.476563,62.267923,-26.367188,-52.268157],#Atlantic
		[62.929688,-60.759160,62.929688,10.141932],#Indian
		]))
#---END DEFINING BLOCKING LINES-

#---WAYPOINTS---
wayPoints=map(lambda (x,y): (y*2*np.pi/360.0,x*2*np.pi/360.0) ,[(64,177),
(30,31),
(41,28),
(11,104),
(54,-130)])
triangularToFull = lambda x : x+x.T 
distanceWayPoints = triangularToFull(np.array([[0,get_dist_2pt(wayPoints[0],wayPoints[1]),get_dist_2pt(wayPoints[0],wayPoints[2]),get_dist_2pt(wayPoints[0],wayPoints[3]),get_dist_2pt(wayPoints[0],wayPoints[4])],
			[0,0, get_dist_2pt(wayPoints[1],wayPoints[2]) , get_dist_2pt(wayPoints[1],wayPoints[3]), get_dist_2pt(wayPoints[1],wayPoints[0]) + get_dist_2pt(wayPoints[0],wayPoints[4])],
			[0,0,0,get_dist_2pt(wayPoints[2],wayPoints[3]),get_dist_2pt(wayPoints[0],wayPoints[2]) + get_dist_2pt(wayPoints[0],wayPoints[4])],
			[0,0,0,0,get_dist_2pt(wayPoints[3],wayPoints[0]) + get_dist_2pt(wayPoints[0],wayPoints[4])],
			[0,0,0,0,0]]))

print distanceWayPoints
#---END WAYPOINTS

#The arc is from A to B and C is the point. Here A, B are point on sphere(or vector starting at origin of sphere) in euclidean 3D space. C is a vector of the same nature
def get_dist_point_arc(A,B,C):
	#print "cross",A, B, C
	t0= time.time()
	D=np.cross(A,B)
	E=np.cross(D,C)
	U=np.cross(E,D)
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
