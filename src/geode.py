import numpy as np

import time


NEGLIGEABLE_RADIAN=0.005

EARTH_RADIUS=6371


def parse_file_geo_dist(filep):
        temp=[]
        c=0
        for line in filep:
                temp.append(map(float,line.split()))
                c+=1
        return np.array(temp)  


def toDegrees(list):
	return list*360/(2*np.pi)


def toRadians(list):
	return list*np.pi/180.0


def get_norm(x):	
	return  np.sqrt(np.sum(np.power(x,2),axis=1))


def get_angle(p1,p2):
	p1=np.array(p1)
	p2=np.array(p2)
#	print "get_angle", p1, p2, np.dot(p1,np.array(p2).transpose())
	return  np.abs(np.arccos(np.sum(p1*p2,axis=1)/(get_norm(np.array(p1))*get_norm(np.array(p2))))) #always positive


#Convert points (list of (x,y,z) ) from euclidean coordinates to spherical coordinates in gradients.
def get_spherical_coor(points):
	xs, ys, zs = map(np.array,zip(*points))
	lons=np.mod(np.arctan(ys/xs) + (xs<0)*np.pi, 2*np.pi )
	lats=np.arccos(zs/np.sqrt(xs*xs+ys*ys+zs*zs))
	result = np.array(zip(*[ lons*(lons <= np.pi) + (lons-2*np.pi)*(lons>np.pi),(np.pi/2-lats) ]))
	return result


#input in lons, lats format, returns in array of point format.
def get_euclidean_coor(lons,lats):
	"""
	lons: 1D array in radiants of longitudes
	lats: 1D array in radiants of latitudes
	"""
	adjusted_lats = np.copy(lats)
	if isinstance(lons,float):
		return EARTH_RADIUS * np.cos(lons) * np.cos(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.cos(adjusted_lats), EARTH_RADIUS*np.sin(adjusted_lats)
	else:
		return np.array([EARTH_RADIUS * np.cos(lons) * np.cos(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.cos(adjusted_lats), EARTH_RADIUS*np.sin(adjusted_lats)]).T



def get_dist_2pt(p1,ps2):
	ps2_lons, ps2_lats = map(np.array, zip(*ps2))
	#print 'ps2_', ps2_lons, ps2_lats
	return get_dist_one_many(p1[0], p1[1], ps2_lons, ps2_lats) #  + 1000000*arcIntersectProhibited(p1[0],p1[1],ps2_lons,ps2_lats)

def get_dist_one_many(lon1,lat1,lons,lats): # great circle distance.
	haversin = lambda theta : np.power(np.sin(theta/2.0),2) #(1-np.cos(theta))/2.0
	arg = haversin(lats-lat1) + np.cos(lats) * np.cos(lat1) * haversin( lons - lon1 )
	arg = 2 * EARTH_RADIUS * np.arcsin(np.sqrt(arg)) 
	return arg

#The arc is from A to B and C is the point. Here A, B are point on sphere(or vector starting at origin of sphere) in euclidean 3D space. C is a vector of the same nature
def get_dist_point_arc(A,B,C):
	D=np.cross(A,B)
	E=np.cross(D,C)
	U=np.cross(E,D)
	U=U*1./np.array([get_norm(U)]).transpose() * EARTH_RADIUS
	is_on_arc =  np.abs(get_angle([A],[B]) - ((get_angle([A],U)) + (get_angle([B],U)))) < NEGLIGEABLE_RADIAN
	t1=time.time()
	sphere_U, sphere_C, sphere_A, sphere_B = map(lambda x: np.array(get_spherical_coor(x)), [U,C,[A],[B]] )
	sphere_C_unzipped = np.array(zip(*sphere_C))
	sphere_U_unzipped = np.array(zip(*sphere_U))
	dist_A = get_dist_2pt(sphere_A[0], sphere_C)
	dist_B = get_dist_2pt(sphere_B[0], sphere_C)
	dist_U = get_dist_2pt_array(sphere_U, sphere_C)
	temp  = np.where(is_on_arc, dist_U, np.minimum(dist_A, dist_B))
	return temp

#Input are numpy arrAY OF PT IN RADIANS(LONS,LATS)
def get_dist_2pt_array(ps1,ps2):
	ps1_lons, ps1_lats = ps1.T
	ps2_lons, ps2_lats = ps2.T #map(np.array, zip(*ps2))
	return get_dist_one_many(ps1_lons, ps1_lats, ps2_lons, ps2_lats)

def git_dist_many_many(lons1,lats1,lons2,lats2):
 	arg = haversin(lats2-lats1) + np.cos(lats2) * np.cos(lats1) * haversin( lons2 - lons1 )
	arg = 2 * EARTH_RADIUS * np.arcsin(np.sqrt(arg)) 
	return arg

#Gets position in in format lons and lats (radiants) and returs the triplet xs,ys,zs
def get_dist(lons,lats): # great circle distance.
	matrix_lats, matrix_lons = map(lambda v : np.array([v for i in range(len(v))]) , [lats,lons])
	arg = np.power(np.sin((matrix_lats-np.array([lats]).transpose())/2.0), 2) + np.cos(np.array([lats]).transpose()) * np.cos(matrix_lats) * np.power(np.sin((matrix_lons-np.array([lons]).transpose())/2),2)
	arg = 2*np.arcsin(np.sqrt(arg)) #np.where( np.fabs(arg) < 1., arg, 0.999999)
	return EARTH_RADIUS * arg


def haversin(theta):
	return np.power(np.sin(np.array(theta)/2.0),2) #(1-np.cos(theta))/2.0

def get_dist_with_train(pos_euclidean,train):
	t,v = map(lambda x: list(get_euclidean_coor(x[0],x[1])), [train.p1, train.p2])
	a = get_dist_point_arc(t,v,pos_euclidean)
	b = np.array( [a] * len(pos_euclidean))
	return  b.T + b

#The option fullOutput gives the binary matrix of only the pairs of points that are using the train (ie that are closer to each other when they take the train
def computeDistances(train,distIndiv,pos_euclidean,fullOutput=False):
	added = get_dist_with_train(pos_euclidean,train)
	tookTrainMask = (added < distIndiv)
	dist_L = np.minimum(distIndiv , added)
	np.fill_diagonal(dist_L,0)
	if fullOutput:
		return dist_L, tookTrainMask
	else:
		return dist_L


class line:
	def __init__(self,lon1,lat1,lon2,lat2,inputformat='radiants'):
		if inputformat == 'degrees':
			self.lon1_d,self.lat1_d,self.lon2_d,self.lat2_d=lon1,lat1,lon2,lat2
			self.lon1,self.lat1,self.lon2,self.lat2=map(toRadians,[lon1,lat1,lon2,lat2])
		elif inputformat == 'radiants':
			self.lon1_d,self.lat1_d,self.lon2_d,self.lat2_d=map(toDegrees,[lon1,lat1,lon2,lat2])
			self.lon1,self.lat1,self.lon2,self.lat2=lon1,lat1,lon2,lat2
		self.p1,self.p2=[self.lon1,self.lat1],[self.lon2, self.lat2]

	
