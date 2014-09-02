import numpy as np
import time
NEGLIGEABLE_RADIAN=0.005
EARTH_RADIUS=6371


def get_norm(x):	
	return  np.sqrt(np.sum(np.power(x,2),axis=1))

def get_angle(p1,p2):
	p1=np.array(p1)
	p2=np.array(p2)
#	print "get_angle", p1, p2, np.dot(p1,np.array(p2).transpose())
	return  np.abs(np.arccos(np.sum(p1*p2,axis=1)/(get_norm(np.array(p1))*get_norm(np.array(p2))))) #always positive

def get_spherical_coor(points):
	xs, ys, zs = map(np.array,zip(*points))
	lons=np.mod(np.arctan(ys/xs) + (xs<0)*np.pi, 2*np.pi )
	lats=np.arccos(zs/np.sqrt(xs*xs+ys*ys+zs*zs))
	result = np.array(zip(*[ lons*(lons <= np.pi) + (lons-2*np.pi)*(lons>np.pi),(np.pi/2-lats) ]))
	return result
#input in lons, lats format, returns in array of point format.
def get_euclidean_coor(lons,lats):
	adjusted_lats = np.copy(lats)
	if isinstance(lons,float):
		return EARTH_RADIUS * np.cos(lons) * np.cos(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.cos(adjusted_lats), EARTH_RADIUS*np.sin(adjusted_lats)
	else:
		return np.array([EARTH_RADIUS * np.cos(lons) * np.cos(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.cos(adjusted_lats), EARTH_RADIUS*np.sin(adjusted_lats)]).T

def get_dist(lons,lats): # great circle distance.
	matrix_lats, matrix_lons = map(lambda v : np.array([v for i in range(len(v))]) , [lats,lons])
	arg = np.power(np.sin((matrix_lats-np.array([lats]).transpose())/2.0), 2) + np.cos(np.array([lats]).transpose()) * np.cos(matrix_lats) * np.power(np.sin((matrix_lons-np.array([lons]).transpose())/2),2)
	arg = 2*np.arcsin(np.sqrt(arg)) #np.where( np.fabs(arg) < 1., arg, 0.999999)
	return EARTH_RADIUS * arg
def haversin(theta):
	np.power(np.sin(theta/2.0),2) #(1-np.cos(theta))/2.0

