from geode import *
import sys
import numpy as np
fout=open("./test_longlat.txt",'w')

try:
	nbpts=int(sys.argv[1])
except:
	print "Usage python randomlonglat.py nbpts"
	exit(1)

rnd_vec = np.random.random(nbpts*3).reshape(nbpts,3)-0.5
 
normalized_vec=np.array([x/np.linalg.norm(x) for x in rnd_vec])*EARTH_RADIUS
#print normalized_vec
spherical_coors=toDegrees(get_spherical_coor(normalized_vec))
#print spherical_coors

np.savetxt(fout,spherical_coors,delimiter="\t")
