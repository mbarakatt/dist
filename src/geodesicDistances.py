import numpy as np
from geode import *
from time import *

NEGLIGEABLE_RADIAN=0.005
longlatfile = "../data/language/longitudelatitude_americasnegative.txt"
outFile="geographicDistancesLanguages.txt"
fOutFile=open(outFile,'w')

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




print "Test"

dist=np.zeros(len(lons)*len(lons)).reshape(len(lons),len(lons))
posZipped=zip(lons,lats)
print posZipped[0][1]
for pos1 in range(len(posZipped)):
	t0=time()
	for pos2 in range(len(posZipped)):
		if  pos1 < pos2 :
			dist[pos1,pos2] = computeDistancesLanguages(posZipped[pos1][0],posZipped[pos1][1],posZipped[pos2][0],posZipped[pos2][1])
			print pos2
	print "pos1 ", pos1, "Time:" , time()-t0
print dist 
	
for i in len(lons):
	for j in len(lons):
		fOutFile.write(str(dist[i,j]) )
		if j == len(lons) -1:
			fOutFile.write("\t")
		else:
			fOutFile.write("\n")
