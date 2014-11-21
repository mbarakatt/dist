#!/usr/bin/env python
from  jaccard_plot import *
import numpy as np
import os
import sys
import geode


#allstuff=os.listdir("../results/2014.09.05.11h09m/")
#all_txt_files = [i for i in allstuff if i[-3:]=="txt"]
#file_path="../results/2014.11.06.13h11m/train_geo_out.txt"
file_path=sys.argv[1] + '/' + 'train_geo_out.txt'
save_folder='/'.join(file_path.split('/')[0:-1]) + '/'
fin=open(file_path,'r')

try:
	mindist = int(sys.argv[1])
	maxdist = int(sys.argv[2])
except:
	print "No min and/or max distance specified. Using 0 to 20000."
	mindist,maxdist = [0,20000]

def parsefile(fp):
	return map(lambda line: map(float,line.split(',')),  fp)

all_lines = np.array(parsefile(fin))
scores = all_lines[:,4]
max = np.max(scores)
min = np.min(scores)

desired_bounds="WORLD"
m=get_map(bounds=desired_bounds, showlines=True)

for line in all_lines[1:]:
	length = geode.get_dist([line[0],line[2]], [line[1],line[3]])[0,1]
	line = [ x*360/(2*np.pi) for x in line[:-1] ] + [ line[-1] ]
	#print "ASDFASDF", min,max,scores[0]
	if length >= mindist and length < maxdist :
		#print line[4],max,min,scores[0],scores[0]<line[4],1-(line[4]-scores[0])/( min-scores[0] )
		#print "LINE", line
		if scores[0] < line[4]:
			plot_great_line([ line[0], line[1] ],[line[2], line[3] ],m , color=( 0. , 0  , 1.  , (line[4]-scores[0])/(max-scores[0]) ) )
		else:
			plot_great_line([ line[0], line[1] ],[line[2], line[3] ],m , color=( 1, 0. , 0 , 1.  -(line[4]-min)/(scores[0]-min)))
		

try: 
	print save_folder + "train_position.txt"
	train = geode.line(*np.loadtxt(save_folder + "train_position.txt", delimiter='\t'),inputformat="degrees")	
except:
	print "No train Specified"
	train = None

if train != None:
	print "Drawing train."
	m = plot_great_line([train.lon1_d, train.lat1_d], [train.lon2_d, train.lat2_d], m)
try:
	os.mkdir("./pic/")
except:
	pass


plt.title('Max: %4f, No: %4f, Min: %4f' % (max,scores[0],min))
savefig(save_folder + str(maxdist) + ".jpg" ,m ,"From " + str(mindist) + "Km To " + str(maxdist)+ "km")
print "DONE"

