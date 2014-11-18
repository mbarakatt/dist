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
import time

NB_BINS=100
EARTH_RADIUS=6371
NEGLIGEABLE_RADIAN=0.005
jaccard_file_path="../data/language/Ruhlen2014jaccard.txt"
parse_jaccard(open(jaccard_file_path,'r'))

relatedness = np.array(parse_jaccard(open(jaccard_file_path,'r')))
#relatedness=np.array([[0,0.5,0.2],[0.5,0,0.1],[0.2,0.1,0]])

longlatfile = "../data/language/longitudelatitude_americasnegative.txt"
f_longlatfile = open(longlatfile, 'r')
lines_longlatfile = f_longlatfile.read().split('\n')
#pos = [np.array([0,10,20]),np.array([0,10,20])]
pos=map(lambda x : np.array(map(float,x)), zip(*[ line.split('\t') for line in lines_longlatfile ]))

lons = pos[0]
lats = pos[1]
lats, lons = map(lambda x: x*2*np.pi/360,[lats,lons]) #correct that

def get_dist(lons,lats): # great circle distance.
	matrix_lats, matrix_lons = map(lambda v : np.array([v for i in range(len(v))]) , [lats,lons])
	arg = np.power(np.sin((matrix_lats-np.array([lats]).transpose())/2.0), 2) + np.cos(np.array([lats]).transpose()) * np.cos(matrix_lats) * np.power(np.sin((matrix_lons-np.array([lons]).transpose())/2),2)
	arg = 2*np.arcsin(np.sqrt(arg))
	return EARTH_RADIUS * arg



def get_euclidean_coor(lons,lats):
	adjusted_lats = np.pi/2 - np.copy(lats)
	return EARTH_RADIUS * np.cos(lons) * np.sin(adjusted_lats), EARTH_RADIUS*np.sin(lons)*np.sin(adjusted_lats), EARTH_RADIUS*np.cos(adjusted_lats)



dist_indiv    = get_dist(lons,lats)
pos_euclidean = np.array(get_euclidean_coor(lons,lats))

maxdist = 20000 #np.max(b)
bins = np.linspace(0,maxdist, NB_BINS+1) #create NB_BINS bins
dc = (bins + (bins[1]-bins[0])/2.0)[0:-1]
best_min = "nan"
best_sigmoid_params=[]



a = np.ndarray.flatten(relatedness)
b = np.ndarray.flatten(dist_indiv)
relatedness = np.histogram(b,bins,weights=a)[0]/2.0
count = np.histogram(b,bins)[0]

plot(dc, relatedness/count)
plt.gca().invert_yaxis()
#ylim(0,1)
savefig("geo_dist.jpg")






