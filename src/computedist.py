from zipdist import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imread


ZC=Zip_Codes("../data/zipcode/zipcode.csv")
print "Done loading the zip codes"

print ZC.get_distance("07304","30006")


findiv_pos=open("/lb/project/gravel/data/SCCS/crossref/ors_124_115.txt",'r')
findiv_pos.readline()#removing the first line

indiv_posx=[]
indiv_posy=[]
for line in findiv_pos:
	sp = line.split()
	indiv_posx.append(float(sp[3]))
	indiv_posy.append(float(sp[4]))


img = imread("../data/usa-map.jpg")
plt.scatter(indiv_posx,indiv_posy,zorder=1)
plt.imshow(img,zorder=0, extent=[-95, -75,26, 40])
plt.show()
