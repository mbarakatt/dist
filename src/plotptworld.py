from pyproj import Proj
import matplotlib.pyplot as plt
import matplotlib
import random,math,scipy,sys
import numpy as np
from pylab import *


latlongfile="../data/language/latitudelongitude_americasnegative.txt"

f_latlongfile=open(latlongfile,'r')

lat=[]
long=[]
for line in f_latlongfile:
	sp=line.replace('\n','').split(',')
	lat.append(sp[0])
	long.append(sp[1])



p = Proj(proj='merc')#proj='utm',zone=10,ellps='WGS84')

x=[]
y=[]
for i in range(len(lat)):
	xi,yi=p(long[i],lat[i])
	x.append(xi)
	y.append(yi)

print long
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.plot(x,y,'.')
title('All Languages with x,y points')
savefig("all_language_xy.jpg")
