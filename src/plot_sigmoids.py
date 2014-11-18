import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from  pylab import *
import random,math,scipy,sys
import numpy as np
np.set_printoptions(threshold=np.nan)
from scipy.stats.stats import pearsonr
from itertools import product

def sigmoidf(d,params):
	return params[0]/(1.0 + math.exp((params[1])*max(d, 0.01)+params[2])) #smaller number means more IBD

#sigmoid_params_search=[np.linspace(0.5,2,6),np.linspace(50,70,6),np.linspace(0, 0.5, 6)]
sigmoid_params_search=[[1],np.linspace(50,70,6),[0.1]]


xs=np.linspace(0.0125,0.25,10)
fig=figure()
ax = fig.add_subplot(1,1,1)
for cur_params in product(*sigmoid_params_search):
	ys=[ sigmoidf(di,cur_params) for di in xs]
	ax.plot(xs,ys)

savefig("sigmoids_onlyp2.jpg")
