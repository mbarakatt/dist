import numpy as np
from scipy.stats import rankdata
from pylab import *
from weirdTest import *
import time 

#pts=np.random.random(2*50).reshape(2,50)
pts = np.array([np.linspace(0,1,2500*2500).reshape(2500,2500)]*2) #np.array([[0.1,.2,.3,.4,0.5,0.6,0],[0.1,.2,.3,.4,.5,.6,0.25]])
print "Done generating data"
LARGE_NUMBER = iinfo('i').max
mytest = weirdTest(pts[1,:],jaccardThreshold=0.5, negateyaxis=True)

def rank(xs):
	return rankdata(xs)

def passThreshold( ys, t ):
	return ( ys <= t )

def computeCorr(xs,ys,t):
	PassT = passThreshold(ys,t)
	rankxs=rank(xs)
	numPassT = np.sum(PassT)
	ys_upper = rank(xs + PassT*LARGE_NUMBER)  + numPassT 
	new_ys = np.amin(np.array([rank(ys) + (~PassT) * LARGE_NUMBER, ys_upper + PassT*LARGE_NUMBER]),axis=0 )
	rankxs, rankys=rank(xs)	, new_ys
	for i in range(len(xs)):
		annotate(str(rankxs[i])+','+ str(rankys[i]), xy=(xs[i],ys[i]) )
	return 1 - 6*np.sum( np.power( rankxs - rankys ,2 ) )/float(np.power(xs.size,3))


for i in np.linspace(0.11,0.71,7):
	figure()
	pts[0,-1]=i
	t0=time.time()
	corr=mytest(2,pts[0,:],3)#computeCorr(pts[0,:],pts[1,:],0.5)
	print "Time:", time.time() - t0,corr
	continue
	plot(pts[0,:], pts[1,:], '.' )
	plot([0,1],[0.5,0.5])
	title(str(corr))
	ylim(0,1)
	xlim(0,1)
	savefig("pic/temp/corr%f.jpg" % i)

