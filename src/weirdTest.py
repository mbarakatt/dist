import numpy as np
from pylab import *
from scipy.stats import spearmanr
import time
class weirdTest:
	def __init__(self, distIndiv,distJaccard, jaccardThreshold=0.22,maximumDistance=20000,nbBins=50):
		self.distIndiv = distIndiv.flatten()
		self.distJaccard = distJaccard
		self.triu_indices = np.triu_indices(self.distJaccard.shape[0],1)
		self.distJaccardArgSort = -np.argsort(self.distJaccard.flatten().reshape(self.distJaccard.shape))/(np.size(self.distJaccard)+1) + 0.5
		self.maximumDistance = maximumDistance
		self.jaccardThreshold=jaccardThreshold
		self.nbBins = nbBins
		self.bins = np.linspace(0.00001, self.maximumDistance, nbBins+1) #create nbBins bins
		self.dc = (self.bins + (self.bins[1]-self.bins[0])/2.0)[0:-1]
		self.jaccardBins = np.linspace(0.0,1.0,10+1)
		self.jaccardDc = (self.jaccardBins+ (self.jaccardBins[1]-self.jaccardBins[0])/2.0)[0:-1]
		self.y=self.get_ut_flatten(self.distJaccard)
		self.nb_y_below_t = np.sum( self.y <= jaccardThreshold )
		self.y_temp = np.argsort( np.argsort(( self.y > jaccardThreshold ) * 10000000 + self.y ))

	def __call__(self,line,bestDistances,tookTrainMask):
		retval=self.test(self.get_ut_flatten(bestDistances))
		#if retval !=0.5:
		#	print "Number to take train:", np.sum(tookTrainMask),tookTrainMask
	
		return retval
		#self.numberTookTrainMask=np.sum(tookTrainMask)
		#self.distIndivTookTrain=self.distIndiv*tookTrainMask
		#self.distJaccardTookTrain=self.distJaccard*tookTrainMask
		#self.bestDistances = bestDistances
		#return np.sum(np.nan_to_num((self.distIndiv-self.bestDistances)/self.distIndiv)*(self.distJaccardArgSort))

	def get_ut_flatten(self,M):
		#flat=np.triu(M,k=1).flatten()
		#return flat[np.nonzero(flat)]
		return M[self.triu_indices]

	def spear(self,xss):
	#return 1-6*np.sum(np.power((np.argsort(np.argsort(xss[xss!=0]))+1)-(np.argsort(np.argsort(yss[xss!=0]))+1),2))/float((np.power(np.sum(xss!=0),3) - np.sum(xss!=0)))
		#txss=xss[0:5]
		#tyRanks=self.yRanks[0:5]
		yRanks=np.amin([self.y_temp, np.argsort(np.argsort( (self.y <= self.jaccardThreshold )*10000000 + xss )) + self.nb_y_below_t],axis=0)
		#print "avg", np.mean(xss),np.mean(yRanks)
		#print "tiny", txss, tyRanks
		retval= 1-6*np.sum(np.power((np.argsort(np.argsort(xss))) - (yRanks),2))/float((np.power(xss.size,3) - xss.size))
		#if retval != 0.5:
		#	print "Rank", np.argsort(np.argsort(xss)), yRanks , "Dist", xss,self.y 
		
		return retval
		#return np.sum( map(lambda (xs,ys) :np.sum(xs!=0)*np.nan_to_num(1 - 6*np.sum(np.power((np.argsort(xs[xs!=0])+1)-(np.argsort(ys[xs!=0])+1),2))/float(np.power(np.sum(xs!=0),3)-np.sum(xs!=0))) , np.array([xss,yss]).swapaxes(0,1)))/float(np.sum(xss[xss!=0]))

	def test(self,  bestDistancesMask):
		#return (self.spear(bestDistancesMask.flatten(),distJaccardMask.flatten())- self.spear(distIndivMask.flatten(), distJaccardMask.flatten()))
		#tinyBestDistancesMask=bestDistancesMask[0:2,0:2]
		#tinyDistJaccardMask=distJacardMask[0:2,0:2]
		#print "TinyDists:", tinyBestDistancesMask, tinyDistJaccardMask
		a = self.spear(bestDistancesMask)
		#b= self.spear(self.distIndiv) 
		#print "ab",a,b
		return a
		#return  np.sum(np.absolute(np.polyval(polyCoefWO,distIndivMask)- distJaccardMask) - np.absolute(np.polyval(polyCoefWO,bestDistancesMask)-distJaccardMask))/self.numberTookTrainMask



#OldStuff
		#--------------To get the times closer plot		
		#countPerJaccardBin=np.histogram(self.distJaccard.flatten() , self.jaccardBins)[0]
		#avgNormalDistPerJaccardBin=np.histogram(self.distJaccard.flatten(), self.jaccardBins,weights=self.distIndiv.flatten())[0] /countPerJaccardBin
		#avgBestDistPerJaccardBin=np.histogram(self.distJaccard.flatten(), self.jaccardBins,weights=self.bestDistances.flatten())[0] /countPerJaccardBin	
		#print "test1", self.jaccardDc
		#print "test2", avgNormalDistPerJaccardBin
		#print "test3", avgBestDistPerJaccardBin
		#plot(self.jaccardDc,avgNormalDistPerJaccardBin/avgBestDistPerJaccardBin)
		#ylabel("Times closer")
		#xlabel("Jaccard Bins")
		#title("Effect on avg distances in Jaccard bins of putting a particular Train")
		#show()
		#return 0
		#--------------
		#r=getRank(self.bestDistances)
		#s1=getRank(self.jaccard[self.jaccard<=self.jaccardThreshold],ignoreZeros=False)
		#s= np.append(b1 ,getRank(self.bestDistances[self.jaccard> self.jaccardThreshold] + len(b1))
		#s[s>len(b1)]=len(b1)

		#countJaccardPerBin=np.histogram((self.distIndiv).flatten(), self.bins,weights=(self.distJaccard).flatten())[0]
		#countPerBin = np.histogram((self.distIndiv).flatten(), self.bins)[0]
		#return self.test(self.distIndiv,self.distJaccard,bestDistances) 
