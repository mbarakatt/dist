import numpy as np
from pylab import *
from scipy.stats import spearmanr
import time
class weirdTest:
	def __init__(self, distIndiv,distJaccard, jaccardThreshold=1.0,maximumDistance=20000,nbBins=50):
		self.distIndiv = distIndiv
		self.distJaccard = distJaccard
		self.distJaccardArgSort = -np.argsort(self.distJaccard.flatten().reshape(self.distJaccard.shape))/(np.size(self.distJaccard)+1) + 0.5
		self.maximumDistance = maximumDistance
		#self.jaccardThreshold
		self.nbBins = nbBins
		self.bins = np.linspace(0.00001, self.maximumDistance, nbBins+1) #create nbBins bins
		self.dc = (self.bins + (self.bins[1]-self.bins[0])/2.0)[0:-1]
		self.jaccardBins = np.linspace(0.0,1.0,10+1)
		self.jaccardDc = (self.jaccardBins+ (self.jaccardBins[1]-self.jaccardBins[0])/2.0)[0:-1]

#	def getRank(M,filterMatrix=[],ignoreZeros=True):#Ignores 0 entries by default
#		if filterMatrix == []:
#			filterMatrix=M
#		return np.argsort((np.argsort(xss[filterMatrix!=0 ) | not ignoreZeros]))+1


	def __call__(self,line,bestDistances,tookTrainMask):
		self.numberTookTrainMask=np.sum(tookTrainMask)
		self.distIndivTookTrain=self.distIndiv*tookTrainMask
		self.distJaccardTookTrain=self.distJaccard*tookTrainMask
		self.bestDistances = bestDistances

		return np.sum(np.nan_to_num((self.distIndiv-self.bestDistances)/self.distIndiv)*(self.distJaccardArgSort))
		
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

		#countJaccardPerBin=np.histogram((self.distIndiv).flatten(),self.bins,weights=(self.distJaccard).flatten())[0]
		#countPerBin = np.histogram((self.distIndiv).flatten(), self.bins)[0]
		#return self.test(self.distIndiv,self.distJaccard,bestDistances) 

	def spear(self,xss,yss):
		return 1-6*np.sum(np.power((np.argsort(np.argsort(xss[xss!=0]))+1)-(np.argsort(np.argsort(yss[xss!=0]))+1),2))/float((np.power(np.sum(xss!=0),3) - np.sum(xss!=0)))
		#return np.sum( map(lambda (xs,ys) :np.sum(xs!=0)*np.nan_to_num(1 - 6*np.sum(np.power((np.argsort(xs[xs!=0])+1)-(np.argsort(ys[xs!=0])+1),2))/float(np.power(np.sum(xs!=0),3)-np.sum(xs!=0))) , np.array([xss,yss]).swapaxes(0,1)))/float(np.sum(xss[xss!=0]))

	def test(self,distIndivMask, distJaccardMask, bestDistancesMask):
		#return (self.spear(bestDistancesMask.flatten(),distJaccardMask.flatten())- self.spear(distIndivMask.flatten(), distJaccardMask.flatten()))
		return self.spear(bestDistancesMask,distJaccardMask)- self.spear(distIndivMask,distJaccardMask)
		#return  np.sum(np.absolute(np.polyval(polyCoefWO,distIndivMask)- distJaccardMask) - np.absolute(np.polyval(polyCoefWO,bestDistancesMask)-distJaccardMask))/self.numberTookTrainMask
