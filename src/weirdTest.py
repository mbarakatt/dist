import numpy as np
from scipy.stats import spearmanr
import time
class weirdTest:
	def __init__(self,distIndiv,distJaccard,maximumDistance=20000,nbBins=50):
		self.distIndiv=distIndiv
		self.tookTrainMask=tookTrainMask
		self.distJaccard=distJaccard
		self.numberTookTrainMask=np.sum(tookTrainMask)
		self.maximumDistance=maximumDistance
		self.nbBins=nbBins
		self.bins=np.linspace(0.00001, maximumDistances, nbBins+1) #create nbBins bins
		self.dc = (nbBins + (nbBins[1]-nbBins[0])/2.0)[0:-1]

	def __call__(self,line,bestDistances,tookTrainMask):
			
		countJaccardPerBin=np.histogram((distJaccard*tookTrainMask).flatten(),self.bins)[0]
		countPerBin = np.histogram((distIndiv*tookTrainMask).flatten(), self.bins)[0]
		polyCoefWO=np.polyfit(dc,self.countJaccardPerBin/self.countPerBin,4)
		
		countPerBin = np.histogram((bestDistances*tookTrainMask).flatten(), self.bins)[0]
		polyCoefW=np.polyfit(dc,self.countJaccardPerBin/self.countPerBin,4)
		return self.test(distIndiv*tookTrainMask,self.distJaccard*tookTrainMask,bestDistances*tookTrainMask)

	def test(self,distIndivMask, distJaccardMask, bestDistancesMask,polyCoefW,polyCoefWO):
		return  np.sum(np.absolute(np.polyval(polyCoefWO,distIndivMask)- distJaccardMask) - np.absolute(np.polyval(polyCoefWO,bestDistancesMask)-distJaccardMask))
