import numpy as np
from scipy.stats import spearmanr
import time
class spearmanMantelTest:
	def __init__(self,distIndiv,distJaccard,numberPermutated=4 ):
		self.distIndiv=distIndiv
		self.distJaccard=distJaccard
		self.permutatedJaccard =np.array([ np.random.permutation(np.random.permutation(self.distJaccard).T).T for _ in range(numberPermutated)])

	def __call__(self,line,bestDistances):
		listCorr=self.test(bestDistances, self.permutatedJaccard)	
		corr=self.test(bestDistances,self.distJaccard)
		return 	(corr - np.mean(listCorr))

	def spear(self,xs,ys):
		return [spearmanr(xs,y)[0] for y in ys]


	def test(self,M1,M2):
	#Performs a weird mantel test on the rank matrices of M1 and M2 (where the rank is taken by row)
	#Note that the result is NOT normalized
		if len(M2.shape)<3:
			M2=np.array([M2])
		value = 1-6*np.sum(np.power(np.argsort(M1)-np.argsort(M2),2),axis=(1,2))/float(np.power(len(M1),3)-(len(M1)))	
		if len(value)==1:
			return value[0]
		else:
			return value
