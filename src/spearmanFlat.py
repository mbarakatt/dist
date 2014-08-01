import numpy as np
from scipy.stats import spearmanr
import time 
class spearmanFlat:
	def __init__(self,distIndiv,distJaccard):
		self.distIndiv=distIndiv
		self.distJaccard=distJaccard
	def __call__(self,line,bestDistances):
		return	self.test(bestDistances,self.distJaccard)	
	def spear(self,xs,ys):
		return [spearmanr(xs,y)[0] for y in ys]
		#xs_argsort=np.argsort(xs)
		#ys_argsort=map(lambda y : np.argsort(y),ys)
		#return  1 - (6* np.sum(np.power(xs_argsort-ys_argsort, 2), axis=1))/(np.power(len(xs),3)-float(len(xs)))
	def test(self,M1,M2):
		M1_flatten=M1.flatten()
		if len(M2.shape) < 3:
			M2=np.array([M2])
		M2_flatten=np.array(map(lambda x : x.flatten(), M2))
		return self.spear(M1_flatten,M2_flatten)[0]
