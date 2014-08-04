import numpy as np
import sys
import time
class mantelTest:
	def __init__(self,distIndiv,distJaccard,numberPermutated=4,distMask=20000):
		self.distIndiv=distIndiv
		self.distJaccard=distJaccard
		self.distMask=distMask
		self.permutatedJaccard =np.array([ np.random.permutation(np.random.permutation(distJaccard).T).T for _ in range(numberPermutated)])
		self.distIndivMask=self.getMask(self.distIndiv)
		self.numberItemInMask=np.sum(self.distIndivMask)
		#self.list_corr=spearmanMantelTest(dist_indiv,permutated_distJaccard)
		#self.corr=spearmanMantelTest(dist_indiv,np.array([distJaccard]))

	def getMask(self,M):
		return (M<=self.distMask)

	def __call__(self,line,bestDistances):
		return self.test(bestDistances, self.distJaccard)

	def meanStdMask(self,M, mask, item_in_mask):
		if len(M.shape) < 3:
				M=np.array([M])
		M_mask=M*mask
		mean=np.sum(M_mask,axis=(1,2))/item_in_mask
		return mean, np.sqrt(np.sum(np.power(((mean/np.sum(mask)).reshape((M_mask.shape[0],1,1)) - M_mask),2),axis=(1,2))/np.sum(mask))

	#M1 is a single matrix, M2 can be an list of matrix
	def test(self,M1,M2):
		if len(M2.shape) < 3:
			M2=np.array([M2])
		mask=self.distIndivMask
		item_in_mask=self.numberItemInMask
		M1_mean_mask, M1_std_mask= self.meanStdMask(M1,mask,item_in_mask)
		M2_mean_mask, M2_std_mask= self.meanStdMask(M2,mask,item_in_mask)
		#raise Exception(map(str,[M1_mean_mask,M1_std_mask,M2_mean_mask,M2_std_mask]))
		results=1.0/(item_in_mask - 1) * np.sum(((M1 - np.mean(M1)))/np.std(M1)*((M2 - np.mean(M2,axis=(1,2)).reshape((M2.shape[0],1,1)))/np.std(M2,axis=(1,2)).reshape((M2.shape[0],1,1))) ,axis=(1,2))
		return results[0]
