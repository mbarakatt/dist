
class mantelTest:
	def __init__(self,distIndiv,distJaccard,numberPermutated=4,distMask=7000):
		self.distIndiv=distIndiv
		self.distJaccard=distJaccard
		self.distMask=distMask
		self.permutatedJaccard =np.array([ np.random.permutation(np.random.permutation(relatedness).T).T for _ in range(numberPermutated)])
		self.distIndivMask=getMask(self.distIndiv)
		self.numberItemInMask=np.sum(self.distIndivMask)
	def get_MASK(M):
		return (M<=self.distMask)
	def __call__(self,line):
		t0=time.time()
		print "Starting mantelTest...",
		bestDistances=computeDistances(line,self.distIndiv,self.distJaccard)
		test(bestDistances, self.distJaccard)
		print "mantelTest test time:", time.time() - t0
	#M1 is a single matrix, M2 can be an array of matrix
	def test(M1,M2):
		if len(M2.shape) < 3:
			M2=np.array([M2])
		mask=self.distIndivMask
		item_in_mask=self.numberItemInMask
		M1_mean_mask, M1_std_mask= mean_std_MASK(M1,mask,item_in_mask)
		M2_mean_mask, M2_std_mask= mean_std_MASK(M2,mask,item_in_mask)
		results=1.0/(item_in_mask - 1) * np.sum(((M1 - M1_mean_mask)*mask)/M1_std_mask*((M2 - M2_mean_mask.reshape((M2.shape[0],1,1)))*mask)/M2_std_mask.reshape((M2.shape[0],1,1)) ,axis=(1,2))
		return -results
