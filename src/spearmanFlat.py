
plass spearmanFlat:
	def __init__(self,distIndiv,distJaccard):
		self.distIndiv=distIndiv
		self.distJaccard=distJaccard
	def __call__(self,line):
		bestDistances=computeDistances(line,self.distIndiv,self.distJaccard)
		return	test(bestDistances,self.distJaccard)	
	def spear(xs,ys):
		return [scipy.stats.spearmanr(xs,y)[0] for y in ys]
		#xs_argsort=np.argsort(xs)
		#ys_argsort=map(lambda y : np.argsort(y),ys)
		#return  1 - (6* np.sum(np.power(xs_argsort-ys_argsort, 2), axis=1))/(np.power(len(xs),3)-float(len(xs)))
	def test(M1,M2):
		t0=time.time()
		print "Starting spearmanFlat", 
		M1_flatten=M1.flatten()
		if len(M2.shape) < 3:
			M2=np.array([M2])
		M2_flatten=np.array(map(lambda x : x.flatten(), M2))
		print "spearmanFlat test time: ", time.time()-t0
		return spear(M1_flatten,M2_flatten)
