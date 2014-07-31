
class spearmanMantelTest:
	def __init__(self,distIndiv,distJaccard):
		self.distIndiv=distIndiv
		self.distJaccard=distJaccard
	def __call__(self,line,bestDistances):
		return test(bestDistances,self.distJaccard)	
	def spear(xs,ys):
		return [scipy.stats.spearmanr(xs,y)[0] for y in ys]
	def test(M1,M2):
		if len(M2.shape) < 3:
			M2=np.array([M2])
		M2_flatten=np.array(map(lambda x : x.flatten(), M2))
		return spear(M1_flatten,M2_flatten)
