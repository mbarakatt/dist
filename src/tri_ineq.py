import numpy as np

dist_Matrix=np.loadtxt("../data/language/Ruhlen2014jaccard.txt",dtype=float)


for i in range(len(dist_Matrix)):
	for j in range(len(dist_Matrix)):
		result=  dist_Matrix[i,:] + dist_Matrix[:,j] - dist_Matrix[i,j] 
		if np.sum(result < 0) != 0 :
			print "Doesn't respect the triangle inequality "
			print i,j, result[(result <0).nonzero()]
