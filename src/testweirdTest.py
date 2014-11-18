#The goal of this code is to test if the weirdTest class is working properly
from weirdTest import weirdTest
import numpy as np
import matplotlib.pyplot as plt

#dist = np.arange(144).reshape(12,12)
#related = np.fliplr([np.linspace(0,144,144)])[0].reshape(12,12)
temp_dist=np.loadtxt('../SCCS/geodesicDist.txt',delimiter='\t')
temp_related=np.loadtxt('../SCCS/cM_Matrix_18.txt',delimiter='\t')
dist = temp_dist[0::20,0::20]
related = temp_related[0::20,0::20]
corr=[]
for i in range(7):
	#dist[0,2]=i
	#related[0,2]=80
	myweirdTest = weirdTest(dist,related,threshold=40,negateyaxis=True)
	f = plt.figure()
	ax=f.add_subplot(111)
	print "before:", dist[(related>79) & (related<82)]
	dist[(related>79) & (related<82)] -= 10*i
	print "after:", dist[(related>79) & (related<82)]
	ax.plot(dist[myweirdTest.triu_indices],related[myweirdTest.triu_indices],'.')
	print len(list(related[myweirdTest.triu_indices]))
	corr.append(myweirdTest([0],dist,[0]))
	#ax.set_title(corr[-1])
	print corr[-1]
	#plt.show()

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(range(7),corr)
plt.show()
