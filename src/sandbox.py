import numpy as np
import matplotlib.pyplot as plt 

def point_density():
	dist_data=np.loadtxt('../SCCS/geodesicDist.txt',delimiter='\t')
	fig, ax = plt.subplots()
	dist_list=dist_data[np.triu_indices(dist_data.shape[0]),1]
	y_values,bins=np.histogram(dist_list)
	new_bins=bins[0:-1]+(bins[1]-bins[0])/2.0
	#print y_values, bins
	ax.set_title("Point Density")
	ax.set_xlabel("Distance (km)")
	ax.set_ylabel("nb of points")
	ax.set_xlim((0,1800))
	ax.bar(bins[0:-1],y_values,width=bins[1]-bins[0])
	plt.show()
	

point_density()

def some_function():
	f=plt.figure()
	ax=f.add_subplot(111)
	ax.set_ylim((-2,2))
	ax.set_xlim((-2,2))
	ax.plot(1,1,'.')
	ax.quiver([.5],[0.5],[1],[1], angles='xy', scale_units='xy' , scale=1) #Important note, all the params here must be there if you want the naturaly expected behaviour. Bigger scale makes arrow shorter
	plt.show()
