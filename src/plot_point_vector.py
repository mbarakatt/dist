#!/usr/bin/env python
#import jaccard_plot
import sys
import numpy as np
from matplotlib import use
#use('Agg')
import matplotlib.pyplot as plt

QUIVER_SCALE=80000
data = np.loadtxt("temp/all.txt", delimiter='\t' , usecols = (0,1,2,3,4))
data2 = np.genfromtxt("temp/all.txt", delimiter='\t' , usecols = range(5,data.shape[0]+5),dtype=None)
xss = np.loadtxt("../SCCS/geodesicDist.txt", delimiter='\t')
yss = np.loadtxt("../SCCS/cM_Matrix_18.txt", delimiter='\t')

def linearize_distance_matrix(M):
	return M[np.triu_indices(M.shape[0],1)]

xs = linearize_distance_matrix(xss)
ys = linearize_distance_matrix(yss)
sorted_data = data[data[:,0].argsort()].T
N, X ,Y ,U ,V = sorted_data
#print "b"
#indices_matrix= np.array([[str(i)+ ',' + str(j)  for j in range(len(ys))] for i in range(len(ys))])  
#print "e"
#index_above=(yss[i] > 40 )* indices_matrix

def scatterxy(ax,xs,ys,xlim,ylim):
	ax.plot(xs, ys, '.')
	ax.set_xlim(*xlim)
	ax.set_ylim(*ylim)
	return ax

def simple_quiver(ax,X,Y,U,V):
	ax.plot(X,Y,'.')
	ax.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=QUIVER_SCALE,width=1/700.0)#smaller scale makes arrow longer
	return ax

def produce_simple_quiver():
	print "producing simple quiver"
	f=plt.figure(figsize=(10,10))
	ax=f.add_subplot(111)
	ax=simple_quiver(ax,X,Y,U,V)
	#plt.show()
	return f,ax
f,ax=produce_simple_quiver()

def ask_particular_point():
	global data2
	data2=data2[data[:,0].argsort()]
	
	above=True
	if above:#Above is good for theshold.
		indices_above = np.array(np.nonzero(yss > 40)).T
		outname_prefix='a'
	else:
		indices_above = np.array(np.nonzero(yss < 40)).T
		outname_prefix='b'
	print "running for arg" , sys.argv[1]	
	run_for = tuple(indices_above[int(sys.argv[1])]) #np.array(np.max(yss,axis=0)> 40 )
	
	print run_for
	if (X[run_for[1]],Y[run_for[1]]) == (X[run_for[0]],Y[run_for[0]]):
		exit(1)	
	f=plt.figure(figsize=(12, 6))
	ax1 = f.add_subplot(121)
	ax2 = f.add_subplot(122,autoscale_on=False)
	ax1.plot(X,Y,'.')
	f,b = map(float,data2[run_for].split(','))
	
	vect=np.array([X[run_for[1]]-X[run_for[0]], Y[run_for[1]]-Y[run_for[0]]])
	unit_vect= vect / np.linalg.norm(vect)
	
	ax1.quiver(X[run_for[0]],Y[run_for[0]], unit_vect[0]*(f-b),+unit_vect[1]* (f-b),angles='xy',scale_units='xy',scale=10000,width=1/400.0)
	ax1.plot(X[run_for[1]],Y[run_for[1]],'.',color='red')
	ax1.plot(X[run_for[0]],Y[run_for[0]],'.',color='green')
	ax=scatterxy(ax2,xs,ys,(0,1800),(0,150))
	ax2.plot([xss[run_for]-200,xss[run_for]+200], [yss[run_for],yss[run_for]],color='red')
	ax1.set_title(map(str,[f,b]))
	ax2.set_title(map(str,[xss[run_for],yss[run_for]]))
	#plt.draw()
	print "Saving figure to ", "arrows/%s%d" % (outname_prefix,int(sys.argv[1]))
	plt.savefig("arrows/%s%d" % (outname_prefix,int(sys.argv[1])))
#ask_particular_point()


#print clinic_pos_list,unique_clinic_pos_list
def show_clinic_average(ax):
	#Clinic Arrows:
	clinic_pos_list=np.array([X,Y]).T
	unique_clinic_pos_list=np.unique(map(tuple,clinic_pos_list))
	for clinic in unique_clinic_pos_list:
		print 'clinic',clinic
		same_clinic=np.array([(np.absolute(item[0] - clinic[0]) < 0.0001 ) and (np.absolute(item[1] - clinic[1]) < 0.0001 ) for item in clinic_pos_list])
		print np.sum(same_clinic)
		print "UV", len(U[same_clinic]),len(V[same_clinic])
		vc= np.array([np.sum(U[same_clinic]) , np.sum(V[same_clinic])])*1./np.sum(same_clinic)
		print "UV, details", U[same_clinic],"avg", vc[0], V[same_clinic], U,'avg', vc[1]
		print "VC" , vc
		ax.quiver(clinic[0],clinic[1],vc[0],vc[1],angles='xy',scale_units='xy',scale=80000/10.0,width=1/500.0,color='red')
	return ax


ax=show_clinic_average(ax)
plt.show()

