import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

xss=np.loadtxt("../SCCS/geodesicDist.txt",delimiter='\t')
yss=np.loadtxt("../SCCS/cM_Matrix_18.txt",delimiter='\t')

def linearize_distance_matrix(M):
	return M[np.triu_indices(M.shape[0],1)]

xs = linearize_distance_matrix(xss)
ys = linearize_distance_matrix(yss)



f=plt.figure(figsize=(10,10))

ax= f.add_subplot(111)

y_values=[]
x_values=[]
over_y_values=[]
under_y_values=[]

def get_val(s,e):
	index_in_bin=np.nonzero((xs < e) & (xs >= s))
	print index_in_bin[0].size
	val=np.sum(ys[index_in_bin])/index_in_bin[0].size
	var_over_sqrt_n= np.std(ys[index_in_bin])/np.sqrt(index_in_bin[0].size)
	print s,e,val
	y_values.append(val)
	x_values.append((e-s)/2.0+s)
	over_y_values.append(val + 2*var_over_sqrt_n)
	under_y_values.append(val - 2*var_over_sqrt_n)
	

def gen_intervals(start,end,step):
	a=[]
	for i in range(start,end,step):
		a.append([i,i+step])
	return a

for s,e in [[0,0.1]] + gen_intervals(1,1701,100):
	get_val(s,e)

ax.set_yscale('log')
def plot_one_curve(ax,label):
	ax.plot(x_values,over_y_values,'--',color="grey",linewidth=3)
	ax.plot(x_values,under_y_values,'--',color="grey",linewidth=3)
	ax.plot(x_values, y_values,linewidth=4,label=label)
	print "X then Y values"
	print x_values
	print y_values
	ax.set_ylabel("Average pairwise IBD (cM)",fontsize=30)
	ax.set_ylim((0,2.35*max(y_values)))
	ax.set_xlabel("Distance (km)",fontsize=30)
	return ax

ax=plot_one_curve(ax,r'$\geq 18$cM')

y_values=[]
x_values=[]
over_y_values=[]
under_y_values=[]

yss=np.loadtxt("../SCCS/cM_Matrix.txt",delimiter='\t')
ys = linearize_distance_matrix(yss)

for s,e in [[0,0.1]] + gen_intervals(1,1701,100):
	get_val(s,e)

ax=plot_one_curve(ax,r'$\geq 3$cM')

#ax.set_yticklabels(map(str,np.log(range(0,10))),fontsize=17)
ax.tick_params(axis='y', labelsize=25)
ax.tick_params(axis='x', labelsize=25)
plt.yticks([0.1,1],["0.1","1.0"])
plt.xticks(np.arange(0,1501,500))
handles, labels = ax.get_legend_handles_labels()
handles.append(mpl.lines.Line2D([0],[0],linestyle='--',color='grey', linewidth=3))
labels.append(r'$2 \times$ s.e.m.')
ax.legend(handles, labels,fontsize=30)

plt.savefig("temp/avg18.jpg")

