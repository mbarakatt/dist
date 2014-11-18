import numpy as np
from jaccard_plot import * 


def scatterxy(ax,xss,yss,ylim=(0,0),xlim=(0,0)):
	"""This methods returns an axe object with a scatter plot of the data in the upper triangle of xss and yss."""
	triu=np.triu_indices(len(xss),1)
	if ylim != (0,0):
		ax.set_ylim(ylim)
	if xlim != (0,0):
		ax.set_xlim(xlim)
	ax.scatter(xss[triu],yss[triu])
	return ax	

def show_scatterxy(geodesic_dist_matrix,IBD_matrix,out_folder):
	"""This method makes use of scatterxy define in the same module to show the actual picture."""
	f = plt.figure()
	ax = f.add_subplot(111)
	ax = scatterxy(ax, geodesic_dist_matrix, IBD_matrix)
	plt.savefig(out_folder + 'scatterxy.png')

def draw_point_IBD(positions,jaccard, output_filename):
	m=plot_map_relatedness(get_map(),jaccard,positions)
	savefig(output_filename,m)
