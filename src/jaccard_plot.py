from __future__ import print_function

import matplotlib
#matplotlib.use('Agg')

from mpl_toolkits.basemap import Basemap
from itertools import imap
import matplotlib.pyplot as plt
import numpy as np
from math import atan2, cos, sin, sqrt, asin

import sys

corner_names			= ['llcrnrlon', 'llcrnrlat', 'urcrnrlon', 'urcrnrlat']
make_bbox			   = lambda llcrnr, urcrnr: dict(zip(corner_names, llcrnr + urcrnr))
BOUNDS_				 = [("WORLD",		 ((-180, -90), (180, 90))),
						   ("MAINLAND_USA",  ((-130,  25), (-60, 50))),
						   ("AMERICAS",	  ((-170, -60), (-30, 75))),
						   ("SCCS",	  ((-100, 25), (-70, 42))),
						   ("NORTH_AMERICA", ((-170,  30), (-30, 75))),
						   ("SOUTH_AMERICA", ((-90,  -60), (-30, 15))),
						   ("AFRICA",		((-30,  -45), (60,  45))),
						   ("EUROPE",		((-15,   30), (45,  50)))]
BOUNDS				  = dict(map(lambda (k, v): (k, make_bbox(*v)), BOUNDS_))

latlon			  = True
threshold		   = 0.2
desired_bounds	  = "SCCS"

def is_in_bounds((lon, lat), bounds):
	return (lon > bounds['llcrnrlon'] and lon < bounds['urcrnrlon'] and
			lat > bounds['llcrnrlat'] and lat < bounds['urcrnrlat'])

def show_usage():
	p = lambda *m: print(*m, file=sys.stderr)
	p("jaccard_plot.py -- display a relatedness graph of many languages.")
	p("usage: jaccard_plot.py -p <positions file> -j <jaccard file> [-o <output file>] [-l] [-t <threshold>] [-b <bounds>]")
	p("  The positions file is expected to be in the format 'latitude,longitude'. If it is")
	p("stored in the reverse format, then the -l switch can be used to flip the interpretation.")
	p("Only when two language are distanced less than `threshold` will their distance-line be drawn.")
	p("The default threshold is 0.22")
	p("More line can be added on the graph with the -m switch, beware the input format is weird. By default the first 20 trains are drawn.")
	p("The desired bounds to examine can be indicated with the -b switch. The default bounds contain the whole world.")
	p("Valid bounds consist of:")
	map(lambda (name, ((lon1, lat1), (lon2, lat2))): p("\t", name, " from (", lon1, ",", lat1, ") to (", lon2, ",", lat2, ")"), BOUNDS_)

def parse_locations(handle, bounds=desired_bounds, flip=False):
	return map(lambda line: ((lambda (x, y): (y, x)) if flip else (lambda x: x))(tuple(imap(float, line.replace('\t',',').split(',')))), handle)

def parse_jaccard(handle):
	return map(lambda line: map(float, line.split()), handle)


def parse_best_trains(handle):
	temp= map(lambda l :l.split(',') ,handle)
	return temp


def savefig(output_filename,m,title=""):
	#plt.title(title)
	if output_filename:
		plt.savefig(output_filename,bbox_inches='tight', dpi=200)
	else:
		print('The figure is shown...')
		plt.show()

def get_map(bounds=desired_bounds,showlines=True):
	m = Basemap(projection='cyl', resolution='i', **BOUNDS[desired_bounds])
	if showlines:
		m.drawcoastlines()
		m.drawstates()
		m.drawcountries()
	return m

def plot_map_relatedness(m,jaccard,positions, bounds=desired_bounds,t=threshold):
	"""
	m: can be obtained by calling the get_map method
	jaccard: a matrix
	positions: a list of (lon,lat) in degrees
	t: will draw the value contained in jaccard only if jaccard < t
	"""

	good_ones = [] # :: (distance, (start, end))
	greatest_distance = 0
	total = len(jaccard) * (len(jaccard) - 1)
	approved = 0
	denied = 0
	print(len(positions), len(jaccard))
	for (i, others) in enumerate(jaccard):
		position_i = positions[i]
		for (j, distance) in enumerate(others):
			if (approved + denied) % 5000 == 0:
				pass
				#print(approved + denied, '/', total, " (", approved, " approved, ", denied, " denied)", sep='')
			if i == j :
				continue
			if distance > t:
				denied += 1
				continue
			position_j = positions[j]
			good_ones.append( (distance, (position_i, position_j)) )
			if distance > greatest_distance:
				greatest_distance = distance
			approved += 1

	greatest_distance *= 1.05

	print("Done!")

	print("Plotting... ", end='')

	for (distance, (start, end)) in good_ones:
		if not (is_in_bounds(start, BOUNDS[desired_bounds]) and is_in_bounds(end, BOUNDS[desired_bounds])):
			print("skipping undesired line")
			continue
		m = plot_great_line(start, end ,m , color=(1., 0.5, 0.5, min(1-(distance / greatest_distance),1)))
	xs, ys = m(*zip(*filter(lambda x: is_in_bounds(x, BOUNDS[desired_bounds]), positions)))
	m.scatter(xs, ys,marker='o',s=10.0,color=(0.,0.,1.,1.))
	#plt.gcf().set_size_inches(11, 8.5)
	#plt.savefig("test.jpg", dpi=600)
	return m
	
def plot_great_line(start, end, m, color=(0.0,0.6,0.1,1)):
	"""
	Draws a line on the map between the points start and end.
	start : starting point in radiant [lon,lat]
	end : same as start
	"""
	#color=(0.27,0.,0.,1.)
	if start[0] > 179.8:
		start = (179.8, start[1])
	if end[0] > 179.8:
		end=(179.8, end[1])
	#print("start", start, "end", end)
	line, = m.drawgreatcircle(start[0], start[1], end[0], end[1], color=color,linewidth=0.6)
	p = line.get_path()
	# find the index which crosses the dateline (the delta is large)
	cut_point = np.where(np.abs(np.diff(p.vertices[:, 0])) > 200)[0]
	#print("CUTPOINT",cut_point,len(cut_point)>1)
	if len(cut_point) >= 1:
		cut_point = cut_point[0]
		print("YES")
		#cut_point = cut_point[0]
		# create new vertices with a nan in between and set those as the path's vertices
		if cut_point ==0 and 1==0:
			new_verts = np.concatenate([[np.nan,np.nan],p.vertices[1:,:]])
		else:
			new_verts = np.concatenate([p.vertices[:cut_point, :],
									[[np.nan, np.nan]],
									p.vertices[cut_point+1:, :]])
		p.codes = None
		p.vertices = new_verts
	return m

if __name__ == "__main__":
	args = sys.argv
	positions_filename  = ""
	jaccard_filename	= ""
	output_filename	 = ""
	more_lines_filename = ""
	flipyaxis=False

	try:
		i = 1
		while i < len(args):
			arg = args[i]
			if i + 1 < len(args):
				nextarg = args[i+1]

			if arg == "-p":
				positions_filename = nextarg
				i += 1
			elif arg == "-j":
				jaccard_filename = nextarg
				i += 1
			elif arg == "-o":
				output_filename = nextarg
				i += 1
			elif arg == "-b":
				desired_bounds = nextarg
			elif arg == "-m":
				more_lines_filename=nextarg
			elif arg == '-t':
				threshold = float(nextarg)
				i += 1
			elif arg == '-l':
				latlon = False
			elif arg == '-y':
				flipyaxis = True
			elif arg == "-h" or arg == "--help" or arg == "help":
				show_usage()
				sys.exit(1)
			i += 1
			
	except Exception as e:
		print("An error occurred while parsing commandline arguments:", e, file=sys.stderr)
		show_usage()
		sys.exit(1)

	if jaccard_filename == "":
		print("No Jaccard data specified.", file=sys.stderr)
		show_usage()
		sys.exit(1)
	elif positions_filename == "":
		print("No position data specified.", file=sys.stderr)
		show_usage()
		sys.exit(1)
	elif not desired_bounds in BOUNDS:
		print("Invalid bounds specified.", file=sys.stderr)
		show_usage()
		sys.exit(1)

	print("Parsing input files... ", end='')

	try:
		with open(positions_filename) as f:
			positions = parse_locations(f, bounds=desired_bounds, flip=False)
	except Exception as e:
		print("An error occurred while parsing the positions file:", e)
		sys.exit(2)

	try:
		with open(jaccard_filename) as f:
			jaccard = np.array(parse_jaccard(f))
			if flipyaxis:
				max = 150 #np.max(jaccard)
				jaccard = (jaccard - max)/(np.min(jaccard)-max)
				threshold =(threshold - max)/(np.min(jaccard)-max)
 
				#print jaccard[0,0]

	except Exception as e:
		print("An error occurred while parsing the Jaccard distance file:", e)
		sys.exit(2)
	if more_lines_filename != "":
		try:
			with open(more_lines_filename) as f:
				best_trains=parse_best_trains(f)
		except Exception as e:
			print("An error occurred while parsing the more line file:",e)


	print("Done!")

	print("Figuring out what to plot... ", end='')

	m=plot_map_relatedness(get_map(),jaccard,positions)
	if more_lines_filename != "":
		for item in best_trains[0:20]:
			plot_great_line(item[0],item[1],m ,(1.,0,1.,1.))
	savefig(output_filename,m)

	print("Done!")


