import sys
from pylab import *
import numpy as np

def usage():
	print "python scattergendist.py xs.txt (Matrix from tab delimiter) ys.txt (Same) [-o outfile.jpg]" 

args=sys.argv		

try:
	f_gendist=open(args[2],'r')
	f_geodist=open(args[1],'r')
except:
	usage()
	exit(2)

try:
	outpic = args[args.index("-o")+1]
except IOError as e:
	raise(e)
except:
	outpic=""


def parse_file_matrix(filep):
	temp=[]
	c=0
	for line in filep:
		temp.append(map(float,line.split()))
		sys.stdout.write("\r%d" %c)
		sys.stdout.flush()
		c+=1
	return np.concatenate(temp)

sys.stdout.write("\n")
gendist_data=parse_file_matrix(f_gendist)
geodist_data=parse_file_matrix(f_geodist)

plot(geodist_data,gendist_data,'.')
ylim(0,150)
xlim(0,1400)
xlabel(args[2].split('/')[-1])
ylabel(args[1].split('/')[-1])
#title('With  5 waypoints.')
if outpic != "":
	savefig(outpic)
else:
	show()



