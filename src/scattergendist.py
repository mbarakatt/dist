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

try:
	myxlim = map(float,args[args.index("--xlim")+1].split(','))
except:
	myxlim=(0,0)


try:
	myylim = map(float, args[args.index("--ylim")+1].split(','))
except:
	myylim=(0,0)

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
gendist_data=np.loadtxt(f_gendist)
geodist_data=np.loadtxt(f_geodist,delimiter='\t')

print "Sizes:" , gendist_data.shape , geodist_data.shape
plot(geodist_data,gendist_data,'.')

#
if False:
	bins=np.linspace(0,1800,num=18+1)
	ibd=np.histogram(geodist_data,bins,weights=gendist_data)[0]
	counts=np.histogram(geodist_data,bins)[0]
	print ibd/counts 
	print f
	plot(ibd/counts,(bins+50)[:-1])

if myxlim != (0,0):
	xlim(*myxlim)

if myylim != (0,0):
	ylim(*myylim)

xfilename=args[1].split('/')[-1]
yfilename=args[2].split('/')[-1]
xlabel(xfilename)
ylabel(yfilename)
title (yfilename + " vs " + xfilename)
triu = np.triu_indices(len(gendist_data),k=1)
print geodist_data[triu]
scatter(geodist_data[triu],gendist_data[triu])
#title('With  5 waypoints.')
print "Saving figure"

if outpic != "":
	savefig(outpic)
else:
	show()



