import sys
import numpy as np
import getopt

def usage():
	print "Usage: \n\t CM_outfile NB_oufile < infile.txt"
	print "You entered:", sys.argv
try:
	poutcm=sys.argv[1]
	poutnb=sys.argv[2]
except:
	usage()
	exit(2)
foutcm = open(poutcm, 'w')
foutnb = open(poutnb, 'w')
list_indivs=open("../SCCS/list_indivs.txt",'r').read().split('\n')[:-1]

class lineO:
	def __init__(self,sp_line):
		self.indiv1 = indivs[sp_line[0]] 
		self.indiv2 = indivs[sp_line[2]]
		self.gen_len=float(sp_line[sp_line.index("cM")-1])

cm_Matrix=np.zeros([len(list_indivs)]*2)
nb_Matrix=np.zeros([len(list_indivs)]*2)
		
indivs={}

c=0
for i in list_indivs:
	indivs[i]=c
	c+=1

print "asdf"
c=0
for line in sys.stdin:
	sp=line.split()
	thisline=lineO(sp)
	if thisline.gen_len > 18:
		cm_Matrix[thisline.indiv1,thisline.indiv2]+=thisline.gen_len	
		cm_Matrix[thisline.indiv2,thisline.indiv1]+=thisline.gen_len	
		nb_Matrix[thisline.indiv1,thisline.indiv2]+=1
		nb_Matrix[thisline.indiv2,thisline.indiv1]+=1
	if c % 10000 == 0:
		print c
	c+=1

np.savetxt(foutcm,cm_Matrix,delimiter='\t')
np.savetxt(foutnb,nb_Matrix,delimiter='\t')

print "Done"
