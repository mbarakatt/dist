from zipdist import *

ZC=Zip_Codes("../data/zipcode/zipcode.csv")
print "Done loading the zip codes"

print ZC.get_distance("07304","30006")


findiv_pos=open("/lb/project/gravel/data/SCCS/crossref/ors_124_115.txt",'r')

indiv_pos=[]
for line in findiv_pos:
	sp = line.split()
	indiv_pos.append((float(sp[4]),float(sp[5])))