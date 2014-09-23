from geode import *
import numpy as np

fdist=open("geographicGeoDistLanguagesTest.txt",'r')
distIndiv=parse_file_geo_dist(fdist)
fout=open("test_jaccard.txt",'w')

def get_jaccard(distM):
	distM[distM<100]=100
	f=lambda x : np.sqrt(x)/700.0+0.44
	#(distM - 100)/( 20000 - 100)*(0.6 - 0.5) + 0.5 
	jac=np.random.normal(f(distM),0.11)
	jac[jac<0]=0
	jac[jac>1]=1
	jac=np.triu(jac,k=1)
	return jac + jac.T

jaccardM= get_jaccard(distIndiv)
np.savetxt(fout,jaccardM,delimiter="\t")
