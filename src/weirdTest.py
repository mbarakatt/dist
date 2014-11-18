import numpy as np
from pylab import *
from scipy.stats import rankdata, spearmanr
import time
LARGE_NUMBER=iinfo('i').max


class weirdTest:
	def __init__(self,distIndiv,distRelated, threshold=0.22 , negateyaxis=False):
		self.sizeMatrix = distRelated.shape[0]
		self.distRelated_matrix = distRelated
		self.triu_indices = np.triu_indices(distRelated.shape[0],1)
		self.matrixDistJaccard = distRelated
		self.distRelated=distRelated[self.triu_indices]
		self.threshold=threshold
		if negateyaxis:
			self.threshold *= -1
			self.distRelated *= -1
		self.maskDistJaccard = (self.distRelated <= self.threshold)
		self.nb_y_below_t = np.sum(self.maskDistJaccard )
		#self.y_temp = rankdata((~self.maskDistJaccard) * LARGE_NUMBER + self.distRelated )
		self.sizexss = float(np.power(self.distRelated.size,3))
		######
		self.distIndiv=distIndiv[self.triu_indices]
		self.distargsort=np.argsort(self.distIndiv)
		self.sortdistIndiv=self.distIndiv[self.distargsort]
		self.sortdistRelated = self.distRelated[self.distargsort]
		#print np.unique(self.distIndiv),np.unique(self.distIndiv).size
		self.sortMaskDistIndiv=self.maskDistJaccard[self.distargsort]#List of bool in same order a sortdistIndiv that says if dist is above or below jaccard T.
		#self.sorty_temp=self.y_temp[self.distargsort]
		#print np.unique(self.y_temp),np.unique(self.sorty_temp)
		self.vert_line_index=self.get_line_index()
	
	#This function returns the index of the first and last time value appear in the list	
	def get_index_first_and_last_such_value(self,value):
		index=np.searchsorted(self.sortdistIndiv,value)
		number=np.sum(self.sortdistIndiv == value)
		return index,index + number -1
	
	def get_line_index(self):
		val=[0]
		for i in range(self.sizeMatrix-1):
			val.append(self.sizeMatrix - 1 - len(val)+val[-1])
		return np.array(val)
	def spearman(self,bestDistances):
		return spearmanr(bestDistances[self.triu_indices], self.distRelated )[0]

	#get_index_pts 
	def get_index_pts(self,i):
		vert_indices = (i-1) + self.vert_line_index[ 0:i ]
		hor_indices =  range(i + self.vert_line_index[i], self.sizeMatrix - 1 + self.vert_line_index[i])
		#print "Vert", vert_indices,"Hor",hor_indices
		#print "iVert" , np.argsort(self.distargsort)[vert_indices], "iHor", np.argsort(self.distargsort)[hor_indices]
		return np.concatenate([np.argsort(self.distargsort)[vert_indices] , [0] , np.argsort(self.distargsort)[hor_indices] ])

	def BL(self,a1,a2):
		temp= (self.sortdistRelated[a2:a1] < self.sortdistRelated[a1])
		return -np.sum(self.sortMaskDistIndiv[a2:a1] & temp) + np.sum(self.sortMaskDistIndiv[a2:a1] & ~temp) + np.sum(~self.sortMaskDistIndiv[a2:a1])

	def BR(self,a1,a2):
		temp=(self.sorty_temp[a1+1:a2+1] < self.sorty_temp[a1])
		return np.sum(self.sortMaskDistIndiv[a1+1:a2+1] & temp) - np.sum(self.sortMaskDistIndiv[a1+1:a2+1] & ~temp) - np.sum(~self.sortMaskDistIndiv[a1+1:a2+1])

		
	def TR(self,a1,a2):
		return np.sum(self.sortMaskDistIndiv[a1+1:a2+1]) 

	def TL(self,a1,a2):
		return -np.sum(self.sortMaskDistIndiv[a2:a1]) 

	def tempname(self,a1,a2):
		#print "between:", a1, a2, self.sorty_temp[a1+1:a2], self.sortMaskDistIndiv[a1+1:a2]
		if a1 == -1 : 
			a1 = 0
		if a2 == -1 : 
			a2 = 0

		if self.sortMaskDistIndiv[a1]:
			#Below
			if a1 < a2:
				return self.BR(a1,a2)	
			else:
				return self.BL(a1,a2)
		else:
			#Above
			if a1 < a2:
				return self.TR(a1,a2)
			else:
				return self.TL(a1,a2)
			

	def __call__(self,line,bestDistances,tookTrainMask):
		bestDistancesFlatten=bestDistances[self.triu_indices]
		retval=self.spear(bestDistancesFlatten)
		return retval

	def spear(self, xs):
		self.y_temp = rankdata(self.distRelated) + (~self.maskDistJaccard)*LARGE_NUMBER
		yRanks = np.amin([self.y_temp, rankdata((self.maskDistJaccard) * LARGE_NUMBER + xs)  + self.nb_y_below_t],axis=0)
		#print zip(yRanks[0::20],(self.y_temp)[0::20],(rankdata((self.maskDistJaccard) * LARGE_NUMBER + xs)  + self.nb_y_below_t)[0::20], rankdata(xs)[0::20])
		#print "numberptsbelowt: ", np.sum(self.maskDistJaccard)
		retval = 1 -np.sum(np.power(rankdata(xs) - yRanks,2))
		return retval


	

