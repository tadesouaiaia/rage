#!/usr/bin/env python

import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
import scipy.stats as stats
from scipy.stats import variation as coVar 

from random import random
import numpy as np
from random import randrange
from math import fabs
from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
import seaborn
from math import log
import math
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from random import shuffle
from sklearn.cluster import KMeans	
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import MinMaxScaler
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle as Rect
from matplotlib.lines import Line2D
import os
import sys
from collections import defaultdict as dd
from collections import Counter as cc
from scipy import stats
from math import log 
from scipy.stats import chisquare
import numpy as np 
fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}




#matplotlib.rcParams['xtick.labelsize'] = 8.5
matplotlib.rcParams['ytick.labelsize'] = 6.5


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))




def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation

#    print get_colors(modified_z_score, plt.cm.jet)
    return modified_z_score > thresh, modified_z_score, get_colors(modified_z_score, plt.cm.seismic) 
#    return modified_z_score > thresh, modified_z_score, get_colors(modified_z_score, plt.cm.jet) 


















class mplot:
        def __init__(self,xLen,yLen,key={}):



		
		if 'r_key' in key: self.r_key = key['r_key']
		if 'p_key' in key: self.p_key = key['p_key']
	
		sns.set(rc={'axes.facecolor':'lightgrey', 'figure.facecolor':'white'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
		self.fig.set_facecolor('skyblue') 
    		self.fig.patch.set_facecolor('skyblue')
    		matplotlib.rcParams['savefig.facecolor'] = 'white'
		matplotlib.rcParams['ytick.labelsize'] = 7.5
		
		seaborn.set(rc={'axes.facecolor':'lightgrey', 'figure.facecolor':'white'})

                self.fig.patch.set_facecolor('white')
		self.xLen, self.yLen = xLen,yLen 
		self.xLoc, self.yLoc = 0,0 
		self.options = options

		if self.options and self.options.prefix != None: self.prefix = options.prefix
		else:						 self.prefix = 'model_plot'
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		self.pv_key = dd(int) 


	def add_pts(self,pts,names,coefs,varE,cords=[0,1]): 

		x,y = cords[0], cords[1]
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		for i,p in enumerate(pts): 
			

			size = 5 
			if len(names[i]) == 4: 
				name,current,order,spikes = names[i] 
			else:
				name,order,current,spikes,activity = names[i] 
			if name[0:2] == 'EB':   	color='red'
			elif name[0:2] == 'ES': 	color = 'cyan'
			elif name[0:1] =='T': 			    	color = 'green'
			elif name[0:2] == 'OB': 				color == 'purple'
			else: color= 'k' 

			text = name.split('~')[0]+'_'+str(order+1)+'/'+str(spikes)
			self.ax.text(p[x],p[y],text,color=color,fontsize=9,alpha=0.75) 
			
			#self.ax.text(p[x],p[y],text,color=color,fontsize=9,alpha=0.75) 
			size = 5 + int(current / 12.0)
			self.ax.scatter(p[x],p[y],color=color,s=size) 

		self.ax.set_xlabel('PCA'+str(x+1)+' = '+str(round(varE[x],3)))
		self.ax.set_ylabel('PCA'+str(y+1)+' = '+str(round(varE[y],3)))
			
		self.ax.set_xticks([]) 
		self.ax.set_yticks([]) 

 
		self.label_axis(coefs[x]) 
		self.label_axis(coefs[y],'Y') 



		self.xLoc +=1
		return 



	def label_axis(self,coefs,axis='X'): 
		x_pos =  sorted([(c[0],c[-1]) for c in coefs if c[1] > 0], reverse=True) 
		x_neg =  sorted([(c[0],c[-1]) for c in coefs if c[1] < 0], reverse=True) 
	
		gMax = 50.0
	
		x0,x1 = self.ax.get_xlim() 
		y0,y1 = self.ax.get_ylim() 
		
		yStep = (y1-y0)/gMax
		xStep = 2*((x1-x0)/gMax)
		fMax = 12
		fReal = 11
		fMin  = 6 
		maxItems = 100
		if axis.upper() != 'Y':

			

			gInit,xInit,yInit,k = gMax-3,x1,y1,0   
			for xS,xN in x_pos: 
				yInit-=yStep	
				fS = xS/float(x_pos[0][0])
				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit += (gMax-3)
					xInit += xStep
					yInit = y1-yStep 
					
 
				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin

				if fS > 0.90:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='red',verticalalignment='top',horizontalalignment='left') 
				elif fS > 0.80:	self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='orange',verticalalignment='top',horizontalalignment='left') 
				elif fS > 0.60:	self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='yellow',verticalalignment='top',horizontalalignment='left') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='green',verticalalignment='top',horizontalalignment='left') 
				else:		self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='purple',verticalalignment='top',horizontalalignment='left') 


			gInit,xInit,yInit,k = gMax-3,x0,y1,0   
			for xS,xN in x_neg: 
				yInit-=yStep					
				fS = xS/float(x_neg[0][0])
				 

				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit+=(gMax-3)
					xInit -= (1.1*xStep) 
					yInit = y1-yStep 	
				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin
				
				
				if fS > 0.85:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='k',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.65:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='grey',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='purple',verticalalignment='top',horizontalalignment='right') 
				else:		self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='blue',verticalalignment='top',horizontalalignment='right') 

				


		else: 
			
			gInit,xInit,yInit,k = gMax*1.1,x0,y1,0 
			for xS,xN in x_pos: 
				xInit+=(yStep*1.3)	
				fS = xS/float(x_pos[0][0])
				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit += gMax
					xInit = x0+(2*yStep)
					yInit = y1+yStep 
					

				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin

				if fS > 0.85:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='k',verticalalignment='bottom',horizontalalignment='right') 
				elif fS > 0.65:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='grey',verticalalignment='bottom',horizontalalignment='right') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='purple',verticalalignment='bottom',horizontalalignment='right') 
				else:		self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='blue',verticalalignment='bottom',horizontalalignment='right') 

			gInit,xInit,yInit,k = gMax*1.1,x0,y0,0 
			for xS,xN in x_neg: 
				xInit+=(yStep*1.3) 	
				fS = xS/float(x_neg[0][0])
				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit += gMax
					xInit =  x0+(2*yStep) 
					yInit = y0-yStep 
				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin

				if fS > 0.85:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='k',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.65:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='grey',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='purple',verticalalignment='top',horizontalalignment='right') 
				else:		self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='blue',verticalalignment='top',horizontalalignment='right') 







































































def qc(name): 
	if name[0:2] == 'EB': return 'red'
	elif name[0:2] == 'ES': return 'cyan'
	elif name[0:2] == 'OB': return 'purple'
	elif name[0] == 'T': return 'green'
	else:			return 'k'

def run_pca(data,names): 

	#vals = [[log(x+1.0) for x in sv] for sv in s_vals] 
	fit_mat = np.matrix(data) 
	fit_pca = PCA(n_components = 5).fit(fit_mat)
	varE =  fit_pca.explained_variance_ratio_ 
	pts = fit_pca.transform(fit_mat)
	coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,range(len(comps)))],reverse=True) for comps in fit_pca.components_]



	plot = mplot(1,1) 
	plot.add_pts(pts,names,coefs,varE,cords=[0,1]) 	

	return pts,coefs












def scan_key(spike_file,options):

	pkey = {} 
	skey = {} 
	nkey = {}
	nick =  spike_file.split('@')[0] 
 
	for line in open(spike_file):
		line = line.split() 	
		if line[1] == 'vmon1': 
			pkey[int(line[0])] = [float(x) for x in line[3::]]
		elif line[1] == 'features' and line[2] == 'spikes' and len(line)>3:
			skey[int(line[0])] = int(line[3].split(',')[3])
			nkey[int(line[0])] = len(line) - 3 		



	ki = 1
	for k in sorted(skey.keys()):
		sp = skey[k] 	
		pts = pkey[k]
		ht = pts[sp] 
		offset = 100-ht 

		print nick+'@'+str(ki)+'@0@'+str(nkey[k]),",".join([str(round(pts[i]+offset,1)) for i in range(sp-299,sp+301)])

		ki+=1
#		plt.plot(range(10000),pts[0:10000]) 
#		plt.scatter(sp,pts[sp],color='red',marker='*',s=100) 
#		plt.show() 
		 






				






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-m", "--masks", default = None, type='string', help="Output Filename Prefix")
	parser.add_option( "--show", default = False,action='store_true' , help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	scan_key(args[0],options)













