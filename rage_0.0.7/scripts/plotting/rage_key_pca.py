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
	
		sns.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'cornflowerblue'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
		self.fig.set_facecolor('skyblue') 
    		self.fig.patch.set_facecolor('skyblue')
    		matplotlib.rcParams['savefig.facecolor'] = 'skyblue'
		matplotlib.rcParams['ytick.labelsize'] = 7.5
		
		seaborn.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'cornflowerblue'})

                self.fig.patch.set_facecolor('skyblue')
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
		

			if names[i][0:2] == 'EB':   color='red'
			elif names[i][0:2] == 'ES': color = 'cyan'
			else: 			    color = 'green' 
			self.ax.scatter(p[x],p[y],color=color) 

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
		
		x0,x1 = self.ax.get_xlim() 
		y0,y1 = self.ax.get_ylim() 
		yStep = (y1-y0)/7.0
		xStep = (x1-x0)/7.0
		fMax = 25
		fReal = 16 
		if axis.upper() != 'Y': 
			yInit,fInit,k = y1,x_pos[0][0],0   

			for xS,xN in x_pos[0:5]: 
				yInit-=yStep	
				fS = xS/float(fInit)
				fSize = int(fS*fMax) 
				if fSize > 18: fSize = fReal
				if k > 6: break 
				elif k > 3 and fSize < 7: break 
				elif fSize < 8: fSize = 8 
				self.ax.text(x1,yInit,xN,fontweight='bold',fontsize=fSize,verticalalignment='top',horizontalalignment='left') 

			yInit,fInit,k = y1,x_neg[0][0],0
			for xS,xN in x_neg[0:5]: 
				yInit-=yStep	
				fS = xS/float(fInit)
				fSize = int(fS*fMax) 
				
				if fSize > 18: fSize = fReal 
				if k > 6: break 
				elif k > 3 and fSize < 7: break 
				elif fSize < 8: fSize = 8 
				self.ax.text(x0,yInit,xN,fontweight='bold',fontsize=fSize,verticalalignment='top',horizontalalignment='right') 

		else: 
			xInit,fInit,k = x0,x_pos[0][0],0 
			for xS,xN in x_pos[0:5]: 
				xInit+=yStep	
				fS = xS/float(fInit)
				fSize = int(fS*fMax) 
				if fSize > 18: fSize = fReal
				if k > 6: break 
				elif k > 3 and fSize < 7: break 
				elif fSize < 8: fSize = 8 
				self.ax.text(xInit,y1,xN,rotation=30,fontweight='bold',fontsize=fSize,verticalalignment='bottom',horizontalalignment='left') 

			xInit,fInit,k = x0,x_neg[0][0],0  
			for xS,xN in x_neg[0:5]: 
				xInit+=yStep	
				fS = xS/float(fInit)
				fSize = int(fS*fMax) 
				if fSize > 18: fSize = fReal 
				if k > 6: break 
				elif k > 3 and fSize < 7: break 
				elif fSize < 8: fSize = 8 
				self.ax.text(xInit,y0,xN,rotation=-30,fontweight='bold',fontsize=fSize,verticalalignment='top',horizontalalignment='left') 




 
				#print xS,N
			


 
				#print xS,N
			












































































def scan_key(key_file,options):
	ktype = {} 
	kidx  = {} 
	names = [] 
	for line in open(key_file):
		line = line.split()
		if line[0] == '---':
			cats = line[1::] 
			key  = dd(list)
		else:
			#if line[0][0:2] != 'EB': continue 
			data = line[1::]
			names.append(line[0]) 
			
	    		for i in range(len(data)): 
				key[cats[i]].append(data[i]) 

	scaler = MinMaxScaler() 				
	mod_names, mod_data = [], []  


	for i in range(len(cats)):
		ktype[cats[i]],cd  = 'continuous', key[cats[i]] 

		kLen = len([x for x in key[cats[i]] if x != 'NA'])
		if cats[i] == 'REB_TYPE': continue  
		if 'ISI' in cats[i].split('_'): continue 
		#if kLen / float(len(key[cats[i]])) < 0.85: continue 
		if kLen / float(len(key[cats[i]])) < 0.5: continue 

		try: 		
			c_data = [] 	
			key[cats[i]] = [float(x) if x!= 'NA' else 'NA' for x in key[cats[i]]] 
			cn = [float(x) for x in cd if x != 'NA'] 
			scaler.fit(np.array(cn).reshape(-1,1))  
			for j,v in enumerate(cd): 
				if v == 'NA': 
					c_data.append(0.5) 
				else:
					c_data.append(scaler.transform(np.array([v]).reshape(1,-1))[0][0])

			mod_names.append(cats[i])
			mod_data.append(c_data) 
		except ValueError: 		
			ktype[cats[i]] = 'binary' 
			opts =  cc(cd)
			mod_opts =  [x[0] for x in sorted(opts.items(),key=lambda X: X[1]) if (x[0] != 'NA' and opts[x[0]] > 20) ] 
			if len(mod_opts) == len(opts): mod_opts = mod_opts[1::] 
			for m in  mod_opts: 
				mod_names.append(cats[i]+'='+m) 
				mod_data.append([1.0 if v == m else 0 for v in cd])
				
	
	fit_mat = np.matrix(mod_data) 
	fit_pca = PCA(n_components = 5).fit(fit_mat.getT())

	varE =  fit_pca.explained_variance_ratio_
		 
	pts = fit_pca.transform(fit_mat.getT())

	coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,mod_names)],reverse=True) for comps in fit_pca.components_]

	plot = mplot(2,1) 
	plot.add_pts(pts,names,coefs,varE,cords=[0,1]) 	
	plot.add_pts(pts,names,coefs,varE,cords=[2,3]) 	
	plt.subplots_adjust(left=0.15, bottom=0.25, right=0.85, top=0.75,wspace=0.5,hspace=2.5)
	plt.show() 





		


				






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()


	scan_key(args[0],options)













