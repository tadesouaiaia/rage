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

sns.set(color_codes=True)
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





		self.xLen, self.yLen = 5,5 
		self.xLoc, self.yLoc = 0,0 
		self.options = options

		if self.options and self.options.prefix != None: self.prefix = options.prefix
		else:						 self.prefix = 'model_plot'
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 3)
		self.pv_key = dd(int) 


	def add_pts(self,pts,SHOW='HEIGHT'): 


		for h,w,p,a in pts: 
			if p == 'F': color='blue'
			else: 	     color = 'red'

 
			if SHOW=='HEIGHT': self.ax.scatter(a,0,color=color) 
			elif SHOW=='BOTH': self.ax.scatter(w,h,color=color) 
			elif SHOW == 'DENSITY': self.ax.scatter(w/float(h),a,color=color) 

			else:   	self.ax.scatter(a,0,color=color) 

		if SHOW == 'BOTH': 
	
			self.ax.set_xlabel('weight') 
			self.ax.set_ylabel('height') 

		if SHOW != 'BOTH': 	
			self.ax.set_title('AGE')
			if SHOW == 'DENSITY': self.ax.set_title('DENSITY')  
			#self.ax.set_ylim(-0.1,0.1) 
			self.ax.axis('off')
 	
		else:
			 
			self.ax.set_xticks([]) 
			self.ax.set_yticks([]) 
		


		self.xLoc +=1
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 3)
		return 







		





	def save_mfig(self,mt,preds,covars,key={}):


		p_key, f_key, m_str, x = dd(int), {},  mt+':  y  \propto  ', 1
		for p in preds+covars:
			if p == 'intercept': m_str+=' B_o +'
			else:
				p_key[p.split('=')[0]]+=1
				f_key[p.split('=')[0]] = p 

		for i,(n,p) in enumerate(sorted([(b,a) for (a,b) in p_key.items()])):
			if n == 1: p = f_key[p] 
			if len(p.split('_')) > 1: p = '\_'.join(p.split('_')) 
			if i == 0:
				if n == 1:  m_str+= ' B_'+str(x)+'('+p+')'
				else:       m_str+= ' B_'+str(x)+'('+p+')_'+str(n)
			else:
				if n == 1:  m_str+= ' + B_'+str(x)+'('+p+')'
				else:       m_str+= ' + B_'+str(x)+'('+p+')_'+str(n)
			x+=1


 		plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.90,wspace=0.25,hspace=1.05)
		plt.suptitle('$'+m_str+'$',fontsize=25, fontweight='bold')
		if 'axis_off' in key: plt.axis('off') 

		fig_name = self.options.prefix+"_"+mt+"_predictorModelAnalysis_"+"_".join(self.options.predictors)+"_covars_"+"_".join(covars)+".pdf" 
		self.fig.savefig(str(fig_name),dpi=200)
	
	
		if self.options.show: 
			plt.show() 





def scan(opts): 

	pts = [] 
	bb,fb= 0,0 
	while True: 
		h,w = randrange(100,200), randrange(100,200)  
		d =  h/float(w) 
		fa,ba = 0.2+ random() * 0.8, 0.1+random() * 1.2

		if d < 0.99: 

			bb +=1

 
			if bb < 101:
				pts.append([h,w,'B',fa]) 
		elif d > 1.01: 
		  
			fb +=1
			if fb == 10:
				 
				pts.append([199.5,200,'B',fa]) 
 
			elif fb < 101:
				pts.append([h,w,'F',ba])


		if fb >= 101 and bb >= 101: 
			break 

	plot = mplot(3,3) 
	plot.add_pts(pts,SHOW='BOTH') 	
	plot.add_pts(pts) 	
	plot.add_pts(pts,SHOW='WEIGHT') 	
	plot.add_pts(pts,SHOW='DENSITY') 	
 	plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.90,wspace=0.25,hspace=1.05)
	plt.show() 
	























































































if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)


    parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
    parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

    (options, args) = parser.parse_args()



    

    scan(options)













