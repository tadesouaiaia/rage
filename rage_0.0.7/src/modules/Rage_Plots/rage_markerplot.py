#!/usr/bin/env python
import pandas as pd
import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
import scipy.stats as stats
from scipy.stats import variation as coVar 

from random import random
import numpy as np

import random
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
sns.set(color_codes=True)



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




class subplot:
        def __init__(self,xLen,yLen,options=None):

	
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(18.5, 9.5)
		self.xLen, self.yLen = xLen,yLen
		self.xLoc, self.yLoc = 0,0 
		self.options = options
		self.trend = False



		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		

	def update(self,key={}):

		if 'clear_axes' in key: 
			self.ax.set_xticks([]) 
			self.ax.set_yticks([]) 

		


		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])

		if not self.trend:
			if 'title' in key: 	self.ax.set_title(key['title'],fontweight='bold',y=1.05,x=0.0,horizontalalignment='left')


		self.yLoc += 1 
		if self.yLoc == self.yLen: 
			self.xLoc,self.yLoc = self.xLoc+1, 0 
		if self.xLoc >= self.xLen: return False
		
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)


		self.trend=False

		return True 



	def save(self,f_name,key={}):

		if 'axis_off' in key: plt.axis('off') 
		if 'title' in key: plt.suptitle(key['title'],fontsize=20,fontweight='bold')
		elif self.options and 'title' in vars(self.options).keys() and self.options.title != None:	plt.suptitle(self.options.title)
               	self.fig.savefig(f_name, dpi=100) 

		if self.options.show: 
			plt.show() 



	def add_hist(self,h_data):


		x,y = h_data.keys(), sorted(h_data.values())
	#	a = np.hstack(h_data)
		if len(self.ax.hist(y,bins='auto')[0]) < 2:
			self.ax.clear() 
			self.ax.hist(np.hstack(y),bins=min(4,len(cc(y)))) 
		yMin, yLim = self.ax.get_ylim()
		xMin, xLim = self.ax.get_xlim()
		HI,LO=False,False
		if yMin == 0: 
			h_srt = sorted(y) 
			out_bool, out_scrs, out_colors =  mad_based_outlier(np.array(h_srt)) 
			g_color = out_colors[sorted([(out_scrs[oi],oi) for oi in range(len(out_scrs))])[0][1]]
			for i,h in enumerate(h_srt): 
				if not out_bool[i]:
					self.ax.scatter(h,yLim*1.025,alpha=0.7,color=g_color,clip_on=False) 
				else: 
					self.ax.scatter(h,yLim*1.025,alpha=0.7,color=out_colors[i],clip_on=False) 
				#	if i > 5 and (out_scrs[i]/(out_scrs[i-1]+out_scrs[i-2])) > 0.66:  HI=True 
				#	if HI: 			self.ax.text(h,yLim*1.030,x[i].name.split(";")[-1])


			self.ax.set_ylim(yMin,yLim) 
			self.ax.set_xlim(xMin,xLim) 
		return self




	def add_pca_data(self,pca_pts,key={}):

		#sns.set(rc={'axes.facecolor':'black', 'figure.facecolor':'cornflowerblue'})
		if 'components' in key: comps = key['components']
		else:			comps = [0,1] 
		if 'type' in key:
			if key['type'] == 'tsne':	
				axis_labs,title = ['TS'+str(c+1) for c in comps], 'TSNE'
			elif key['type'] == 'kca': 
				axis_labs,title = ['KC'+str(c+1) for c in comps], 'KCA'
			elif key['type'] == 'mds': 
				axis_labs,title = ['MD'+str(c+1) for c in comps], 'MDS'
			else:
				axis_labs,title = ['X'+str(c+1) for c in comps], key['type'].upper() 
		else:
			title = 'PCA'
			if 'vars' in key: 				axis_labs = ['PC'+str(c+1)+':  '+str(round(100*key['vars'][c],2))+'%' for c in comps]
			else:						axis_labs = ['PC'+str(c+1) for c in comps]

		for i,p in enumerate(pca_pts):
			if 'colors' in key and key['colors'][i] != 'white': 
				if 'sizes' in key: self.ax.scatter(p[0],p[1],s = key['sizes'][i],color=key['colors'][i],alpha=0.6)
				else: 		   self.ax.scatter(p[0],p[1],color=key['colors'][i])
			else: 		    self.ax.scatter(p[0],p[1])

		self.ax.set_xlabel(axis_labs[0]) 	
		self.ax.set_ylabel(axis_labs[1])

		if 'zoom' in key: 
			pX = sorted([p[0] for p in pca_pts])
			pY = sorted([p[1] for p in pca_pts])
			xQ1,xQ3 = int(len(pX)*0.25), int(len(pX)*0.75)
			yQ1,yQ3 = int(len(pY)*0.25), int(len(pY)*0.75)
			iqX = pX[xQ3]-pX[xQ1]
			iqY = pY[yQ3]-pY[yQ1]  
			self.ax.set_xlim(pX[xQ1]-(iqX*2.5),pX[xQ3]+(iqX*2.5))
			self.ax.set_ylim(pY[yQ1]-(iqY*2.5),pY[yQ3]+(iqY*2.5))



		if 'title' in key: self.ax.set_title(key['title']) 
		else:		   self.ax.set_title(title) 
		return self		
		


	def change_limits(self,key):
		x0,x1 = self.ax.get_xlim()
		y0,y1 = self.ax.get_ylim()
		if 'x0' in key: x0 = key['x0'] 
		if 'x1' in key: x1 = key['x1'] 
		if 'y0' in key: y0 = key['y0'] 
		if 'y1' in key: y1 = key['y1'] 
		self.ax.set_xlim(x0,x1)
		self.ax.set_ylim(y0,y1)


	def add_lines(self,x,y,xlabel=None,ylabel=None,clr=None):
		if not clr: self.ax.plot(x,y,zorder=1)
		else: 	    self.ax.plot(x,y,color=clr,zorder=1) 
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 
		return self

	def add_line(self,x,y,key={}):
		lw,alpha = 1.0, 1.0 
		if 'lw' in key.keys():    lw = key['lw'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'color' in key: self.ax.plot(x,y,linewidth=lw,color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.plot(x,y,linewidth=lw,zorder=1,alpha=alpha)
		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])

		return self

	def scatter_pts(self,X,Y,key={}):
		alpha,size,mark =  1.0, 20,'o'
		if 'size' in key.keys():    size = key['size'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'mark' in key: mark= key['mark'] 
		if 'yjitter' in key: Y = [np.random.normal(y,(0.0195)) for y in Y]
		

		if 'color' in key: self.ax.scatter(X,Y,s=size,marker = mark, color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.scatter(X,Y,zorder=1,alpha=alpha,s=size)
		return self

	def add_scatter_trend(self,x,y,xlabel=None,ylabel=None,clr=None):



		if clr and len(clr.split(',')) == 2: 
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),scatter_kws={"color": clr.split(',')[0]}, line_kws={"color": clr.split(',')[1]})
		else:
			if clr == None: clr = 'k'
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),color=clr)






#		 pd.DataFrame(np.array([np.array([values[j][i] for j in range(len(values))]) for i in range(len(idxs))]), columns = labels)

#		self.ax = sns.lmplot(x="x", y="y", data=data);

		return self


	def add_scatter_data(self,x,y,xlabel=None,ylabel=None,clr=None):




		if clr and len(clr.split(',')) == 2: 
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),scatter_kws={"color": clr.split(',')[0]}, line_kws={"color": clr.split(',')[1]})
		else:
			if clr == None: clr = 'k'
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),color=clr)
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

		return self
	
	def add_scatter_pts(self,X,Y,notes='None'):
		notes = notes.split(',') 
		if 'RED' in notes: clr = 'red' 
		elif 'PURPLE' in notes: clr = 'purple' 
		elif 'BLUE' in notes: clr = 'blue'
		else:		      clr = 'k'  
		if 'XY_JITTERS' in notes:
			xJ = [np.random.normal(x,(0.0000001+x*0.1)) for x in X]
			yJ = [np.random.normal(y,(0.0000001+y*0.2)) for y in Y]
			self.ax.scatter(xJ,yJ,color=clr)

		else:
			self.ax.scatter(X,Y,color=clr) 


	




	def add_labels(self,title,xlabel,ylabel):
		self.ax.set_title(title) 
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

	def add_legend(self,labels,colors):
		labs,items = [],[]

		print 'yo'  
		for a,b in zip(labels,colors):
			labs.append(a) 
			items.append(Rect((0,0),1,1,fc=b))

		plt.legend(items,labs,loc = 'upper left',ncol = len(labs)/1.0,bbox_to_anchor=(-0.1,1.16),fontsize=10)
		
	def add_skree(self,pca_vars,key={}):

		if 'color' in key: 	
			self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))],color=key['color'])
		else:
			self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))])

		return self	



		
