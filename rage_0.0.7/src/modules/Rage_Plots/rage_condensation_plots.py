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



class condense_plot:
        def __init__(self,xLen,yLen,options=None,progress=None):



      		seaborn.set(rc={'axes.facecolor':'black', 'figure.facecolor':'lightblue'})
	
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(18.5, 9.5)
		self.xLen, self.yLen = xLen,yLen
		self.xLoc, self.yLoc = 0,0 
		self.options = options
		self.axes = [] 



	def add_data(self,dim_members,dim_run,key={}):



		pts = dim_run['pts']
		axes_labels = dim_run['axes'] 
		
		x_comp, y_comp = self.xLoc, self.yLoc + 1
		legend = {} 
		while True:
			self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))
		
			pt_id, pt_color, pt_mark = '', 'w', 'o' 



			for p,m in zip(pts,dim_members):

			
				if 'color' in m.notes: 
					pt_id = m.notes['labels']
					pt_color = m.notes['color']

				
				if pt_id != 'NA':
					if 'size' in key:	self.axes[-1].scatter(p[x_comp],p[y_comp],color=pt_color,alpha=0.6,s=key['size'])
					elif m.notes['size']: 	self.axes[-1].scatter(p[x_comp],p[y_comp],color=pt_color,s=m.notes['size'],alpha=0.6)
					else:			self.axes[-1].scatter(p[x_comp],p[y_comp],color=pt_color,alpha=0.85)


					if pt_id: legend[pt_id[0]] = (pt_color, pt_mark) 
				
				if self.options.plotnames:
					if 'color' in m.notes: self.axes[-1].text(p[x_comp],p[y_comp],m.name.split(';')[-1],color='white') 
					else: self.axes[-1].text(p[x_comp],p[y_comp],m.name.split(';')[-1],color='cyan') 

			

			if 'zoom' in key: 
				pX = sorted([p[x_comp] for p in pts])
				pY = sorted([p[y_comp] for p in pts])
				xQ1,xQ3 = int(len(pX)*0.25), int(len(pX)*0.75)
				yQ1,yQ3 = int(len(pY)*0.25), int(len(pY)*0.75)
				iqX = pX[xQ3]-pX[xQ1]
				iqY = pY[yQ3]-pY[yQ1]  
				self.axes[-1].set_xlim(pX[xQ1]-(iqX*2.5),pX[xQ3]+(iqX*2.5))
				self.axes[-1].set_ylim(pY[yQ1]-(iqY*2.5),pY[yQ3]+(iqY*2.5))

			self.axes[-1].set_xlabel(axes_labels[x_comp],fontweight='bold') 	
			self.axes[-1].set_ylabel(axes_labels[y_comp],fontweight='bold')
			self.axes[-1].set_xticks([]) 
			self.axes[-1].set_yticks([]) 
			if self.yLoc +1 == self.yLen: self.xLoc, self.yLoc = self.xLoc+1, 0 
			else: 			      self.yLoc +=1



			if y_comp == x_comp+1:   y_comp+=1
			else: 		    x_comp+=1 
#			x_comp, y_comp = x_comp +2, y_comp +2

			if self.xLoc == self.xLen: 
				break 

		if 'title' in key: plt.suptitle(key['title'],fontsize=20,fontweight='bold')
		items,labels  = [Rect((0,0),1,1,fc=b[0]) for a,b in legend.items() if a!='NA'],[a for a in legend.keys() if a != 'NA'] 

		self.axes[0].legend(items,labels,loc = 'upper left', ncol = len(legend.values()), bbox_to_anchor=(1.7,1.1)) 

		if 'out' in key: self.fig.savefig(key['out'], dpi=300)
		if self.options.show: plt.show()
		return self




		

		

























		

	def update(self,key={}):

		if 'clear_axes' in key: 
			self.ax.set_xticks([]) 
			self.ax.set_yticks([]) 

		


		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])
		if 'title' in key: 	self.ax.set_title(key['title'],fontweight='bold')
		self.yLoc += 1 
		if self.yLoc == self.yLen: 
			self.xLoc,self.yLoc = self.xLoc+1, 0 
		if self.xLoc >= self.xLen: return False
		
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		return True 



	def save(self,f_name,key={}):

		if 'axis_off' in key: plt.axis('off') 
		if 'title' in key: plt.suptitle(key['title'],fontsize=20,fontweight='bold')
		elif self.options and 'title' in vars(self.options).keys() and self.options.title != None:	plt.suptitle(self.options.title)
               	self.fig.savefig(f_name, dpi=100) 

		if self.options.show: 
			plt.show() 

		sys.exit() 

	def add_hist(self,h_data):
		a = np.hstack(h_data)
		if len(self.ax.hist(a,bins='auto')[0]) < 2:
			self.ax.clear() 
			self.ax.hist(a,bins=min(4,len(cc(h_data)))) 
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

	def add_line(self,x,y,key={}):
		lw,alpha = 1.0, 1.0 
		if 'lw' in key.keys():    lw = key['lw'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'color' in key: self.ax.plot(x,y,linewidth=lw,color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.plot(x,y,linewidth=lw,zorder=1,alpha=alpha)
		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])


	def scatter_pts(self,X,Y,key={}):
		alpha,size,mark =  1.0, 20,'o'
		if 'size' in key.keys():    size = key['size'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'mark' in key: mark= key['mark'] 
		if 'yjitter' in key: Y = [np.random.normal(y,(0.0195)) for y in Y]
		

		if 'color' in key: self.ax.scatter(X,Y,s=size,marker = mark, color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.scatter(X,Y,zorder=1,alpha=alpha,s=size)


	def add_scatter_data(self,x,y,xlabel=None,ylabel=None,clr=None):
		if clr and len(clr.split(',')) == 2: 
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),scatter_kws={"color": clr.split(',')[0]}, line_kws={"color": clr.split(',')[1]})
		else:
			if clr == None: clr = 'k'
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),color=clr)
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

	
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
		for a,b in zip(labels,colors):
			labs.append(a) 
			items.append(Rect((0,0),1,1,fc=b))

		plt.legend(items,labs,loc = 'upper left',ncol = len(labs)/3,bbox_to_anchor=(-0.1,1.16),fontsize=10)
		

		
