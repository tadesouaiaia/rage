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



class dimplot:
        def __init__(self,options=None,progress=None,xLen=1,yLen=1):





      		#seaborn.set(rc={'axes.facecolor':'black', 'figure.facecolor':'lightblue'})
      		seaborn.set(rc={'axes.facecolor': options.plotColors['AXIS'],  'figure.facecolor':options.plotColors['FACE']})



	
                self.fig = matplotlib.pyplot.gcf()
                #self.fig.patch.set_facecolor('thistle')
                self.fig.patch.set_facecolor(options.plotColors['FACE']) 
               	matplotlib.rcParams['savefig.facecolor'] = options.plotColors['FACE']
                matplotlib.rcParams['ytick.labelsize'] = 7.5



                self.fig.set_size_inches(18.5, 9.5)
		self.xLen, self.yLen = 1,1
		self.xLoc, self.yLoc = 0,0 
		self.options = options
		self.axes = [] 






	def add_data(self,dim_members,dim_runs,dim_comps=[],NAMES=False,NAMEOUTLIERS=False):



		if len(dim_comps) == 0:	dim_comps = [(0,1) for d in dim_runs]
		elif len(dim_runs) == 1 and len(dim_comps) > 1: dim_runs += [dim_runs[0] for d in range(len(dim_comps)-len(dim_runs))] 
		
			

		self.xLen, self.yLen, drs,dcs = 1,1,[],[] 
		for i,d in enumerate(dim_runs):
			if NAMES: drs.append([d,True]) 
			else:     drs.append([d,False]) 	
			dcs.append(dim_comps[i]) 
		
		if len(drs) > 1:
			self.xLen = len(drs) / 2
			if len(drs) < 6:   self.yLen = 2
			elif len(drs) < 9: self.yLen = 3 
			elif len(drs) == 9: 
				self.xLen,self.yLen = 3,3
			else: 		   
				self.yLen = 4 


		if self.xLen * self.yLen < len(drs): self.yLen +=1	
		for dI,(dim_run,SHOWNAMES) in enumerate(drs): 

			pts = dim_run['pts']
			axes_labels = dim_run['axes']
			x_comp, y_comp = dcs[dI][0], dcs[dI][1]


			outlier_key = dd(lambda: dd(list)) 
			legend = {}
			items,labels,OBS,parents,children = [],[],dd(bool),[],[] 
			while True:
				self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))
				for p,m in zip(pts,dim_members):

					m_labs = ['M'] 

					if m.label != None:
						self.axes[-1].plot(p[x_comp],p[y_comp],alpha=0.8,**m.label.vals)
						m_labs.append(m.label.vals['color']) 
						m_labs.append(m.label.vals['marker']) 
						if 'markerfacecoloralt' in m.label.vals: m_labs.append(m.label.vals['markerfacecoloralt'])


						if not OBS[m.label.id]:
							sLab = m.label.id
							sParents  = [sl.split('~')[0].split('=')[0] for sl in sLab.split('|')]
							sChildren = [sl.split('~')[-1].split('=')[-1] for sl in sLab.split('|')]
							sChildren = [",".join(sc.split(',')[0:2])+'...' if len(sc.split(','))>2 else sc for sc in sChildren ]

							parents.append(' & '.join(sParents))
							labels.append('\n'.join(sChildren))
							items.append(m.label.proxy)
							OBS[m.label.id] = True

					else:	
						self.axes[-1].scatter(p[x_comp],p[y_comp],color='silver',alpha=0.8)
					
					if self.options.plotnames or SHOWNAMES:

						
						#if m.label.id.split('|')[0].split('~')[-1] not in ['CBC','STR','MD','HIP']:  continue


 
						if 'color' in m.notes: 	self.axes[-1].text(p[x_comp],p[y_comp],m.name.split(';')[-1],color='black',fontsize=6,alpha=0.7) 
						else: 			self.axes[-1].text(p[x_comp],p[y_comp],m.name.split(';')[-1],color='black',fontsize=6,alpha=0.7) 

					elif NAMEOUTLIERS:
						m_lab_str = ','.join(m_labs) 
						sx,sy = p[x_comp],p[y_comp] 			
						outlier_key[m_lab_str]['x'].append(sx)
						outlier_key[m_lab_str]['y'].append(sy)
						outlier_key['ALL']['x'].append(sx)
						outlier_key['ALL']['y'].append(sy)
						outlier_key['samples'][m.name.split(';')[-1]] = (sx,sy,m_lab_str)


				
				if NAMEOUTLIERS: 
					GLOBAL_HI,LOCAL_HI = 70,95
					extrema = dd(lambda: dd(list))  
					for k in [kk for kk in outlier_key.keys() if kk != 'samples']:
						x_srt, y_srt = sorted(outlier_key[k]['x']), sorted(outlier_key[k]['y'])
						extrema['means']['x'].append((np.mean(x_srt),k))
						extrema['means']['y'].append((np.mean(y_srt),k)) 
						extrema['x']['hi'].append(np.percentile(x_srt,GLOBAL_HI))
						extrema['y']['hi'].append(np.percentile(y_srt,GLOBAL_HI))
						extrema['x']['lo'].append(np.percentile(x_srt,100-GLOBAL_HI))
						extrema['y']['lo'].append(np.percentile(y_srt,100-GLOBAL_HI))
						extrema[k]['x']= (np.percentile(x_srt,100-LOCAL_HI),np.percentile(x_srt,LOCAL_HI))
						extrema[k]['y']= (np.percentile(y_srt,100-LOCAL_HI),np.percentile(y_srt,LOCAL_HI))
						 
					
					xMax,xMin,yMax,yMin = max(extrema['x']['hi']),min(extrema['x']['lo']), max(extrema['y']['hi']),min(extrema['y']['lo'])

					xGMIN,xGMAX =sorted(extrema['means']['x'])[0][1],sorted(extrema['means']['x'])[-1][1]
					yGMIN,yGMAX = sorted(extrema['means']['y'])[0][1],sorted(extrema['means']['y'])[-1][1]

					for name in outlier_key['samples']: 
						
						nx,ny,ns = outlier_key['samples'][name] 
						xE = extrema[ns]['x']						
						yE = extrema[ns]['y']
						ANNO=0



						if   xGMAX != ns and nx > xMax and nx > xE[1]: ANNO=1
						elif xGMIN != ns and nx < xMin and nx < xE[0]: ANNO=2
						elif yGMAX != ns and ny > yMax and ny > yE[1]: ANNO=3 
						elif yGMIN != ns and ny < yMin and ny < yE[0]: ANNO=4
						else: 			       continue 



						self.axes[-1].text(nx,ny,name,color='black',fontsize=7,alpha=0.7)
				
				self.axes[-1].set_xlabel(axes_labels[x_comp],fontweight='bold'); self.axes[-1].set_ylabel(axes_labels[y_comp],fontweight='bold')

				self.yLoc += 1
				if self.yLoc == self.yLen:
					self.xLoc +=1
					self.yLoc = 0 

				break 
#		if self.options.title:	plt.suptitle(self.options.title,fontsize=20,fontweight='bold')

		AXES, LOC, fs, nc = len(self.axes), 'upper center', 8, 1+len(labels)/2 
		if self.yLen == 1:     leg = self.axes[0].legend(items,labels,handletextpad=-0.1,columnspacing=0.1,ncol= nc ,  fontsize=fs,bbox_to_anchor=(0.5,1.30),loc=LOC)
		elif self.yLen == 2:   leg = self.axes[0].legend(items,labels,handletextpad=-0.1,columnspacing=0.1,ncol= nc ,  fontsize=fs,bbox_to_anchor=(1.1,1.30),loc=LOC)
		elif self.yLen == 3:   leg = self.axes[1].legend(items,labels,handletextpad=-0.1,columnspacing=0.1,ncol= nc ,  fontsize=fs,bbox_to_anchor=(0.5,1.30),loc=LOC)
		else: 		       leg = self.axes[1].legend(items,labels,handletextpad=-0.1,columnspacing=0.1,ncol= nc ,  fontsize=fs,bbox_to_anchor=(0.5,1.30),loc=LOC)
		


		return self

	def finish(self,outname=None):

		plt.subplots_adjust(top=0.825)
 
		if outname != None: self.fig.savefig(outname+'.png',dpi=300) 
		#if outname != None: self.fig.savefig(outname+'.pdf',dpi=300) 
		if self.options.show: plt.show()

		self.axes = [] 
		plt.clf() 
		return self




		
class iterplot:
        def __init__(self,samples,features,options=None,progress=None,ITER=4):





      		seaborn.set(rc={'axes.facecolor': options.plotColors['AXIS'],  'figure.facecolor':options.plotColors['FACE']})


	
                self.fig = matplotlib.pyplot.gcf()
                self.fig.patch.set_facecolor(options.plotColors['FACE']) 
               	matplotlib.rcParams['savefig.facecolor'] = options.plotColors['FACE']
                matplotlib.rcParams['ytick.labelsize'] = 7.5



                self.fig.set_size_inches(18.5, 9.5)
		self.xLen, self.yLen = ITER,2
		self.xLoc, self.yLoc = 0,0 
		self.options = options
		self.axes = [] 

		self.features = features 
		self.samples = samples 
		self.sample_key = {s.name: s for s in samples} 





	def add_data(self,iter_dim,SAMPLES,DIMS=3,P_LEN=100):


		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		pts, coefs, axes_labels = iter_dim['pts'], iter_dim['coefs'], iter_dim['axes']



		for p,s in zip(pts,SAMPLES):

			samples = [self.sample_key[sx] for sx in s.split('@')]
		
			if len(samples) == 1 and samples[0].label != None: 
				self.ax.plot(p[0],p[1], alpha=0.8, **samples[0].label.vals) 
			elif len(samples) > 1:
				labels = [s.label for s in samples] 

				colors = list(set([label.vals['color'] for label in labels]))
				msize  = sum([label.vals['markersize'] for label in labels])/1.5

				if len(colors) == 1:
					self.ax.plot(p[0],p[1], alpha=0.8, color = colors[0], markersize = msize)
					
				else:
					c1,c2 = colors[0],colors[1] 
					self.ax.plot(p[0],p[1], alpha=0.8, marker='o',markersize=msize,c=c1, markerfacecoloralt=c2, fillstyle='left')

		self.yLoc +=1	
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		sc = MinMaxScaler(feature_range=(1,100))
                #sY = [[x[0] for x in sc.fit_transform(np.array(y).reshape(-1,1))] for y in Y]


		for s in self.samples: 

			s_vals = sorted(s.cnts.items(), key=lambda X: X[1],reverse=True)[0:P_LEN] 
			my_vals =  sc.fit_transform(np.array([sv[1] for sv in s_vals]).reshape(-1,1))

			s_projections = [[] for j in range(DIMS)] 
			for i,(fi,fv) in enumerate(s_vals): 
				#for j in range(DIMS): s_projections[j].append(fv * coefs[j][fi][1])
				for j in range(DIMS): s_projections[j].append(my_vals[i] * coefs[j][fi][1])

			s_means = [np.mean(x) for x in s_projections]
			if s.label != None: 
				self.ax.plot(s_means[0],s_means[1], alpha=0.8, **s.label.vals) 
			else:
				print 'wtf'
				sys.exit() 

		self.xLoc +=1 
		self.yLoc = 0 





		

		

























			


		

		














		

























			


		

		













