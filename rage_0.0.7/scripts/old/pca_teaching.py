#!/usr/bin/env python

import random
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict as dd
from collections import Counter as cc
import sys
import os
import scipy.stats as stats
from scipy.stats import variation as coVar 

from random import random
import statsmodels.api as sm
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import fdrcorrection as fdr
import random
from math import fabs
#from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
import pickle
from math import log
import numpy as np 
import pylab 
from matplotlib.patches import Rectangle as Rect
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib._png import read_png
from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import Circle as Circ

from sklearn.decomposition import PCA

from sklearn.manifold import TSNE

from random import shuffle
				
from sklearn.cluster import KMeans	


import seaborn
from sklearn.cluster import KMeans

from sklearn.neighbors import KernelDensity
































def get_subplot_size(X): 
	
	
	if X == 1: return 1,1 
	elif X < 4: return 1,X
	elif X == 4: return 2,2 
	elif X < 7:  return 3,2 
	elif X < 10: return 3,3 
	elif X < 13: return 4,3 
	elif X < 17: return 4,4 
	elif X < 21: return 5,4 
	elif X < 26: return 5,5 
	else:        return 6,6 
























def make_size_map(sizes):
	tm,tv = 5,10
	minX,maxX = min(sizes),max(sizes)
	#offsets = [s-minX for x in sizes]
	#print offsets
	return [tm+(s-minX)*(tv/(maxX-minX)) for s in sizes]
	
def make_pv_color_map(vals,scrs,color_key):

	cols = []
	weights = []  
	PV_THRES = 0.05
	PV_2 = 0.005
	if color_key: 
		for val in vals:
			kind,pvX,fcX = scrs[val[1]][0][0]
 			weight = False 
			my_color = []
			wt,cl = 'soft','k'
			if pvX < 0.05:	cl = color_key[kind]
			if pvX < 0.01:  wt = 'bold'


			weights.append(wt) 
			cols.append(cl) 

		return weights,cols
	
	else:
		for val in vals:
			my_color = [] 
			for (a,b,c) in val:
				if b < PV_THRES and c > 0.3: 
					my_color.append('red') 
				elif b<PV_THRES and c < -0.3: my_color.append('blue') 

			if len(my_color) == 0: cols.append('k') 
			elif len(my_color) == 1: cols.append(my_color[0]) 
			else:
				cols.append('k')

		return weights,cols

#def sample_creature([w,wv,b,bv,s,sv,c,cv]):
def sample_creature(stat_list,mp,do,noise = None):
	w,wv,b,bv,s,sv,c,cv = stat_list
	weight = round(np.random.normal(w,wv),0)
	brain  = round(np.random.normal(b,bv),0) 
	speed  = round(np.random.normal(s,sv),2) 
	color  = round(np.random.normal(c,cv),2) 


	creature_noise = [] 

	spot = random.randrange(4) 



	if noise != None:
		for n,trials in noise:
			if trials[0][0] == 'U':
				tLow,tHi = trials[0][1],trials[0][2]
				for i in range(n):
					npt = np.random.uniform(tLow,tHi)
					#if npt < 0: npt = 0
					creature_noise.append(npt) 
			else:

				for i in range(n): 
					mx,sx = trials[random.randrange(len(trials))]
					npt = np.random.normal(mx,sx) 
					#if npt < 1: npt = 0.0 
					creature_noise.append(npt)

#		noise = noise + noise 
#		creature_noise = [] 
#
#		for n,trials in noise: 
#			for i in range(n):
#				spot+=1 
#				if spot >= len(trials): spot = 0
	#			mx,sx = trials[spot]
	#			npt = np.random.normal(mx,sx) 
	#			if npt < 1: npt = 1 
	#			creature_noise.append(npt) 


#		spot = noise[0] 
#		for n,trials in noise:
#			mx,sx = trials[random.randrange(len(trials))]
#			for i in range(n): 
#				#mx,sx = trials[random.randrange(len(trials))]
#				creature_noise.append(np.random.normal(mx,sx))
				

 

#	sys.exit() 

	if weight < 30: weight = 30 
	if speed > 20: speed = 20 
	if speed < 5: speed = 5
	if speed > 13 and s < 15: speed = 13.0 
	if color < 1: color = 1
	if color > 99: color = 99 
	
	if brain > 1000:
		if color / weight < 0.10: 
			brain *= 1.05
			color *= 1.05


	if mp < 1.0:
		weight += np.random.normal(0,weight*mp) 
		brain  += np.random.normal(0,brain*mp) 
		speed  += np.random.normal(0,speed*mp) 
		color  += np.random.normal(0,color*mp) 

	if do > 0: 		
		if random.random() < do: weight = 0.0 
		elif random.random() < do: brain = 0.0 
 		elif random.random() < do: speed = 0.0 
		elif random.random() < do: color = 0.0 	
	'''
	if weight < 1: weight = 1
	if brain < 1:  brain = 1 
	if speed < 1: speed = 1 
	if color < 1: color =1 
	'''
	creature_data = [weight,brain,speed,color]+creature_noise 

	creature_data = [c if c > 0 else 0 for c in creature_data]
	return creature_data 

#	return creature_noise
#	return [weight,brain,speed,color]+creature_noise  
#	return [weight,brain,speed,color]+creature_noise  








def scale_diffs(diffs):
	minD,maxD, max2 =min(diffs),max(diffs), max(diffs) - min(diffs)  
	maxM = 2.0 / (max2) 
	sDiffs =[ -1.0 + ((d - minD) * maxM) for d in diffs]
	return sDiffs
	


def kde_sklearn(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)


def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs): 
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins, 
                      y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)




def kdeMD(vals, bandwidth, xbins = 100j, ybins = 100j, **kwargs):

	
	#i_vals = [np.array([vals[i][j] for i in range(len(vals))]) for j in range(len(vals[0]))]
	i_vals = [np.array([vals[i][j] for i in range(len(vals))]) for j in range(len(vals[0]))]
	v_stats = [] 
	my_grids = [np.mgrid[v.min():v.max():xbins] for v in i_vals]


	#my_sample = np.vstack([iv.ravel() for iv in i_vals]).T 
	#my_train  = np.vstack(i_vals).T

	my_sample = np.vstack([iv.ravel() for iv in i_vals])
	my_train  = np.vstack(i_vals)


    	kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    	kde_skl.fit(my_train)
    	z = np.exp(kde_skl.score_samples(my_sample))	

	return z  




def pca_analyze(feats,vals,pts,run):


	A=kdeMD(vals[0:7],1.0) 
	B=kdeMD(vals[7::],1.0) 

#	for i in range(len(A)):
#		plt.scatter(A[i],B[i]) 
#	plt.show() 

	pca_vals = [[] for i in range(len(pts[0]))]  
	my_pts = [[pts[j][i] for j in range(len(pts[i]))] for i in range(len(pts[0]))]
	my_means =  [np.mean(m) for m in my_pts]

	xLen, yLen, xLoc,yLoc = 3,2,0,0 
	OPTION_A = True 
	#OPTION_A = False 
	ax1 = plt.subplot2grid((xLen,yLen), (0,0), rowspan = 1, colspan = 1) 
	ax2 = plt.subplot2grid((xLen,yLen), (0,1), rowspan = 1, colspan = 1) 
	ax3 = plt.subplot2grid((xLen,yLen), (1,0), rowspan = 1, colspan = 1) 
	ax4 = plt.subplot2grid((xLen,yLen), (1,1), rowspan = 1, colspan = 1) 
	ax5 = plt.subplot2grid((xLen,yLen), (2,0), rowspan = 1, colspan = 1) 
	ax6 = plt.subplot2grid((xLen,yLen), (2,1), rowspan = 1, colspan = 1) 
	scaled_vals = [scale_diffs(v) for v in vals] 

	f_top,f_tops,rel_dict,scr_dict,s_vals, mix_vals, kMins = {},dd(lambda: {}),{}, {}, {}, dd(list),dd(int)  
	km = [KMeans(n_clusters=p) for p in range(1,30)] 

	for i in range(len(vals)):
		f1,clr, lw = feats[i],'purple', 0.9 
		s_vals[f1] = scaled_vals[i] 
		if OPTION_A: my_vals = np.array(scaled_vals[i]).reshape(-1,1) 
		else:	     my_vals = np.array(vals[i]).reshape(-1,1) 
		scores = [-1*round(km[p].fit(my_vals).score(my_vals),2) for p in range(len(km))]
		rel_rates = [(scores[0]-v) / scores[0] for v in scores]
		scr_dict[f1] = scores 
		rel_dict[f1] = rel_rates 
		if f1[0] != 'r': clr = 'cyan'
		kMins[f1] = [r for r in range(len(rel_rates)) if rel_rates[r] > 0.99 or r > 20][0] 
	#	ax1.plot(range(len(rel_rates)), rel_rates, color = clr,linewidth=lw)

	r_ax = 0 
	for i in range(len(vals)):
		x_grid = np.linspace(-4.5, 3.5, 1000)
		m_vals = np.array(vals[i]).reshape(-1,1)

		kde = KernelDensity(kernel='linear', bandwidth=0.5).fit(m_vals)
		kde = KernelDensity(kernel='tophat', bandwidth=0.2).fit(m_vals)
		kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(m_vals)

		minV,maxV = min(vals[i]),max(vals[i]) 
		X_plot = np.linspace(int(minV-1),int(maxV+1), 3000)[:, np.newaxis]
		m_scr = np.exp(kde.score_samples(X_plot))
		print f1
		wows = 0
		maxM = float(max(m_scr))
		for j in range(1,len(m_scr)):
			if j < maxM * 0.10: continue 
			if m_scr[j] > m_scr[j-1] and m_scr[j] > m_scr[j+1]:
				print 'wow',m_scr[j],j 
				wows+=1
		f1 = feats[i] 
		if f1[0] == 'r' and wows > 3: 
			if r_ax == 0: 
				ax1.plot(X_plot[:,0],m_scr,color='orange')
				ax1.set_title(f1)
				r_ax = 1 
			elif r_ax == 1: 
				ax2.plot(X_plot[:,0],m_scr,color='orange')
				ax2.set_title(f1)
				r_ax = 2
			else:
				continue  
		elif f1 == 'weight': 
			ax3.plot(X_plot[:,0],m_scr,color='orange')
			ax3.set_title(f1)
		elif f1 == 'shade':
			ax4.plot(X_plot[:,0],m_scr,color='orange')
			ax4.set_title(f1)
		elif f1 == 'brain':
			ax5.plot(X_plot[:,0],m_scr,color='orange')
			ax5.set_title(f1)
		elif f1 == 'speed':
			ax6.plot(X_plot[:,0],m_scr,color='orange')
			ax6.set_title(f1)
		else:
			continue 
	
	plt.show()
	sys.exit()

	sys.exit() 

		
	kMin,kMax = min(kMins.values()),max(kMins.values())
	clr, lw = 'red', 0.2 
	for i in range(len(vals)):
		f1 = feats[i] 
		f_scrs = scr_dict[f1][0:kMins[f1]+2]
		f_check =  sorted([[round(((f_scrs[s-1]-f_scrs[s])/(f_scrs[s]-f_scrs[s+1])),3),s+1]  for s in range(1,len(f_scrs)-1)],reverse=True)
		f_top[f1] = sum(rel_dict[f1][kMin:kMax+1]) / (kMax-kMin+1.0) 
		for k,(a,b) in enumerate(f_check): 
			f_tops[f1][b] = [k+1,b,a]	
	cands = sorted([(b,a) for (a,b) in f_top.items()],reverse=True) 

	for i in range(len(cands)): 
		for j in range(i+1,len(cands)):
			scr1,f1 = cands[i]
			scr2,f2 = cands[j]
			diffs = [s_vals[f1][a]-s_vals[f2][a] for a in range(len(s_vals[f1]))]
			if OPTION_A: 	sdiffs = np.array(scale_diffs(diffs)).reshape(-1,1) 
			else:		sdiffs = np.array(diffs).reshape(-1,1) 
			scores = [-1*round(km[p].fit(sdiffs).score(sdiffs),2) for p in range(len(km))]
			rel_rates = [(scores[0]-v) / scores[0] for v in scores]
			for p in range(len(rel_rates)): mix_vals[p].append(rel_rates[p])
			rel_dict[(f1,f2)] = rel_rates
			if f1[0] != 'r' or f2[0]  != 'r': clr = 'blue' 
			if f1[0] != 'r' and f2[0]  != 'r': 
				clr,lw = 'green',1.0
			ax2.plot(range(len(rel_rates)), rel_rates, color=clr,linewidth=lw)
			#ax1.plot(range(len(rel_rates)), rel_rates, color=clr,linewidth=lw)
			r1,r2 = rel_dict[f1],rel_dict[f2]
			t1,t2 = f_tops[f1],f_tops[f2] 
			v1,v2 = sorted(t1.values()),sorted(t2.values()) 

			f_scrs = scores[0:max(kMins[f1],kMins[f2])+2]
			f_check =  sorted([[round(((f_scrs[s-1]-f_scrs[s])/(f_scrs[s]-f_scrs[s+1])),3),s+1]  for s in range(1,len(f_scrs)-1)],reverse=True)
			pair_top = {} 
			for k,(a,b) in enumerate(f_check): 
				pair_top[b] = [k+1,b,a]

		#	k_vals = [v1[0][1],v2[0][1],f_check[0][1],f_check[0][2]

			if v1[0][1] == v2[0][1]: 
				print 'whooooo'	

			if v1[0][1] == 2 and v2[1][1] == 2 and f_check[0][1] == 2: 
				print 'state1'		
				print f_check[0] 
				print f_tops[f1]	
				print f_tops[f2]
				for k in range(1,5):
					min_rate = min(r1[k],r2[k]) 
					print f1,f2,k+1,r1[k],r2[k],rel_rates[k],rel_rates[k]/min_rate
			else:
			
				print 'state_2' 
				print f1,f2,v1[0],v2[0],f_check[0]
				print f1,v1 
				print f2,v2
			



	ax1.set_title('K_MEANS_SCORES') 
	ax2.set_title('K_MEANS_SCORE_DIFFERENCES') 
	plt.suptitle(OPTION_A) 
	plt.show()


class animal_plot:
        def __init__(self,args):

		self.tFS = 25
		self.log = False 
		self.color_key = {'PANTHERS': 'darkblue','GORILLAS': 'green','ELEPHANTS': 'brown','HUMANS':'red','POLAR_BEARS': 'purple'}
		self.coord_key = {'weight': 0, 'brain': 1, 'speed':2,'shade':3} 
		self.label_key = {'weight': 'Body Mass (g)', 'brain': 'Brain Mass (g)', 'speed': 'Top Speed (m/s)','shade': 'Grayscale Shade'} 
		self.color_key['ORCAS'] = 'k' 
		self.color_key['DOLPHINS'] = 'gray'
		self.color_key['TIGERS'] = 'orange'

		self.noise = [] 
		self.noise.append([5, [(100,50),(1000,500)]])
		self.noise.append([5, [('U',0,5000)]])
		self.noise.append([10,[(100,25)]])
		self.noise.append([10,[(1000,250)]])
#		self.noise.append([5,[(25,25)]])
#		self.noise.append([5,[(250,250)]])
#		self.noise.append([5,[(2500,2500)]])
#		self.noise.append([5,[(100,200)]])
#		self.noise.append([5,[(1000,2000)]])
#		self.noise.append([1,[(1000,2000)]])
#		self.noise.append([5,[(1000,1000)]])
#		self.noise.append([10,[(1,250)]])
#		self.noise.append([5,[(50,25)]])
#		self.noise.append([10,[(1,250)]])
	#	self.noise.append([1,[(1000,750)]])
#		self.noise.append([5, [(50,10),(500,100)]])
#		self.noise += [[5,[(20000,10000)]]]
#		self.noise = [[5,[(2000,1000)]]]
#		self.noise = [[5,[(2000,1000)]]]
		self.er,self.dp = 0.05,0.025
#		self.er,self.dp = 0.1,0.000001		

		self.ncol = len(self.color_key) 
		
		self.labels,self.items = [],[] 
		for xI,yI in self.color_key.items(): 
			self.labels.append(xI) 
			self.items.append(Rect((0,0),1,1,fc=yI))
	
		self.stats = {} 	
		self.stats['HUMANS']         =      		[  70,  20, 1400, 30,      8.5,   1.0, 50, 20] 
		self.stats['ELEPHANTS']    =      [ 3000, 100,  500, 20,      7.0, 1.0, 15,  8] 
		self.stats['GORILLAS']      =      		[ 200,  50,  500, 20,    9.75, 0.75, 88, 2] 
		self.stats['PANTHERS']     =      	[  60,  14,  150, 20,     18.5, 1, 90, 2] 
		self.stats['POLAR_BEARS']          =      	[ 400,  75,  480, 20,    10.25, 0.75, 8, 2] 
		self.stats['ORCAS']         =       		[ 10000,  1000,  5000, 20,    0.5, 0.01, 90, 2] 
	 	self.stats['DOLPHINS']       =      		[  250,  25, 1600, 100,      0.5,   0.1, 85, 10] 
		self.stats['TIGERS']     =       		[  300,  50, 250, 20,      15.5,   1, 60, 20] 

	

		self.tXLoc,self.tYLoc = 0.94,1.035

		ds = ['mass(kg)','','brain_size(g)','','Top_Land_Speed(m/s)','', 'Grayscale_Shade(%)'] 
#		for s,S in self.stats.items(): 
#			print s, 
#			for i in range(0,len(S),2):
#				print ds[i],S[i],S[i+1],
#			print "" 
			

#		sys.exit() 

	def sample_data(self,num,randvar=None):

		self.num = num 
		SAMPLES = num 
		animals = {} 
		
		for x in ['HUMANS','ELEPHANTS','GORILLAS','PANTHERS','POLAR_BEARS','ORCAS','DOLPHINS','TIGERS']:

			if randvar: 
				myNum = int(np.random.normal(num,num/8.0))
				if myNum < 10: myNum = 10
				animals[x] = [sample_creature(self.stats[x],self.er,self.dp,self.noise) for i in range(myNum)] 

			else:
				animals[x] = [sample_creature(self.stats[x],self.er,self.dp,self.noise) for i in range(SAMPLES)] 


		self.key = animals 	
		self.labels,self.items = [],[] 

		for xI,yI in self.color_key.items():
			 
			if len(self.key[xI]) < 100: 
				self.labels.append(xI + ' ('+str(len(self.key[xI]))+' )')
			else: 
				self.labels.append(xI + ' ('+str(len(self.key[xI]))+')') 
			self.items.append(Rect((0,0),1,1,fc=yI))
	
		self.merge_data() 

	def merge_data(self): 

		self.samples, self.sample_colors, self.features, self.feature_vals, self.feature_log_vals = [], [], [], [], []  

	
		first_ex =  self.key[self.key.keys()[0]][0] 
		nd = len(self.coord_key.values()) 



		k = 1  
		for j in range(len(self.coord_key.values()),len(first_ex)):
			self.coord_key['r'+str(k)] = j 
			k+=1


		for lab,labI in self.coord_key.items(): 
			f_data, f_log_data = [] , [] 
			f_log_data = [] 
			for animal,A in self.key.items():
				for i,a in enumerate(A):
					pt = a[labI] 
					f_data.append(pt) 
					f_log_data.append(log(pt+0.01)) 
					if labI == 0: 
						self.samples.append(animal.strip('S')+'_'+str(i+1))
						self.sample_colors.append(self.color_key[animal]) 	

						

			self.features.append(lab) 
			self.feature_vals.append(f_data) 
			self.feature_log_vals.append(f_log_data) 
		


	def run_tsne(self):

		pp,ni = 200, 4000
		tsne = TSNE(n_components=2, verbose=0, perplexity=pp, n_iter=ni)
		ts = tsne.fit_transform(self.pca_pts)		

		ax = plt.subplot2grid((1,1), (0,0), rowspan = 1, colspan = 1)

		for i in range(len(ts)): 
			ax.scatter(ts[i][0],ts[i][1],c=self.sample_colors[i]) 

		plt.legend(self.items,self.labels,ncol=self.ncol,bbox_to_anchor=(self.tXLoc,self.tYLoc),loc='upper right',fontsize=12) 

		if self.log: plt.title('TSNE for log feature values',fontsize=self.tFS,fontweight='bold',y=1.05) 
		else:		plt.title('TSNE for raw features',fontsize=self.tFS,y=1.05,fontweight='bold') 

		ax.set_xticks([]) 
		ax.set_yticks([]) 


	def run_multi_transform(self,rand=False):

		runs=[25,50,100,250] #250]

		if rand: 
			runs = [10,25,50,100,250,500] 
			runs=[25,50,100,250] #250]
			runs = [50,100] 

			runs = [100] 

		self.multi_transform = {} 
		self.multi_tsne = {} 
		pp,ni = 200, 4000
		pcaC = len(self.features) - 1 

		s_out = open('easy_data.txt','w') 
#		s_co = open('easy_coeffs.out','w') 
		s_co = sys.stdout
		for r in runs: 
			self.sample_data(r,rand) 
			if self.log:    vals = self.feature_log_vals
			else:		vals = self.feature_vals 


			s_out.write('%s\n' % " ".join(["---"]+self.samples))

			for f,V in zip(self.features,self.feature_vals):
				s_out.write('%s\n' % " ".join([f]+[str(v) for v in V]))
				



			fitMat = np.matrix(vals)
   			fitMat = fitMat - fitMat.mean(axis=0) 
       			run = PCA(n_components = pcaC).fit(fitMat.getT())
       			fitPts = run.transform(fitMat.getT())	
			#pca_analyze(self.features,vals,fitPts,run) 

			coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,self.features)],reverse=True) for comps in run.components_]

			for cz,c in enumerate(coefs):
				cP,cN = [x for x in c if x[1] > 0 ]  , [x for x  in c if x[1] < 0 ]  
				for xz,x in enumerate(c): 
					s_co.write('%s %d %f %f %d\n' % (x[-1],cz+1,round(x[0],4),round(x[1],4),xz+1))

			tsne = TSNE(n_components=2, verbose=0, perplexity=pp, n_iter=ni)
			ts = tsne.fit_transform(fitPts)
			self.multi_transform[r] = [fitPts,ts,self.sample_colors,self.items,self.labels]

		self.xLen, self.yLen = 2,2 
		xLoc,yLoc = 0,0 	

		for r,m in sorted(self.multi_transform.items()): 
			pca,tsne,clrs,items,labs = m 
			ax = plt.subplot2grid((self.xLen,self.yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
			for p,c in zip(pca,clrs): 
				ax.scatter(p[0],p[1],c=c) 

			ax.legend(items,labs,ncol=4,bbox_to_anchor=(0.95,1.12),loc='upper right',fontsize=11) 
			ax.set_xticks([]) 
			ax.set_yticks([]) 
			if yLoc == 0: yLoc +=1 
			else: 
				xLoc += 1
				yLoc = 0
		
			xMin,xMax = ax.get_xlim()
			yMin,yMax = ax.get_ylim()
	
			if len(runs) < 3: ax.text(xMin+xMax*0.1,yMax*1.1,'PCA:',fontsize=20) 


		if len(runs) > 2: 
			plt.suptitle('PCA (noiseRate= '+str(2*self.er)+', dropRate= '+str(self.dp)+' )',fontsize=22,fontweight='bold') 
			plt.show() 
			xLoc,yLoc = 0,0 	

		else: 
			nstr = [] 
			for nx,tx in self.noise:

				if tx[0][0] == 'U': 
					nstr.append(str(nx)+'@U('+str(tx[0][1])+','+str(tx[0][2])+')')
				elif len(tx) == 1: 
					nstr.append(str(nx)+'@N('+str(tx[0][0])+','+str(tx[0][1])+')')
 				else:
					nstr.append(str(nx)+'@BiModal( '+", ".join(['N{'+str(z1)+','+str(z2)+'}' for z1,z2 in tx])+' )')
						


		
			supt = 'noiseR= '+str(2*self.er)+', dropR= '+str(self.dp)+', '+"noiseDims= "+", ".join(nstr) 
			plt.suptitle(supt,fontsize=16,fontweight='bold') 


		for r,m in sorted(self.multi_transform.items()): 
			pca,tsne,clrs,items,labs = m 
			ax = plt.subplot2grid((self.xLen,self.yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
			for p,c in zip(tsne,clrs): 
				ax.scatter(p[0],p[1],c=c) 

			ax.legend(items,labs,ncol=4,bbox_to_anchor=(0.95,1.15),loc='upper right',fontsize=11) 
			ax.set_xticks([]) 
			ax.set_yticks([]) 
			if yLoc == 0: yLoc +=1 
			else: 
				xLoc += 1
				yLoc = 0
			
			xMin,xMax = ax.get_xlim()
			yMin,yMax = ax.get_ylim()
			if len(runs) < 3: ax.text(xMin+xMax*0.10,yMax*1.1,'TSNE:',fontsize=20) 

		if len(runs) > 2: 
			plt.suptitle('TSNE (noiseRate= '+str(2*self.er)+', dropRate= '+str(self.dp)+' )',fontsize=22,fontweight='bold') 


		else:
 			plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.88,wspace=0.05,hspace=0.18)
		plt.show() 





		




	def run_pca(self): 
	
		if self.log:    vals = self.feature_log_vals
		else:		vals = self.feature_vals 

		rvals = self.feature_vals
		
		fitMat = np.matrix(vals)
   		fitMat = fitMat - fitMat.mean(axis=0) 
       		run = PCA(n_components = 4).fit(fitMat.getT())	
		fitPts = run.transform(fitMat.getT())
       		fitPts = PCA(n_components = 4).fit_transform(fitMat.getT())	

		self.pca_pts = fitPts 

		coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,self.features)],reverse=True) for comps in run.components_]
		p1,p2 = run.explained_variance_ratio_[0],run.explained_variance_ratio_[1]
		ax = plt.subplot2grid((1,1), (0,0), rowspan = 1, colspan = 1)
		myx = [] 
		self.coefX = {cx[-1]: cx[1] for cx in coefs[0]} 
		self.coefY = {cx[-1]: cx[1] for cx in coefs[1]} 

	

		for i in range(len(fitPts)): 
			x,y = fitPts[i][0],fitPts[i][1]   
			#ax.scatter(x,y,c=self.sample_colors[i])
			myx.append(x) 
			xSpot,ySpot = 0,0 
			for j,f in enumerate(self.features):
				#print self.coefX[f],vals[j][i] 	
				xSpot += self.coefX[f]*vals[j][i] 
				ySpot += self.coefY[f]*vals[j][i] 

			ax.scatter(xSpot,ySpot,c=self.sample_colors[i]) 
#			if self.samples[i] == 'HUMAN_50': 
#				ax.scatter(xSpot,ySpot,c=self.sample_colors[i],marker='*',s=200) 

			#print x,self.sample_colors[i] 


		cstr,fs,c_key  = '',None, [] 


		c_key = [] 
		for c in coefs[0]: 			
			vs,vl = str(fabs(int(100*c[1]))),c[-1] 
			c_add = [c[-1],c[1]]
			if c[1] < 0: 
				cstr+= ' - '+vs+'*'+vl
				c_key.append(c_add) 
			elif not fs: 
				fs = vs+'*'+vl
				c_first = [c_add] 
			else: 	   
				cstr+= ' + '+vs+'*'+vl
				c_key.append(c_add) 
 
		cstr1 = fs + cstr 
		self.c1_key = c_first + c_key 



		cstr,fs,c_key  = '',None, [] 
		for c in coefs[1]: 			
			vs,vl = str(fabs(int(100*c[1]))),c[-1] 
			c_add = [c[-1],c[1]]
			if c[1] < 0: 
				cstr+= ' - '+vs+'*'+vl
				c_key.append(c_add) 
			elif not fs: 
				fs = vs+'*'+vl
				c_first = [c_add] 
			else: 	   
				cstr+= ' + '+vs+'*'+vl
				c_key.append(c_add) 
		cstr2 = fs + cstr 
		self.c2_key = c_first + c_key 

		ax.set_xlabel('PC1 = '+cstr1 + ' ( '+str(round(p1,2))+' var explained )',fontweight='bold')
		ax.set_ylabel('PC2 = '+cstr2 + ' ( '+str(round(p2,2))+' var explained )',fontweight='bold')

		plt.legend(self.items,self.labels,ncol=self.ncol,bbox_to_anchor=(self.tXLoc,self.tYLoc),loc='upper right',fontsize=12) 
	
		if self.log:
			plt.title('PCA for log feature values',fontsize=self.tFS,fontweight='bold',y=1.05) 
		else:
			plt.title('PCA for raw features',fontsize=self.tFS,y=1.05,fontweight='bold') 

		ax.set_xticks([]) 
		ax.set_yticks([]) 

	def build_pca(self): 

		self.xLen,self.yLen = 4,2
		self.build_comp(self.c1_key) 
		self.build_comp(self.c2_key,'Y') 
		
		plt.show() 	


	def build_comp(self,c_key,AXIS='X'):
		if AXIS == 'X': 	xLoc,yLoc = 0,0 
		else:			xLoc,yLoc = 0,1 
		cstr = '' 
		ax = plt.subplot2grid((self.xLen,self.yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
		for i in range(len(c_key)):
 
			ax = plt.subplot2grid((self.xLen,self.yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
			mk =  c_key[0:i+1] 
 

			if cstr == '':
				if AXIS == 'X':
					cstr = 'X= '+str(100*round(mk[-1][1],3))+'*'+mk[-1][0] 
				else: cstr = 'Y= '+str(100*round(mk[-1][1],3))+'*'+mk[-1][0] 
			else:
				cadd =  str(100*fabs(round(mk[-1][1],3)))+'*'+mk[-1][0] 
				if mk[-1][1] > 0: cstr += ' + '+cadd
				else:		  cstr += ' - '+cadd
			my_data = []

 
			for animal,A in self.key.items(): 
				clr = self.color_key[animal]  
				for a in A:
					myVal = 0 

					for x,y in mk:
							
						if self.log: myVal +=  log(a[self.coord_key[x]]+0.01) * y 
						else: 	     myVal +=  a[self.coord_key[x]] * y 
					my_data.append([myVal,0,clr]) 
			shuffle(my_data) 
			for x,y,clr in my_data:
				#if AXIS == 'X':
				y = np.random.normal(0,0.1) 

				if AXIS == 'X': 
					ax.scatter(x,y,c=clr,alpha=0.6,s=50) 
					#if i == 3: print x,clr
				else:
					ax.scatter(y,x,c=clr,alpha=0.6,s=50) 
			ax.set_xticks([]) 
			ax.set_xticks([]) 
			ax.set_yticks([])

		

			if AXIS == 'X': 
				ax.set_ylim(-1,1)   
				if i == 0: ax.set_title('PC1 Build',fontsize=20,fontweight='bold',y=1.1) 


			else: 		
				ax.set_xlim(-1,1)   
				if i == 0: ax.set_title('PC2 Build',fontsize=20,fontweight='bold',y=1.1) 

		

			ax.set_xlabel(cstr,fontweight='bold',fontsize=18)
			xLoc+=1

		if AXIS == 'X':
			plt.legend(self.items,self.labels,ncol=self.ncol,bbox_to_anchor=(2.04,4.74),loc='upper right',fontsize=12) 


	def draw_2d_data(self):
		self.xLen,self.yLen = 4,4
		self.draw_2d('weight','weight',0,0) 
		self.draw_2d('weight','brain',0,1) 
		self.draw_2d('weight','speed',0,2) 
		self.draw_2d('weight','shade',0,3) 
		self.draw_2d('brain','weight',1,0) 
		self.draw_2d('brain','brain',1,1) 
		self.draw_2d('brain','speed',1,2) 
		self.draw_2d('brain','shade',1,3) 	
		self.draw_2d('speed','weight',2,0) 
		self.draw_2d('speed','brain',2,1) 
		self.draw_2d('speed','speed',2,2) 
		self.draw_2d('speed','color',2,3) 
		self.draw_2d('shade','weight',3,0) 
		self.draw_2d('shade','brain',3,1) 
		self.draw_2d('shade','speed',3,2) 
		self.draw_2d('shade','shade',3,3) 

	def draw_2d_brief(self):
		self.xLen,self.yLen = 2,3
		self.draw_2d('weight','brain',0,0) 
		plt.legend(self.items,self.labels,ncol=self.ncol,bbox_to_anchor=(3.05,1.10),loc='upper right',fontsize=12) 
		self.draw_2d('weight','speed',0,1) 
		self.draw_2d('weight','shade',0,2) 
		self.draw_2d('brain','speed',1,0) 
		self.draw_2d('brain','shade',1,1) 	
		self.draw_2d('speed','shade',1,2) 
		
		if self.log: 		plt.suptitle('Log 2d Values for Animal Data') 
		else: 			plt.suptitle('Raw 2d Values for Animal Data') 

	def draw_1d_data(self):
		self.xLen,self.yLen = 2,2
		plot.draw_1d('weight',0,0) 
		plt.legend(self.items,self.labels,ncol=self.ncol,bbox_to_anchor=(1.97,1.10),loc='upper right',fontsize=12) 
		plot.draw_1d('brain',0,1) 
		plot.draw_1d('speed',1,0) 
		plot.draw_1d('shade',1,1)

		if self.log: 		plt.suptitle('Log Values for Animal Data') 
		else: 			plt.suptitle('Raw Values for Animal Data') 


	def draw_2d(self,xC,yC,xLoc,yLoc):

		ax = plt.subplot2grid((self.xLen,self.yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
		my_data = [] 
		for animal,A in self.key.items():
			clr = self.color_key[animal]  
			for a in A: 
				x,y = a[self.coord_key[xC]],a[self.coord_key[yC]]	  
				if self.log: 
					x,y = log(x+0.01) , log(y+0.01) 
				#ax.scatter(x,y,c=clr) 
				my_data.append([x,y,clr]) 
		shuffle(my_data) 
		for x,y,clr in my_data: 
			ax.scatter(x,y,c=clr,alpha=0.8) 
		ax.set_xticks([]) 
		ax.set_yticks([]) 
		ax.set_xlabel(self.label_key[xC],fontweight='bold') 
		ax.set_ylabel(self.label_key[yC],fontweight='bold')


	def draw_1d(self,xC,xLoc,yLoc):
		ax = plt.subplot2grid((self.xLen,self.yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
		my_data = [] 
		for animal,A in self.key.items():
			clr = self.color_key[animal]  
			for a in A: 
				x = a[self.coord_key[xC]]

				if self.log: x = log(x+0.01) 
				#ax.scatter(x,0,c=clr,alpha=0.5)
				my_data.append([x,np.random.normal(0,0.01),clr]) 

		shuffle(my_data) 
		for x,y,clr in my_data: 
			ax.scatter(x,y,c=clr,alpha=0.7,s=300) 

		ax.set_ylim(-0.07,0.07)  
		ax.set_xticks([]) 
		ax.set_yticks([]) 
		ax.set_xlabel(self.label_key[xC],fontweight='bold') 





if __name__ == '__main__':

	import sys,os
	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)
	parser.add_option('--title',  dest= "title", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--pca',  dest= "pca", type = 'str' , default = None, help = "horizontal data")
	(options, args) = parser.parse_args()

	#def sample_creature(w,wv,b,bv,s,sv,c,cv):
	
	plot = animal_plot(options) 
	plot.sample_data(50) 

	plot.log = True
	'''
	
	#plot.draw_1d_data() 
	#plt.show() 
	
#	plot.log = True 
#	plot.draw_1d_data() 
#	plt.show() 

	
	#plot.log = False 
	#plot.draw_2d_brief() 
	#plt.show()
	plot.log = True 
	plot.draw_2d_brief()
	plt.show() 
	 
	
#	plot.log = False 
#	plot.run_pca() 
#	plt.show()
	plot.log = True
	'''
#	plot.run_multi_transform() 
	plot.run_multi_transform(True) 
#	plot.multi_transform() 
	sys.exit() 
#	sys.exit() 

	plot.run_pca() 
	plt.show()

	plot.run_tsne()
	plt.show() 

#	plot.build_pca() 
	sys.exit() 




	for animal,A in animals.items(): 
		weights = [a[0] for a in A] 
		brains  = [a[1] for a in A] 
		speeds  = [a[2] for a in A]
		shades  = [a[3] for a in A] 
		print animal+':','Cnt',len(A)
		print '     STATS (m/std)',
		print 'WEIGHT (kg):',round(np.mean(weights)/1000.0,1),'('+str(round(np.std(weights)/1000.0,1))+')',
		print 'BRAINSIZE (g):',np.mean(brains),'('+str(round(np.std(brains),1))+')', 
		print 'TOPSPEED (m/s):',np.mean(speeds),'('+str(round(np.std(speeds),1))+')', 
		print 'SHADE (percent):',np.mean(shades),'('+str(round(np.std(colors),1))+')' 
	print "" 















