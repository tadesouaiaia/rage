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
import math
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

from sklearn.preprocessing import MinMaxScaler






























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





def sim_uni(a,b,c,x=1):
	return [np.random.uniform(a,b,x) for i in range(c)]
def sim_norm(a,b,c,x=1):
	return [np.random.normal(a,b,x) for i in range(c)]



class K_Means:
        #def __init__(self,samples,colors,ids,bw=0.2):
        def __init__(self,data):

		self.tm = 0.80
		scaler = MinMaxScaler()


		if np.isscalar(data[0]): 
			self.type = '1d'
			self.vals = scaler.fit_transform(np.array(data).reshape(-1,1)) 
		
		else:
			if np.isscalar(data[0][0]):
				if len(data[0]) == 2: 
					self.type = '2d'
					self.vals = scaler.fit_transform(data) 
				elif len(data[0]) == 1:
					self.type = '1d'
					self.vals = scaler.fit_transform(data)
				else:
					print 'here is where i be' 
					sys.exit()
			
			else:
				print 'whoa weird',data[0]
				print data
				sys.exit() 

			 
		if self.type == '1d': self.standard = 8 
		else:		      self.standard = 2.5


	def test(self,test_num=10):
		
		self.km = [KMeans(n_clusters=p) for p in range(1,test_num)]

		self.test_ks()
		i_scrs = [sum([x[0]*x[1] for x in R[-1]])/self.samples for R in self.result.values()]
		g_scrs = [(R[0]+R[1]+R[2])/3.0 for R in self.result.values()]

		return np.mean(i_scrs), np.mean(g_scrs) 


	def test_ks(self):
		self.result = {} 
		self.standard = 4
		self.run = [self.km[p].fit(self.vals) for p in range(1,len(self.km))]	
		k_centers, k_labels = [k.cluster_centers_ for k in self.run], [k.labels_ for k in self.run]
		self.samples = len(k_labels[0]) 
 
		for i,centers,labels in zip(range(len(self.run)),k_centers,k_labels):
			center_key, val_key = dd(list), dd(lambda: dd(list))
			for k,c in enumerate(centers):
				center_key[k] = sorted([(dist(centers[m],c),m) for m in range(len(centers)) if m != k])
				for j in range(len(self.vals)):
					val_key[labels[j]][k].append(dist(self.vals[j],c))
			if i < 0:	continue  
			else:
				center_stats = {} 
				for myC,altCs in sorted(center_key.items()):
					mDists, nAvgs, nScrs  = val_key[myC][myC], [], []
					mArray = [[mD,[]] for mD in mDists]
					mAvg = np.mean(sorted(mDists)[0:int(1.0+len(mDists)*self.tm)])
					for (jScr,jNum) in altCs:
						jDists = val_key[myC][jNum]
						for n in range(len(jDists)): mArray[n][-1].append(jDists[n])
						jScrs   = sorted([pJ/(pJ+pM) for pM,pJ in zip(mDists,jDists)],reverse=True)[0:int(1.0+len(mDists)*self.tm)]
						jAvg = np.mean(sorted(jDists[0:int(len(mDists)*self.tm)]))
						if len(jScrs) < 1: nScrs.append(0.0)
						else:		   nScrs.append(sum(jScrs)/float(len(jScrs)))
						nAvgs.append(jAvg)  	
					center_stats[myC] = [mAvg,[nAvgs,nScrs],mArray]
				center_sizes = cc(labels).values() 
				c_min,c_avg,c_max,n_clusts = [], [], [],[]  
				for cent in sorted(center_stats):
					mAvg,[cAvgs,cScrs],cArray = center_stats[cent]
					c_dists = [round(cAvg/(mAvg+cAvg),3) for cAvg in cAvgs]
					c_min.append(min(c_dists))
					c_avg.append(round(np.mean(c_dists),3))
					c_max.append(max(c_dists))	
					n_rates = sorted([min([d/(m+0.0000000001) for d in D]) for m,D in cArray])
					if len(n_rates) == 0:   n_scr = 0.0 
					else: 			n_scr = len([n for n in n_rates if n>self.standard])/float(len(n_rates))
					n_clusts.append([len(n_rates),round(n_scr,3)])
				self.result[i+2] = [round(sum(c_avg)/len(c_avg),3),round(sum(c_min)/len(c_min),3),round(sum(c_max)/len(c_max),3),n_clusts]
	
					
				
 
def dist(d1,d2):
	
	return np.linalg.norm(d1-d2)*np.linalg.norm(d1-d2)
	print sum([(d1[i]-d2[i])*(d1[i]-d2[i]) for i in range(len(d1))])**0.5
	print np.linalg.norm(d1-d2)


def k_test():
	km = [KMeans(n_clusters=p) for p in range(1,10)]


	
	
	
	vals = [[5],[5],[3],[3],[1],[1],[9],[9]]
	
	u1 = sim_uni(10,20,20) 
	for a,b,c in [[50,60,20],[90,110,20],[130,140,20]]:
		u1+= sim_uni(a,b,c)

	n1 = sim_norm(40,20,20)+sim_norm(60,20,20)+sim_norm(80,20,20)+sim_norm(250,25,20)+sim_norm(400,25,20)  
#	for a,b,c in [[50,10,20],[90,10,20],[130,10,20]]:
#		n1+= sim_norm(a,b,c)
	
	
	k_means = K_Means([n[0]*1.0 for n in n1]) 
	k_means.test()
	return 
	sys.exit() 
	### 
 
	scaler = MinMaxScaler()
	vals = u1
	vals = n1

	vals = scaler.fit_transform(vals) 

	kmeans = [km[p].fit(vals) for p in range(len(km))]
	k_centers = [k.cluster_centers_ for k in kmeans]
	k_scores  = [-1*k.score(vals) for k in kmeans]

	rel_rates = [(k_scores[0]-v) / k_scores[0] for v in k_scores]
	rel_divs = [k_scores[i-1]/k_scores[i] for i in range(1,len(k_scores))]
	rel_divs = [rel_divs[i-1]/rel_divs[i] for i in range(1,len(rel_divs))]

	print k_scores
	print k_centers
	print rel_rates
	print rel_divs
	print rel_divs.index(max(rel_divs))+2








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


def kde1D(vals, bw):
	m_vals = np.array(vals).reshape(-1,1) 
	kde = KernelDensity(kernel='linear', bandwidth=bw).fit(m_vals)
	kde = KernelDensity(kernel='tophat', bandwidth=bw).fit(m_vals)
	kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(m_vals)
	minV,maxV = min(vals),max(vals)
	X_plot = np.linspace(int(minV)-1,int(maxV)+1, 1000)[:, np.newaxis]
	m_scr = np.exp(kde.score_samples(X_plot))	
	x_pts = X_plot[:,0] 
	maxM,minM, modes = float(max(m_scr)), float(min(m_scr)),[] 
	return X_plot[:,0],m_scr
	

def find_local_maxima(x,y,p1=0.01,p2=0.05,bl=75,wl=50,sl=25,jl=10):

	total_cands, local_cands , scored_result, obsMax, endPass, midPass = [], [], [], [0,0], 0, 0
	eMax,eMean,eMed,eLoc1,eLoc2 = 0.10,  0.50, 0.50, 0.5, 0.5 
	fMax,fMean,fMed,fLoc1,fLoc2 = 0.05,  0.30, 0.30, 0.5,  0.5
	xy = [(x[i],y[i]) for i in range(len(x))]
	yRel =  [xy[i][1] for i in range(len(xy)) if xy[i][0] > 0.1 and xy[i][0] < 0.9]
	gMax,maxY,meanY,medY = max(y),max(yRel), np.mean(yRel), np.median(yRel) 
	a = np.array(y)
	max_pts = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
	maxima = [(i,x[i],y[i]) for i in range(wl,len(max_pts)-wl) if max_pts[i]]
	max_pass = [m for m in maxima if m[2]/maxY > p1]
	try: 
		max_groups = [[max_pass[0]]]
	except IndexError: 
		return scored_result,[0,0,0,0,0],0.5
	for i in range(1,len(max_pass)):
		if max_pass[i][0] - max_groups[-1][-1][0] < wl: max_groups[-1].append(max_pass[i])
		else:						max_groups.append([max_pass[i]])
	for mg in max_groups:
		memP =  sorted(list(set([m for m in mg if m[1] < 0.02 or m[1] > 0.98] + [sorted([(m[-1],m) for m in mg])[-1][-1]])))
		for maxG in memP: 
			maxV,mL = maxG[-1] , maxG[0] 
			if maxG[1] > 0.02 and maxG[2] > obsMax[1]: obsMax = [maxG[1],maxG[2]]
			locMin = max(min(y[mL-wl:mL]),min(y[mL:mL+wl]))
			locAvg = sum(y[mL-sl:mL-jl]+y[mL+jl:mL+sl]) / (2.0*(sl-jl))
#			maxScr, locScr1, locScr2, meanScr, medScr  = maxV/gMax,  1-(locAvg/maxV), 1-locMin/maxV,    1-(meanY/maxV),  1-(medY/maxV)
			maxScr, locScr1, locScr2, meanScr, medScr  = maxV/gMax,  1-(locAvg/(locAvg+maxV)), 1-(locMin/(locMin+maxV)),    1-(meanY/(maxV+meanY)),  1-(medY/(medY+maxV))
			local_cands.append([maxG[1],maxG[2],[maxScr,meanScr,medScr,locScr1,locScr2]])
	for m in local_cands: 
 		if m[0] < 0.05 and m[2][0] > eMax and m[2][1] > eMean and m[2][2] > eMed and m[2][3] > eLoc1 and m[2][4] > eLoc2: 
			scored_result.append([m[0],m[1],True,round(sum(m[2])/5.0,3)])
			endPass+=1
		elif m[0] > 0.05 and m[2][0] > fMax and m[2][1] > fMean and m[2][2] > fMed and m[2][3] > fLoc1 and m[2][4] > fLoc2: 
			scored_result.append([m[0],m[1],True,round(sum(m[2])/5.0,3)])
			midPass+=1 		
		else:
			scored_result.append([m[0],m[1],False,round(sum(m[2])/5.0,3)])
	return scored_result, [len(maxima),len(max_pass),len(scored_result),endPass+midPass,midPass],obsMax[0]


def scale_data(data):
	scaler = MinMaxScaler()
	d = np.array(data,dtype=float)
	return scaler.fit_transform(d.astype(float).reshape(d.shape[1],-1)).reshape(d.shape[0],-1) 

def scale_vals(vals):
	scaler = MinMaxScaler() 
	return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0] 


def kscore_pts(vals,r=20): 	
	km = [KMeans(n_clusters=p) for p in range(1,r)]
	my_vals = np.array(vals).reshape(-1,1)
	scores = [-1*km[p].fit(my_vals).score(my_vals) for p in range(len(km))]
	rel_rates = [(scores[0]-v) / scores[0] for v in scores]
	return scores,rel_rates 

def kscore_comps(vals,r=20): 	

	km = [KMeans(n_clusters=p) for p in range(1,r)]
	scores = [-1*km[p].fit(vals).score(vals) for p in range(len(km))]
	rel_rates = [(scores[0]-v) / scores[0] for v in scores]
	return scores,rel_rates




def val_analyze(f_data):
	feats,vals, fname  = f_data.feats, f_data.vals, f_data.f_name 
	xLen, yLen, xLoc,yLoc = 10,10,0,0
	axs = [plt.subplot2grid((xLen,yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)]
	axs = []
	f_vars, f_modes, f_obs, rel_8  = {}, {}, {} , [] 
	f_top,f_tops,rel_dict,scr_dict,s_vals, mix_vals, kMins = {},dd(lambda: {}),{}, {}, {}, dd(list),dd(int)  
	scaled_vals = scale_data(vals) 

	k_dict = {} 
	for i in range(len(vals)):
		f1,clr, lw = feats[i],'purple', 0.9 
		f_obs[f1] = len([v for v in vals[i] if v > 0])/float(len(vals[i]))
		s_vals[f1] = scaled_vals[i] 
		f_vars[f1] =  [stats.variation(scaled_vals[i]),stats.variation(vals[i]) ,stats.variation(scaled_vals[i])/stats.variation(vals[i]) ]
		scr_dict[feats[i]], rel_dict[feats[i]] = kscore_pts(vals[i])	
		k_dict[feats[i]] = K_Means(scaled_vals[i]).test() 

	relMax = min([0.99, max([rd[8] for rd in rel_dict.values()])])
	for f in feats: 
		f1,clr, lw = feats[i],'purple', 0.9 
		if f[0] != 'r': clr = 'cyan' 	
		kMins[f] = [r for r in range(len(rel_dict[f])) if rel_dict[f][r] > relMax or r > 18][0] 
		#axs[-1].plot(range(len(rel_dict[f])), rel_dict[f], color = clr,linewidth=lw)
	kMin,kMax = min(kMins.values()),max(kMins.values())
	f_infos = {f: (rel_dict[f][kMax] - rel_dict[f][1]/float(kMax)) for f in feats}
	xLoc,yLoc, r_cnt, f_cnt = 0,0,0,0 

	f_iter = 1 
	for i in range(len(vals)):	
		x,y = kde1D(scaled_vals[i],0.05) 	
		modes,mode_cnts,mode_shift = find_local_maxima(x,y) 
		f_modes[feats[i]] = [modes,mode_cnts,mode_shift]
		#if (mode_cnts[-1] == 0 and mode_cnts[-2] < 2 and mode_cnts[-3] <4): continue 
		if (mode_cnts[-1] == 0 and mode_cnts[-2] < 2 and mode_cnts[-3] <5): continue 
		if f_cnt > 99:
			for ax in axs: 
				ax.set_xticks([]) 
				ax.set_yticks([]) 
			plt.suptitle(fname+'-'+str(f_iter)) 
			plt.show()
			f_iter+=1
			#plt.clf() 
			axs = [] 
			xLoc,yLoc, r_cnt, f_cnt = 0,0,0,0 
 		if f_cnt <= 99:
			axs.append(plt.subplot2grid((xLen,yLen), (xLoc,yLoc), rowspan = 1, colspan = 1))
			axs[-1].plot(x,y,color='orange',zorder=1)
			for (a,b,VERDICT,scr) in modes:
				if VERDICT: mClr = 'lime' 
				else:	    mClr = 'red' 	
				axs[-1].scatter(a,b,c=mClr,marker='*',s=int(scr*250),zorder=2,edgecolor=mClr) 
			axs[-1].set_title(feats[i].split(';')[-1],fontweight='bold') 
			axs[-1].set_xlim(-0.05,1.05)
			f_cnt, yLoc = f_cnt+1, yLoc + 1 
			if feats[i][0] == 'r': r_cnt+=1 
			if yLoc == yLen: xLoc,yLoc = xLoc+1,0

	for ax in axs: 
		ax.set_xticks([]) 
		ax.set_yticks([]) 
	plt.suptitle(fname+'-'+str(f_iter)) 
	plt.show()

	for i in range(len(vals)):
		f1 = feats[i] 
		f_scrs = scr_dict[f1][0:kMins[f1]+2]
		f_check =  sorted([[(f_scrs[s-1]-f_scrs[s])/(f_scrs[s]-f_scrs[s+1]),s+1]  for s in range(1,len(f_scrs)-1)],reverse=True)
		f_top[f1] = sum(rel_dict[f1][kMin:kMax+1]) / (kMax-kMin+1.0) 
		for k,(a,b) in enumerate(f_check): f_tops[f1][b] = [k+1,b,a]	



	for a,B in f_modes.items():
		print a,'mode-offset',int(100*round(B[-1],3)),'mode-cnts'," ".join([str(s) for s in B[1]]),'|',
		for b in B[0]: print int(100*round(b[0],3)),round(b[1],3),b[2],'|',
		print ""

	for a,B in f_tops.items():
		bX=sorted(B.values())[0:5]
		print a,'k-top-groups',
		for b in bX: print b[1],b[2],
		print ""  
 
	for a,b in f_top.items():
		print a,'ftop_score',b 

	for a,b in f_vars.items():
		print a,'vars',b[0],b[1],b[2] 		
	for a,b in f_infos.items():
		print a,'info',b,'obs',f_obs[a] 
	
	for f in k_dict.keys():
		print f,'ind_grp_kscores',k_dict[f][0],k_dict[f][1]
	for ax in axs: 
		ax.set_xticks([]) 
		ax.set_yticks([]) 
	plt.suptitle(fname) 
	plt.show()
	sys.exit() 



def plot_kde_w_modes(x,y,ax1,j,jLen):	

	color_key = {True: 'lime', False: 'red'}
	modes,mode_cnts,mode_shift = find_local_maxima(x,y) 
	ax1.plot(x,y,color='orange',zorder=1)
	for (a,b,VERDICT,scr) in modes:	ax1.scatter(a,b,c=color_key[VERDICT],marker='*',s=int(scr*250),zorder=2,edgecolor=color_key[VERDICT]) 
	ax1.set_title('pca_iter_'+str(j)+' ( '+str(jLen)+' genes, '+str(len(modes))+' modes)')
	ax1.set_xlim(-1,2)
	ax1.set_xticks([])
	ax1.set_yticks([])
	return modes,mode_cnts,mode_shift 


def plot_paired_pca_pts(ax1,ax2,pts1,pts2,samples,ids,clrs):


	for n in range(len(pts1)):
		if ids[samples[n]] == 'NA': continue 
		try: 				ax1.scatter(pts1[n][0],pts1[n][1],c = clrs[samples[n]])
		except IndexError:		ax1.scatter(pts1[n][0],pts1[n][0],c = clrs[samples[n]])
		try: 				ax2.scatter(pts2[n][0],pts2[n][1],c = clrs[samples[n]])
		except IndexError:		ax2.scatter(pts2[n][0],pts2[n][0],c = clrs[samples[n]])
	
	ax1.set_title('CurrentPCA')
	ax2.set_title('previousPCA')
	ax1.set_xticks([]); ax1.set_yticks([]); ax2.set_xticks([]); ax2.set_yticks([])
	ax1.set_xlim(ax1.get_xlim()[0]-(0.2*ax1.get_xlim()[0]),ax1.get_xlim()[1]-(ax1.get_xlim()[1]*.1)) 
	ax2.set_xlim(ax2.get_xlim()[0]-(0.1*ax2.get_xlim()[0]),ax2.get_xlim()[1]-(ax2.get_xlim()[1]*.1)) 
	return None 


class pca_analyzer:
        #def __init__(self,samples,colors,ids,bw=0.2):
        def __init__(self,data,bw=0.2):
		self.data = data 
		self.bw = bw 
		self.samples, self.colors, self.ids = data.samples, data.sample_color, data.key 
		self.xLen, self.yLen, self.xLoc,self.yLoc = 12,6,0,0
		self.iter = 0 
		self.prev_pts = None 

	def update_prev_pts(self,vPts):
		self.prev_pts = vPts 


	def analyze_component(self,pca_arrays,pca_vals,pca_feats,pca_coeffs,PRINT=True):

		ax1 = plt.subplot2grid((self.xLen,self.yLen),(self.xLoc,self.yLoc))
		ax2 = plt.subplot2grid((self.xLen,self.yLen),(self.xLoc,self.yLoc+1))
		ax3 = plt.subplot2grid((self.xLen,self.yLen),(self.xLoc,self.yLoc+2))
		self.iter += 1
		pca_pts = [sum(p) for p in pca_arrays]
		pca_scale = scale_vals(pca_pts)
		x,y = kde1D(pca_scale,self.bw) 
		self.k_run, self.kPts, self.k_coefs, self.k_vars = run_pca(pca_vals,pca_feats,0) 
		self.kde_modes 		   	= plot_kde_w_modes(x,y,ax1,self.iter,len(pca_feats))
		self.pair_compare 		= plot_paired_pca_pts(ax2,ax3,self.kPts,self.prev_pts,self.samples,self.ids,self.colors)
		if self.iter == 1: plt.legend(self.data.items,self.data.labels,loc = 'upper right',ncol = self.data.ncol,bbox_to_anchor=(1.5,1.8))
			
		if PRINT: self.print_iter_results(pca_feats,pca_pts,pca_vals)

		if self.xLoc + 1   == self.xLen and self.yLoc == 3:  return False 
		elif self.xLoc + 1 == self.xLen: 		     self.xLoc,self.yLoc = 0,3 
		else:						     self.xLoc +=1 

		return True 



	def print_iter_results(self,pca_feats,pca_pts,pca_vals):
		

		modes,mode_cnts,mode_shift = self.kde_modes 
		result = {'pca': K_Means(pca_pts).test(), 'current': K_Means(self.kPts).test(), 'previous': K_Means(self.prev_pts).test(),'modes': [modes,mode_cnts,mode_shift]}
		print 'pcaIter',self.iter,len(pca_feats),'|','pca/k/prev-ig_scrs',
		print round(result['pca'][0],4),round(result['pca'][1],4),
		print round(result['current'][0],3),round(result['current'][1],3),
		print round(result['previous'][0],3),round(result['previous'][1],3),'|',
		print 'mode-offset',int(100*round(mode_shift,3)),'mode-cnts'," ".join([str(s) for s in mode_cnts]),'|',
		for b in modes: print int(100*round(b[0],2)),round(b[1],2),b[2],'|',
		print " ".join(pca_feats)
		return result







def pca_interate(data): 
	feats,vals,samples, sample_clr, id_key = data.feats, data.vals, data.samples,data.sample_color, data.key 
	xLen, yLen, xLoc,yLoc = 12,6,0,0
	my_vals, my_feats, v_out = scale_data(vals), [f for f in feats], [] 
	v_run, vPts, v_coefs, v_vars = run_pca(my_vals,my_feats) 

	pca_analysis = pca_analyzer(data)
	pca_analysis.update_prev_pts(vPts) 
 	
	for j in range(30):
		v_sum,v_last,v_kept = 0,0,{}
		for i,vc in enumerate(v_coefs):
			v_sum += vc[0] 
			v_kept[vc[-1]] = vc[1]
			v_out.append(vc[-1]) 
			if v_sum > 0.90: break 
			elif i > 0 and vc[0] / v_coefs[0][0] < 0.05: break 
			elif i > 0 and vc[0] / v_last < 0.2: break  
			v_last = vc[0] 
		jLen, pca_pts, pca_vals, pca_feats, new_vals, new_feats = len(v_kept.keys()), [[] for i in range(len(vPts))], [], [] , [] , [] 
		for (f,v) in zip(feats,vals):
			if f in v_kept:
				pca_vals.append(v);     pca_feats.append(f) 
				for k in range(len(v)):	pca_pts[k].append(v[k]*v_kept[f])	
			elif f not in v_out:
				new_vals.append(v); new_feats.append(f) 

		status = pca_analysis.analyze_component(pca_pts,pca_vals,pca_feats,v_kept)
		if not status or len(new_vals) == 0: break 

		my_vals, my_feats = [a for a in new_vals], [a for a in new_feats] 
		v_run, vPts, v_coefs, v_vars = run_pca(my_vals,my_feats) 
		pca_analysis.update_prev_pts(vPts) 

	plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.90,wspace=0.05,hspace=0.50)
	plt.suptitle(data.f_name) 
	plt.show() 



def run_pca(vals,my_feats,comps=2):

	fitMat = np.matrix(vals)
   	fitMat = fitMat - fitMat.mean(axis=0) 

	if len(vals) < 3 and comps == 0: 	return False, [[vals[0][n],vals[-1][n]] for n in range(len(vals[0]))], False , [1,1]
	elif len(vals) > 2:			v_run = PCA(n_components = max(comps,2)).fit(fitMat.getT())
	else:					v_run = PCA(n_components = 1).fit(fitMat.getT())
	
       	vPts = v_run.transform(fitMat.getT())
	v_coefs = sorted([(c*c,c,f) for (c,f) in zip(v_run.components_[0],my_feats)],reverse=True) # for comps in v_run.components_]
	var_rates = v_run.explained_variance_ratio_ 
	
	return v_run, vPts, v_coefs, var_rates+[var_rates[0]]


def pca_analyze(data):

	pca_iter_result = pca_interate(data) 
	#pca_iterate(vals,feats,data.samples,data.sample_color) 

	sys.exit() 
	pcaC = len(vals) - 1 
	xLen, yLen, xLoc,yLoc = 2,2,0,0
	ax1 = plt.subplot2grid((xLen,yLen), (0,0), rowspan = 1, colspan = 1) 
	ax2 = plt.subplot2grid((xLen,yLen), (0,1), rowspan = 1, colspan = 1) 
	fitMat = np.matrix(vals)
   	fitMat = fitMat - fitMat.mean(axis=0) 
	scaleMat = np.matrix(scaled_vals) 
	scaleMat = scaleMat - scaleMat.mean(axis=0) 

       	v_run = PCA(n_components = pcaC).fit(fitMat.getT())
       	vPts = v_run.transform(fitMat.getT())	
	v_coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,feats)],reverse=True) for comps in v_run.components_]

       	s_run = PCA(n_components = pcaC).fit(scaleMat.getT())
       	sPts = v_run.transform(scaleMat.getT())	
	s_coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,feats)],reverse=True) for comps in s_run.components_]
	
	s_co = sys.stdout
	
	v_cos,s_cos = dd(lambda: {}), dd(lambda: {}) 
	for cz,c in enumerate(v_coefs):
		v_cos[cz+1] = {x[-1] : [x[0],x[1],xz+1] for xz,x in enumerate(c)}
	for cz,c in enumerate(s_coefs):
		s_cos[cz+1] = {x[-1] : [x[0],x[1],xz+1] for xz,x in enumerate(c)}

	
	for i in range(1,10):
		for f in v_cos[i]:
			print f,i,'v',v_cos[i][f][0],v_cos[i][f][1],v_cos[i][f][2],'s',s_cos[i][f][0],s_cos[i][f][1],s_cos[i][f][2]
	
	sys.exit() 
	for i in range(len(vPts)):
		sample = samples[i]
		v1,v2 = vPts[i][0:2]
		s1,s2 = sPts[i][0:2]
	
		ax1.scatter(v1,v2, c= data.sample_color[sample])
		ax2.scatter(s1,s2, c = data.sample_color[sample]) 

	plt.show()


	sys.exit() 


#def pair_analyze()



class test_data:
        def __init__(self,f_data,f_key,options):

		self.options = options 
		self.colors = ['darkblue','g','brown','r','purple','gray','orange','k','cornflowerblue','magenta','cyan','lime','yellow','pink','blue']
		self.tFS = 25
		self.log = True 
		self.color_key = {'PANTHER': 'darkblue','GORILLA': 'green','ELEPHANT': 'brown','HUMAN':'red','POLAR': 'purple','DOLPHIN': 'gray','TIGER':'orange','ORCA': 'k'}
		self.feats = [] 
		self.vals = [] 
		self.log_vals = []
		self.f_name = f_data 
		for line in open(f_data):
			line = line.split() 
			if line[0] == '---': self.samples = line[1::] 
			else:
				self.feats.append(line[0]) 
				self.vals.append([float(x) for x in line[1::]])
				self.log_vals.append([math.log(v+1.0) for v in self.vals[-1]])
 
		self.key = {} 
		for line in open(f_key): 
			line = line.split() 
			if line[0] == '---': continue 
			self.key[line[0]] = line[1] 

		try: 
			self.sample_color = {k: self.color_key[self.key[k]] for k in self.key.keys()}
		except KeyError:
			self.color_key = {} 
			for i,k in enumerate(list(set(self.key.values()))):
				self.color_key[k] = self.colors[i] 	
			self.sample_color = {k: self.color_key[self.key[k]] for k in self.key.keys()}
		

	
		self.ncol = len(list(set(self.key.values())))
		self.labels,self.items = [],[] 
		for xI,yI in self.color_key.items():
			if xI == 'NA': continue  
			self.labels.append(xI) 
			self.items.append(Rect((0,0),1,1,fc=yI))


	def run_analysis(self,rand=False):
		self.xLen, self.yLen = 2,2 
		ax1 = plt.subplot2grid((self.xLen,self.yLen), (0,0), rowspan = 1, colspan = 1)
		ax2 = plt.subplot2grid((self.xLen,self.yLen), (0,1), rowspan = 1, colspan = 1)
		pcaC,pp,ni = len(self.feats) - 1 , 250, 5000
		if self.log:    vals = self.log_vals
		else:		vals = self.vals 

		if self.options.pca:
			vc = pca_analyze(self)

		else:
			vc = val_analyze(self)

		return 
		sys.exit() 

		fitMat = np.matrix(vals)
   		fitMat = fitMat - fitMat.mean(axis=0) 
       		run = PCA(n_components = pcaC).fit(fitMat.getT())
       		fitPts = run.transform(fitMat.getT())	
		#pca_analyze(self.features,vals,fitPts,run) 


		for i,p in enumerate(fitPts): 
			s,clr = self.samples[i],self.sample_color[self.samples[i]]
			ax1.scatter(p[0],p[1],c=clr) 
		ax1.set_title('PCA') 
		
		s_co = open('my_coeffs.txt','w') 
		coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,self.feats)],reverse=True) for comps in run.components_]

		for cz,c in enumerate(coefs):
			cP,cN = [x for x in c if x[1] > 0 ]  , [x for x  in c if x[1] < 0 ]  
			for xz,x in enumerate(c): 
				s_co.write('%s %d %f %f %d\n' % (x[-1],cz+1,round(x[0],4),round(x[1],4),xz+1))

		tsne = TSNE(n_components=2, verbose=0, perplexity=pp, n_iter=ni)
		ts = tsne.fit_transform(fitPts)

		for i,p in enumerate(ts): 
			s,clr = self.samples[i],self.sample_color[self.samples[i]]
			ax2.scatter(p[0],p[1],c=clr) 

		plt.legend(self.items,self.labels,ncol=self.ncol,bbox_to_anchor=(0.5,1.12),loc='upper right',fontsize=11) 
		ax2.set_title('TSNE') 
		plt.show() 
		sys.exit() 






	def build_pca(self): 

		self.xLen,self.yLen = 4,2
		self.build_comp(self.c1_key) 
		self.build_comp(self.c2_key,'Y') 		
		plt.show() 	






if __name__ == '__main__':

	import sys,os
	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)
	parser.add_option('--title',  dest= "title", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--pca',action='store_true', dest='pca',default=False,help='Show QQ plot')

	(options, args) = parser.parse_args()

	#def sample_creature(w,wv,b,bv,s,sv,c,cv):
#	k_test() 
#	sys.exit() 
#	sys.exit()

	if len(args) == 2: 
		my_data = test_data(args[1],args[0],options) 
		my_data.run_analysis() 
		
