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
	max_groups = [[max_pass[0]]]
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


def kscore_vals(vals,r=20): 
	
	km = [KMeans(n_clusters=p) for p in range(1,r)]
	my_vals = np.array(vals).reshape(-1,1)
	scores = [-1*km[p].fit(my_vals).score(my_vals) for p in range(len(km))]
	rel_rates = [(scores[0]-v) / scores[0] for v in scores]
	return scores,rel_rates 


def val_analyze(feats,vals):
	xLen, yLen, xLoc,yLoc = 8,8,0,0
	axs = [plt.subplot2grid((xLen,yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)]
	f_vars, f_modes, f_obs, rel_8  = {}, {}, {} , [] 
	f_top,f_tops,rel_dict,scr_dict,s_vals, mix_vals, kMins = {},dd(lambda: {}),{}, {}, {}, dd(list),dd(int)  
	scaled_vals = scale_data(vals) 

	for i in range(len(vals)):
		f1,clr, lw = feats[i],'purple', 0.9 
		f_obs[f1] = len([v for v in vals[i] if v > 0])/float(len(vals[i]))
		s_vals[f1] = scaled_vals[i] 
		f_vars[f1] =  [stats.variation(scaled_vals[i]),stats.variation(vals[i]) ,stats.variation(scaled_vals[i])/stats.variation(vals[i]) ]
		scr_dict[feats[i]], rel_dict[feats[i]] = kscore_vals(vals[i])
	relMax = min([0.99, max([rd[8] for rd in rel_dict.values()])])
	for f in feats: 
		f1,clr, lw = feats[i],'purple', 0.9 
		if f[0] != 'r': clr = 'cyan' 	
		kMins[f] = [r for r in range(len(rel_dict[f])) if rel_dict[f][r] > relMax or r > 18][0] 
		axs[-1].plot(range(len(rel_dict[f])), rel_dict[f], color = clr,linewidth=lw)
	kMin,kMax = min(kMins.values()),max(kMins.values())
	f_infos = {f: (rel_dict[f][kMax] - rel_dict[f][1]/float(kMax)) for f in feats}

	yLoc, r_cnt, f_cnt = 1,0,0 

	for i in range(len(vals)):	
		x,y = kde1D(scaled_vals[i],0.05) 	
		modes,mode_cnts,mode_shift = find_local_maxima(x,y) 
		f_modes[feats[i]] = [modes,mode_cnts,mode_shift]
		if f_cnt > 62 or (mode_cnts[-1] == 0 and mode_cnts[-2] < 2 and mode_cnts[-3] <3): continue 
 		else:
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
	
	for ax in axs: 
		ax.set_xticks([]) 
		ax.set_yticks([]) 
	plt.show()
	sys.exit() 



def pca_interate(data): 
	feats,vals,samples, sample_clr = data.feats, data.log_vals, data.samples,data.sample_color
	xLen, yLen, xLoc,yLoc = 12,6,0,0
	ax = plt.subplot2grid((xLen,yLen), (0,0), rowspan = 1, colspan = 1) 
	pre_comps, v_out, color_key = {},[], {True: 'lime', False: 'red'}
	my_vals, my_feats = scale_data(vals), [f for f in feats]
	v_run, vPts, v_coefs, v_vars = run_pca(my_vals,my_feats) 

	for j in range(30):
		v_sum,v_pos,v_neg,v_last,v_kept = 0,0,0.000001,0,{}
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
				pca_vals.append(v) 	
				pca_feats.append(f) 
				for k in range(len(v)):	pca_pts[k].append(v[k]*v_kept[f])	
			elif f not in v_out:
				new_vals.append(v) 
				new_feats.append(f) 

		ax1,ax2,ax3 = plt.subplot2grid((xLen,yLen), (xLoc,yLoc)), plt.subplot2grid((xLen,yLen), (xLoc,yLoc+1)),plt.subplot2grid((xLen,yLen), (xLoc,yLoc+2)) 
		
		pca_pts = [sum(p) for p in pca_pts]
		scores, rel_rates = kscore_vals(pca_pts)
		pca_scale = scale_vals(pca_pts)
		x,y = kde1D(pca_scale,0.2) 
		modes,mode_cnts,mode_shift = find_local_maxima(x,y) 
		ax1.plot(x,y,color='orange',zorder=1)
		for (a,b,VERDICT,scr) in modes:	ax1.scatter(a,b,c=color_key[VERDICT],marker='*',s=int(scr*250),zorder=2,edgecolor=color_key[VERDICT]) 
		
		k_run, kPts, k_coefs, k_vars = run_pca(pca_vals,pca_feats,0) 
		for n in range(len(vPts)):
			try: 			ax2.scatter(kPts[n][0],kPts[n][1],c = sample_clr[samples[n]])
			except IndexError:	ax2.scatter(kPts[n][0],kPts[n][0],c = sample_clr[samples[n]])
			try: 			ax3.scatter(vPts[n][0],vPts[n][1],c = sample_clr[samples[n]])
			except IndexError:	ax3.scatter(vPts[n][0],vPts[n][0],c = sample_clr[samples[n]])
		
		ax1.set_title('pca_iter_'+str(j+1)+' ( '+str(jLen)+' genes, '+str(len(modes))+' modes)')

		ax2.set_title('CurrentPCA')
		ax3.set_title('previousPCA')
#		ax3.set_title('Previous PCA var2='+str(v_vars[-1])) 
#		ax2.set_title('CurrentPCA var2='+str(k_vars[-1]))
#		ax3.set_title('Previous PCA var2='+str(v_vars[-1])) 
		ax1.set_xticks([]); ax1.set_yticks([]); ax2.set_xticks([]); ax2.set_yticks([]); ax3.set_xticks([]); ax3.set_yticks([]) 
		ax1.set_xlim(-1,2)
		ax2.set_xlim(ax2.get_xlim()[0]-(0.1*ax2.get_xlim()[0]),ax2.get_xlim()[1]-(ax2.get_xlim()[1]*.1)) 
		ax3.set_xlim(ax3.get_xlim()[0]-(0.2*ax3.get_xlim()[0]),ax3.get_xlim()[1]-(ax3.get_xlim()[1]*.1)) 
		
		pre_comps[j] = [v_kept.keys(),[scores,[modes,mode_cnts,mode_shift]]]

		print 'pcaIter',j+1," ".join(v_kept.keys()),'|'," ".join([str(s) for s in scores]),'|',
		print 'mode-offset',int(100*round(mode_shift,3)),'mode-cnts'," ".join([str(s) for s in mode_cnts]),'|',
		for b in modes: print int(100*round(b[0],3)),round(b[1],3),b[2],'|',
		print ""

		if j == 0: plt.legend(data.items,data.labels,loc = 'upper right',ncol = data.ncol,bbox_to_anchor=(1.5,1.8))
		xLoc += 1

		if xLoc == xLen and yLoc == 3:  break
		elif xLoc == xLen: 		xLoc,yLoc = 0,3 
		
		if len(new_vals) < 1: break 
		my_vals, my_feats = [a for a in new_vals], [a for a in new_feats] 

		v_run, vPts, v_coefs, v_vars = run_pca(my_vals,my_feats) 

	plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.90,wspace=0.05,hspace=0.50)
	plt.show() 
	return pre_comps 



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

	feats,vals,samples = data.feats, data.log_vals, data.samples

	pcaC = len(vals) - 1 

	xLen, yLen, xLoc,yLoc = 2,2,0,0
	ax1 = plt.subplot2grid((xLen,yLen), (0,0), rowspan = 1, colspan = 1) 
	ax2 = plt.subplot2grid((xLen,yLen), (0,1), rowspan = 1, colspan = 1) 

	scaler = MinMaxScaler()
	scaler.fit(vals) 
	scaled_vals = scaler.transform(vals)  
	pca_iter_result = pca_interate(data) 
	#pca_iterate(vals,feats,data.samples,data.sample_color) 

	sys.exit() 
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
        def __init__(self,f_data):

		self.colors = ['darkblue','g','brown','r','purple','gray','orange','k','cornflowerblue','magenta','cyan','lime','yellow','pink','blue']
		self.tFS = 25
		self.log = True 
		self.color_key = {'PANTHER': 'darkblue','GORILLA': 'green','ELEPHANT': 'brown','HUMAN':'red','POLAR': 'purple','DOLPHIN': 'gray','TIGER':'orange','ORCA': 'k'}
		self.feats = [] 
		self.vals = [] 
		self.log_vals = [] 
		for line in open(f_data):
			line = line.split() 
			if line[0] == '---': self.samples = line[1::] 
			else:
				self.feats.append(line[0]) 
				self.vals.append([float(x) for x in line[1::]])
				self.log_vals.append([math.log(v+1.0) for v in self.vals[-1]])
 
		self.key = {} 

	def run_norms(self,rand=False):

		sample_vals = [[self.vals[i][j] for i in range(len(self.vals))] for j in range(len(self.vals[0]))]
		rpg,lpg = [], [] 
		for V in sample_vals:
			obs_vals = [v for v in V if v>0] 
			if len(obs_vals) > 0:
				rpg.append(sum(obs_vals)/len(obs_vals)) 
				lpg.append(sum([log(x) for x in obs_vals])/len(obs_vals))
		rpg = sorted(rpg)[0:int(len(rpg)*0.95)] 
		lpg = sorted(lpg)[0:int(len(lpg)*0.95)] 

		RPG = sum(rpg)/len(rpg) 
		LPG = sum(lpg)/len(lpg) 
		MFR,MFL = [],[]
		MOD_SAMPLES = []  
		for n,V in enumerate(sample_vals):
			obs_vals = [v for v in V if v>0] 
			log_vals = [log(v+1.0) for v in V]
			if len(obs_vals) > 0:
				MOD_SAMPLES.append(self.samples[n]) 
				rpg = sum(obs_vals)/len(obs_vals) 
				lpg  = sum([log(x+1) for x in obs_vals])/len(obs_vals)
				mR,mL = RPG/rpg, LPG/(lpg+0.001) 
				
				MFR.append([v *mR for v in V])
				MFL.append([v *mL for v in log_vals])

		R_vals = [[MFR[i][j] for i in range(len(MFR))] for j in range(len(MFR[0]))]
		L_vals = [[MFL[i][j] for i in range(len(MFL))] for j in range(len(MFL[0]))]
		
		print '---',' '.join(MOD_SAMPLES)
		#	
		for i in range(len(R_vals)):
			print self.feats[i],' '.join([str(s) for s in R_vals[i]])	




if __name__ == '__main__':

	import sys,os
	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)
	parser.add_option('--title',  dest= "title", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--pca',  dest= "pca", type = 'str' , default = None, help = "horizontal data")
	(options, args) = parser.parse_args()

	#def sample_creature(w,wv,b,bv,s,sv,c,cv):

	cut1, cut2 = 1.0, 1.0
	if len(args) == 1: 
		for line in open(args[0]):
			line = line.split() 
			if line[0] == '---': continue 	
			line_data = line[0:13]
			s1,s2 = int(line[5]),int(line[9])
			kd = dd(lambda: [0,0]) 

			se1, se2 = float(s1)/(s1+s2), float(s2)/(s1+s2) 

			for i in range(13,len(line),5):
				kx,c1,c2  = line[i],float(line[i+2]),float(line[i+3])
				#ct, se1, se2 = c1+c2, float(s1)/(s1+s2), float(s2)/(s1+s2) 
				
				kd[kx][0] +=c1 
				kd[kx][1] +=c2 

				if kx == 'HP' or kx == 'TP': 
					kd['AD'][0] += c1
					kd['AD'][1] += c2 
				elif kx[-1] == '+':
					kd['BIG'][0] += c1
					kd['BIG'][1] += c2 
				elif kx in ['SVZ','IZ','SP']: 
					kd['REG'][0] += c1
					kd['REG'][1] += c2 
										
			for k,(x1,x2) in kd.items():
				xt     = x1+x2 
				e1,e2  = se1 * xt , se2 * xt 

			#	if x1+x2 < 10: continue 

				if x1 > e1 and x2 < e2: 
					scr1,scr2 = (x1/e1) , (e2)/(x2+1.0)
					if scr1 > cut1 and scr2> cut2:	
						print " ".join(line_data),k,x1,x2,'|',round(scr1,3),round(scr2,3)
				elif x1 < e1 and x2 > e2:
					scr1,scr2 = x2/e2, (e1/(x1+1.0))
					if scr1 > cut1 and scr2> cut2:	
						print " ".join(line_data),k,x1,x2,'|',round(scr1,3),round(scr2,3)


	
	








