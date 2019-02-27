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





class dist:
        def __init__(self,bw=0.1):
		self.kde = KernelDensity(kernel='gaussian', bandwidth=bw)
		self.x_size = 500


	def run(self,vals):

#		self.x_size = int(len(vals)/2.0)

		print len(vals) 
		m_vals = np.array(vals).reshape(-1,1) 
		
		print m_vals 



		minV,maxV = min(vals),max(vals)

		print minV,maxV

		X_plot = np.linspace(-0.1,maxV, len(m_vals))[:, np.newaxis]
		#X_plot = np.linspace(int(minV)-1,int(maxV)+1, self.x_size)[:, np.newaxis]
		y = np.exp(self.kde.fit(m_vals).score_samples(X_plot))
		x = X_plot[:,0] 


		return x,y 	





class samples:
        #def __init__(self,samples,colors,ids,bw=0.2):
        def __init__(self,bw=0.2):
		self.kde = KernelDensity(kernel='gaussian', bandwidth=bw)
		self.x_size = 500

	def run(self,vals):

#		self.x_size = int(len(vals)/2.0)
		m_vals = np.array(vals).reshape(-1,1) 
		minV,maxV = min(vals),max(vals)
		X_plot = np.linspace(int(minV)-1,int(maxV)+1, self.x_size)[:, np.newaxis]
		y = np.exp(self.kde.fit(m_vals).score_samples(X_plot))
		x = X_plot[:,0] 
		return x,y 	

	

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

	#scaler
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
