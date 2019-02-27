#!/usr/bin/env python
# formula: response ~ predictor + predictor
from collections import Counter as cc
from collections import defaultdict as dd
from math import exp
from math import fabs
from math import factorial 
from math import log
from operator import mul
from random import random
from random import shuffle
from scipy.signal import savgol_filter as svgf 
from scipy.stats import chisquare
from scipy.stats import pearsonr as pearsonr
from scipy.stats import poisson as PSN 
#from scipy.stats import spearmanr as spearmanr
from scipy.stats import spearmanr as spearmanr
from scipy.stats import ttest_ind 
from scipy.stats import variation as coVar 
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import MinMaxScaler
from statsmodels.stats import power as smp 
from statsmodels.stats.multitest import fdrcorrection as fdr
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif 
# import formula api as alias smf
import itertools
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pylab 
import random
import scipy.spatial.distance as sds 
import scipy.stats as stats
import seaborn
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sys

#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 
import warnings
warnings.filterwarnings("ignore")


def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]

def reg_error(msg):
        sys.stderr.write('RageRegTestError: '+msg+'\n')
        sys.exit()




def vif_test(X):

	vd,vd_out = dd(list),{}

	for i,n in enumerate(X.names):
		try: 
			vd_out[n] = round(vif(X.array,i),4) 
		except:
			vd_out[n] = 0.0
	return vd_out



class ListStats:
	def __init__(self,L,goal='min'):
	
		p1 = max(1,int(len(L)*0.01))
	
		self.mean = np.mean(L) 
		self.std  = np.std(L) 
		self.p1 = np.mean(L[0:p1])
		self.p10 = np.mean(L[0:int(p1*10)])
		




class RegModel:
        def __init__(self,options,X,requirements='dex'):

		self.options = options 
		self.X = X 

		if options.model == 'OLS':   self.regress = RegOLS(self.X,requirements) 

		


		self.stats = dd(lambda: {}) 	
                self.rs_key = [0.01,  0.03,   0.05,    0.10,     0.25]
                self.pv_key = [0.01, 0.001, 0.0001, 0.00001,0.0000001]
		self.resid = [] 


	def run_full(self,Y): 

		

		self.reg = [] 
		for y in Y:
			r = self.regress(y,self.X)

			











	def run2(self, Y, X,req=[]):

		self.vif, self.X, self.req = vif_test(X), X, req
		if 'pwr' in req: 	self.reg = [self.regress(y,X).test_pwr() for y in Y]
		else:		        self.reg = [self.regress(y,X) for y in Y]		

		return self


	def score(self, Y, X,req=['resids']):
		self.vif, self.X,rs,self.pv_mins,self.resid = vif_test(X),X, [], [] , [] 
		for y in Y:
			r = self.regress(y,X)
			rs.append(r.model.rsquared)
			try: self.pv_mins.append(min([r.model.pvalues[idx] for idx in self.X.predictor_idx]))
			except ValueError: self.pv_mins.append(min([r.model.pvalues[idx] for idx in self.X.covariate_idx]))
			if 'resids' in req:	self.resid.append(r.model.resid) 

                self.pv_cnt = [len([p for p in self.pv_mins if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
                self.rs_cnt = [len([p for p in rs if p > self.rs_key[j]]) for j in range(len(self.rs_key))]

                for k,K in zip(['rs','pv'],[rs,sorted(self.pv_mins)]): self.stats[k] = ListStats(K) 
		return self

	

	def summarize(self):
		bic =sorted([r.model.bic for r in self.reg],reverse=True)
		rs, ars = sorted([r.model.rsquared for r in self.reg],reverse=True), sorted([r.model.rsquared_adj for r in self.reg],reverse=True)
		pwr_05,pwr_001 = sorted([r.pwr[0.05] for r in self.reg],reverse=True), sorted([r.pwr[0.001] for r in self.reg],reverse=True)

		try:		   self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.predictor_idx]) for r in self.reg]
		except ValueError: 
			try: self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.covariate_idx]) for r in self.reg]
			except: self.pv_mins = [1 for r in self.reg]

                self.pv_cnt = [len([p for p in self.pv_mins if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
                self.rs_cnt = [len([p for p in rs if p > self.rs_key[j]]) for j in range(len(self.rs_key))]
                for k,K in zip(['bic','rs','ars','pv','pwr_05','pwr_001'],[bic,rs,ars,sorted(self.pv_mins),pwr_05,pwr_001]): self.stats[k] = ListStats(K) 
		for i,n in enumerate(self.X.names): 	self.stats[n] =  ListStats(sorted([r.model.pvalues[i] for r in self.reg]))
		return self


	def resids(self):
		if len(self.resid) > 0: return self.resid
		return [r.model.resid for r in self.reg]



	def run_model(self,Y,X,requirements=['dex']):
 
		self.vif, self.X, self.req = vif_test(X), X, req

		for y in Y: 
			reg = self.regress(y,X)




		if 'pwr' in req: 	self.reg = [self.regress(y,X).test_pwr() for y in Y]

		else:		        self.reg = [self.regress(y,X) for y in Y]		

		return self




        def __str__(self):
                return str('Dex Test')

'''
	def out(self):
		rs, ars = sorted([r.model.rsquared for r in self.reg],reverse=True), sorted([r.model.rsquared_adj for r in self.reg],reverse=True)
		try:		   self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.predictor_idx]) for r in self.reg]
		except ValueError: self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.covariate_idx]) for r in self.reg]
                self.pv_cnt = [len([p for p in self.pv_mins if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
                self.rs_cnt = [len([p for p in rs if p > self.rs_key[j]]) for j in range(len(self.rs_key))]
                for k,K in zip(['rs','pv'],[rs,sorted(self.pv_mins)]): self.stats[k] = ListStats(K) 
		return self
'''








	### EACH TIME IT JUST GETS A TWO DIMENSIONAL X ###


	### BUT WAIT HOW ABOUT THINGS LIKE BALANCE ???  ### 

	### BALANCE AND VIF ARE DONE X vs X ### 

	### CHI SQUARE FOR SEG IS NOT ### --- THAT IS DONE ON THE GENE LEVEL ### 


class RegOLS:
        def __init__(self,require = 'dex'):


		self.pwr = {0.05: 0.50, 0.001: 0.1}

		if    require == 'dex': 	self.run = self.run_dex 
		elif  require == 'full':        self.run = self.run_full
		else: require == 'test': 	self.run = self.run_test
		

		self.y,self.X = y, X 

		self.model = sm.OLS(y,X.array).fit() 
		
		
		

		
		




	def test_pwr(self,alphas=[0.05,0.001]):
		if self.model.rsquared > 0: 
			df_de  =  len(self.X.names) -1 
			df_num =  len(self.y) - len(self.X.names)  
			f_2 = np.sqrt(self.model.rsquared / (1-self.model.rsquared) )
	
			for a in alphas:	self.pwr[a] = smp.FTestPower().solve_power(effect_size=f_2,df_num=df_num,df_denom=df_de,alpha=a)

		return self	





			



