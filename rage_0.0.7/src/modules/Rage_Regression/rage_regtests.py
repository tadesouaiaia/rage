#!/usr/bin/env python
# formula: response ~ predictor + predictor
from collections import Counter as cc
from collections import defaultdict as dd
#from math import exp
#from math import fabs
from scipy.misc import factorial 
#from math import log
#from operator import mul
#from random import random
#from random import shuffle
#from scipy.signal import savgol_filter as svgf 
#from scipy.stats import chisquare
#from scipy.stats import pearsonr as pearsonr
#from scipy.stats import poisson as PSN 
#from scipy.stats import spearmanr as spearmanr
import random
#from scipy.stats import spearmanr as spearmanr
#from scipy.stats import ttest_ind 
#from scipy.stats import variation as coVar 
#from sklearn.linear_model import LogisticRegression
#from sklearn.neighbors import KernelDensity
#from sklearn.preprocessing import MinMaxScaler
from statsmodels.stats import power as smp 
from statsmodels.stats.multitest import fdrcorrection as fdr
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif 
import statsmodels.discrete.count_model as scm 
import statsmodels.discrete.discrete_model as sdm

# import formula api as alias smf
#import itertools
import math
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
import os
#import pandas as pd
#import pylab 
#import random
#import scipy.spatial.distance as sds 
import scipy.stats as stats
#import seaborn
import statsmodels.api as sm
import statsmodels.genmod.families as sfams
import statsmodels.genmod.families.links as slinks
#m.GLM(y,X,family=sm.families.NegativeBinomial(link=statsmodels.genmod.families.links.log)).fit(disp=0)

#import statsmodels.formula.api as smf
import sys
from statsmodels.base.model import GenericLikelihoodModel
from scipy.stats import nbinom
from scipy.stats import poisson

#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 
import warnings
warnings.filterwarnings("ignore")


def scale_vals(vals,f1=0,f2=1):
		scaler = MinMaxScaler(feature_range=(f1,f2))
		return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]





def vif_test(X):

	vd,vd_out = dd(list),{}

	for i,n in enumerate(X.names):
		try: 
			vd_out[n] = round(vif(X.array,i),4) 
		except:
			vd_out[n] = 0.0
	return vd_out









class RegTest:
	def __init__(self,X,dist='OLS',alphas=[0.05,0.01],log=True):

		self.X,self.xLen,self.dist,self.permute,self.zero_infl,self.alphas,self.dfd,self.dfn, self.LOG =  X, len(X.names), dist.upper(), self.permute_REG, 0.0, alphas, len(X.names) -1,  len(X.array) - len(X.names),False
		self.test_history = []  



		if self.dist[0] == 'Z': 
			Z_KEY = {'ZGP':(scm.ZeroInflatedGeneralizedPoisson, scm.GeneralizedPoisson),'ZNB':(scm.ZeroInflatedNegativeBinomialP,scm.NegativeBinomialP),'ZPO':(scm.ZeroInflatedPoisson,scm.Poisson)}
			self.execute, self.permute = self.execute_zif, self.permute_ZIF
			self.zi_model, self.nz_model = Z_KEY[self.dist]	
		elif self.dist in ['NB','GP','PO']:
			F_KEY = {'GP': scm.GeneralizedPoisson,'NB':scm.NegativeBinomialP,'PO':scm.Poisson}
			self.execute, self.nz_model, self.permute = self.execute_cm , F_KEY[self.dist], self.permute_CM

		elif self.dist[0:3] == 'OLS': 
			self.execute, self.LOG, self.permute = self.execute_ols, (dist[-3::].upper() == 'LOG'), self.permute_OLS 
			
		else:
			print "UNSUPPORTED"
			sys.exit() 


	def test(self,y): 

		#print self.LOG,'huh' 

		#if self.LOG: y = [math.log(yi+1.0,2) for yi in y]



		self.output,self.zero_infl, self.rsq, self.rsa, self.bic, self.aic  = [],0.0,'NA','NA','NA','NA'
		self.valid, self.y, self.yA, self.yLen, self.history  = True, y, np.array(y), len(y), '' 
		self.execute() 

		if self.valid: 
			self.v_explained = 1-(np.var(self.res.resid) / np.var(self.yA)) 
			try : self.pwr = {a: smp.FTestPower().solve_power(effect_size=np.sqrt(self.v_explained/(1-self.v_explained)),df_num=self.dfn,df_denom=self.dfd,alpha=a) for a in self.alphas}
			except: self.pwr = {a: 0.5 for a in self.alphas}
			if any([np.isnan(pw) for pw in self.pwr.values()]): self.pwr = {a: 0.5 for a in self.alphas}
		else:
			self.v_explained = 0 
			self.pwr = {a: 0.0 for a in self.alphas} 
			self.output = [(0.5,b,t,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues,self.res.tvalues, self.res.params, self.X.names))]
			self.bic, self.aic, self.rsq, self.rsa = 0, 0, 0, 0 
			 
			self.tvalues = self.res.tvalues

			#self.bic, self.aic, self.rsq, self.rsa = self.res.bic, self.res.aic, self.res.prsquared, 1- (((1-self.res.prsquared)*(self.yLen-1)) / self.dfn) #(self.yLen-self.X.len-1))
			#self.output = [(p,b,x,i in self.X.predictor_idx) for i,(p,b,x) in enumerate(zip(self.res.pvalues, self.res.params, self.X.names))]

		return self


	def execute_ols(self): 

		self.history += 'ols-fit'
		self.res = sm.OLS(self.yA, self.X.array).fit(disp=0) 

		self.fval,self.fpv, self.tvalues = self.res.fvalue, self.res.f_pvalue, self.res.tvalues
		self.rsq, self.rsa, self.bic, self.aic  = round(self.res.rsquared,5), round(self.res.rsquared_adj,3), round(self.res.bic,3), round(self.res.aic)
		 
		self.output = [(p,t,b,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues, self.res.tvalues, self.res.params, self.X.names))]
		self.history += '-sucess'

	def execute_zif(self): 
		self.reg = self.zi_model(self.yA, self.X.array,exog_infl=self.X.zp) 
		self.history += 'zif-fit'

		try: 
			self.res = self.reg.fit(disp=0) 	
			if any([np.isnan(pv) for pv in self.res.pvalues]): 
				self.history += '-fail,zif-regularized'
				self.res = self.reg.fit_regularized(disp=0) 
				if any([np.isnan(pv) for pv in self.res.pvalues]): raise np.linalg.linalg.LinAlgError('NAN Vals') 
		except np.linalg.linalg.LinAlgError:	self.history += '-fail'
		except AssertionError:			self.history += '-fail'
		if self.history[-4::] == 'fail': 
			self.execute_cm() 
			return 
	
		
		self.history += '-sucess'
		self.v_explained = 1-(np.var(self.res.resid) / np.var(self.yA)) 
		self.zero_infl = 1.0 - (1.0 / (1 + np.exp(self.res.params[0])))	
		self.tvalues = self.res.tvalues
		self.bic, self.aic, self.rsq, self.rsa = self.res.bic, self.res.aic, self.res.prsquared, 1- (((1-self.res.prsquared)*(self.yLen-1)) / self.dfn) #(self.yLen-self.X.len-1))
		self.output = [(p,t,b,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues[1::],self.res.tvalues,self.res.params[1::], self.X.names))]

			



	def execute_cm(self):

		self.history += ',cm-fit'
		FAIL = False
		try: 
			self.res = self.nz_model(self.yA, self.X.array).fit(disp=0)
			if any([np.isnan(pv) for pv in self.res.pvalues]): raise np.linalg.linalg.LinAlgError('NAN Vals') 
		except np.linalg.linalg.LinAlgError:	self.history += '-fail'
		except AssertionError:			self.history += '-fail'
			

		if self.history[-4::] == 'fail': 
			self.history += ',cm-regfit'
			try: 
				self.res = self.nz_model(self.yA, self.X.array).fit_regularized(disp=0)
				if any([np.isnan(pv) for pv in self.res.pvalues]): raise np.linalg.linalg.LinAlgError('NAN Vals') 
			except np.linalg.linalg.LinAlgError:	self.history += '-fail'
			except AssertionError:			self.history += '-fail'

		
		if self.history[-4::] == 'fail':
			self.valid = False 
			return 
		
		if self.history[-4::] != 'fail': 
 
			self.history += '-sucess'
			self.tvalues = self.res.tvalues

			self.bic, self.aic, self.rsq, self.rsa = self.res.bic, self.res.aic, self.res.prsquared, 1- (((1-self.res.prsquared)*(self.yLen-1)) / self.dfn) #(self.yLen-self.X.len-1))
			self.output = [(p,t,b,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues, self.res.tvalues,self.res.params, self.X.names))]
 



	def permute_OLS(self,y,Xp,key):



 
		for P in Xp:
			
			#print P.array
 
			#res = sm.OLS(np.array(y),Xp[0].array).fit(disp=0) 
			res = sm.OLS(np.array(y),P.array).fit(disp=0) 
			for i,(p,n) in enumerate(zip(res.pvalues,  self.X.names)):
				if n in key:

 
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1

				
	def permute_CM(self,y,Xp,key): 
		for P in Xp:
			try:  
				res = self.nz_model(np.array(y), P.array).fit(disp=0)
			except:	
				try: res = self.nz_model(np.array(y), P.array).fit_regularized(disp=0)

				except: return 
					

			
			for i,(p,n) in enumerate(zip(res.pvalues, self.X.names)):
				if n in key: 
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1


	def permute_ZIF(self,y,Xp,key): 
		for P in Xp: 
			try: res = self.zi_model(self.yA, P.array,exog_infl=self.X.zp).fit(disp=0) 
			except: 
				try: res = self.nz_model(np.array(y), P.array).fit(disp=0)
				except: 
					try:    res = self.nz_model(np.array(y), P.array).fit_regularized(disp=0)
					except: res.pvalues = [0.5 for x in self.X.names]



			for i,(p,n) in enumerate(zip(res.pvalues, self.X.names)):
				if n in key: 
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1





	def permute_REG(self,y,Xp,key): 		
		for P in Xp: 

			


			model = self.reg(y, Xp[0].array).fit(disp=0)
			for i,(p,b,n) in enumerate(zip(model.pvalues,model.params,self.X.names)):
				if n in key:
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1

	def permute_GLM(self,y,Xp,key): 
		for P in Xp: 
			model = sm.GLM(y, P.array, family= self.family).fit(disp=0)
			for i,(p,b,n) in enumerate(zip(model.pvalues,model.params,self.X.names)):
				if n in key:
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1





