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
import statsmodels.miscmodels.count as msc
import statsmodels.discrete.count_model as scm 
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

import rage_regtests  
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





class RegModel:
	def __init__(self,X,dist,options,progress=False,FULL = False):

		self.X, self.dist, self.options, self.progress, self.vif  = X, dist, options, progress, vif_test(X) 
		self.alp1, self.alp2 = 0.05, 0.001	


		



		self.regress = rage_regtests.RegTest(X,self.dist,alphas=[self.alp1,self.alp2]) 
		
		#if options.model == 'OLS':   self.regress = RegOLS(X,alphas= [self.alp1, self.alp2]) 
	

		self.out = dd(list) 	
		self.rs_key = [0.01,  0.03,   0.05,    0.10,     0.25]
		self.pv_key = [0.05, 0.01, 0.001, 0.0001, 0.00001,0.0000001]
		# self.pv_key = [0.01, 0.001, 0.0001, 0.00001,0.0000001]
		self.resid = [] 

	def run(self, Y, Ynames=None):
		self.Y = Y 
		t_var,d_var = 0, 0
		if self.progress: self.progress.start_minor('Running '+self.dist+' Regression',len(self.Y))
		for yi,y in enumerate(Y):
			if self.progress: self.progress.mark() 
			t = self.regress.test(y) 
			y_var = np.var(y) 
			t_var += y_var
			d_var += y_var * t.rsq
			
			print 'ok' 
			self.out['params'].append(t.output) 
			self.out['zp'].append(t.zero_prob) 
			if len(self.X.predictor_idx)>0:	 self.out['pvs'].append([p[0] for p in sorted(t.output) if p[-1]][0]) 
#			for r,k in zip([t.rsq, t.rsa, t.bic, t.resids,t.n_resids, t.c_resids,t.pwr[self.alp1],t.pwr[self.alp2]],['rsq','rsa','bic','resids','null_resids','covariate_resids','pwr1','pwr2']): 
			for r,k in zip([t.rsq, t.rsa,t.bic,t.pwr[self.alp1],t.pwr[self.alp2]],['rsq','rsa','bic','pwr1','pwr2']): self.out[k].append(r) 		
		
	
		self.total_var = d_var / t_var 

		return self	


	def aggregate(self,PARAMS=False): 

		if len(self.X.predictor_idx) > 0:	pv_list = [sorted([p[0] for p in P if p[-1]])[0] for P in self.out['params']]
		else:					pv_list = [] 

		self.pv_cnt = [len([p for p in pv_list if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
		self.rs_cnt = [len([p for p in self.out['rsq'] if p > self.rs_key[j]]) for j in range(len(self.rs_key))]
		
		if PARAMS: 	
			self.pv_dict = {self.X.names[i] : [p[i][0] for p in self.out['params']] for i in range(len(self.X.names))} 

		return self 		

	def get_resids(self,COVAR=True):
		y_res = [] 

		self.resids, self.covar_resids = [], [] 
		if self.options.dist == 'OLS': 
			for i,y in enumerate(self.Y): 
				self.covar_resids.append([y[s] - sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names)) if not self.out['params'][i][j][-1]]) for s in range(len(self.X.array))])
				self.resids.append([y[s] - sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names))]) for s in range(len(self.X.array))])


		elif self.options.dist[0].upper() != 'Z':
	
			for i,y in enumerate(self.Y):
				self.resids.append([y[s] - np.exp(sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names))])) for s in range(len(self.X.array))])
				self.covar_resids.append([y[s] - np.exp(sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names)) if not self.out['params'][i][j][-1]])) for s in range(len(self.X.array))])

		else:
			for i,y in enumerate(self.Y): 
				params = self.out['params'][i]
				
				m_res,c_res = [], [] 
				for s in range(len(self.X.array)):

					model_pred = (1-self.out['zp'][i]*self.X.zp[s]) * np.exp(sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names))]))
					covar_pred =  (1-self.out['zp'][i]*self.X.zp[s]) * np.exp(sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names)) if not self.out['params'][i][j][-1]]))

			#		model_pred = (1-self.out['zp'][i]) * np.exp(sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names))]))
				#	covar_pred =  (1-self.out['zp'][i]) * np.exp(sum([self.out['params'][i][j][1]*self.X.array[s][j] for j in range(len(self.X.names)) if not self.out['params'][i][j][-1]]))
					m_res.append(y[s] - model_pred) 
					c_res.append(y[s] - covar_pred)

				
				self.resids.append(m_res)
				self.covar_resids.append(c_res) 
				
		if COVAR: 
			return self.resids, self.covar_resids

	#def run_permutation_tests(self,Y,V,minT=199,midT=1000,maxT=10000): 
	def run_permutation_tests(self,Y,V,minT=10,midT=10,maxT=10): 

		minCut, midCut  = 0.05 * minT, 0.01 * (midT+minT)
		self.permutations = dd(list) 
		for i,y in enumerate(Y):


			p_cands,self.p_key = [[a[2],a[0]] for a in [r for r in self.out['params'][i] if r[-1]]],{} 


			for c,pv in p_cands: 			
				if pv > 0.10: self.permutations[c].append(round(pv,2)) 
				else: 	      self.p_key[c] = [pv,0,0]
		
			
			
			if len(self.p_key) > 0: 

				self.regress.permute(y,[ V.select_variables(V.variables, permute=V.predictors) for j in range(minT)],self.p_key) 

				for c,(pv,L,G) in self.p_key.items(): 
					if L > minCut:    self.permutations[c].append((self.p_key.pop(c)[1])/float(minT))



				if len(self.p_key) > 0: 	
					self.regress.permute(y,[ V.select_variables(V.variables, permute=V.predictors) for j in range(midT)],self.p_key) 
					for c,(pv,L,G) in self.p_key.items(): 
						if L> midCut:    self.permutations[c].append((self.p_key.pop(c)[1])/float(minT+midT))

					if len(self.p_key) > 0: 	
						self.regress.permute(y,[ V.select_variables(V.variables, permute=V.predictors) for j in range(maxT)],self.p_key) 
						for c,(pv,L,G) in self.p_key.items(): 
							self.permutations[c].append((self.p_key.pop(c)[1])/float(minT+midT+maxT))



	### EACH TIME IT JUST GETS A TWO DIMENSIONAL X ###


	### BUT WAIT HOW ABOUT THINGS LIKE BALANCE ???  ### 

	### BALANCE AND VIF ARE DONE X vs X ### 

	### CHI SQUARE FOR SEG IS NOT ### --- THAT IS DONE ON THE GENE LEVEL ### 






class RegTest:
	def __init__(self,X,dist='OLS',alphas=[0.05,0.01],log=True):

		self.X, self.xLen, self.dist, self.permute, self.zero_prob, self.alphas, self.dfd, self.dfn =  X, len(X.names), dist, self.permute_REG, 0.0, alphas, len(X.names) -1,  len(X.array) - len(X.names)  



		#F_KEY = {'TW': sfams.Tweedie(link=slinks.log), 'PO': sfams.Poisson(link=slinks.log), 'NB': sfams.NegativeBinomial(link=slinks.log), 'GA': sfams.Gamma(link=slinks.log), 'NO': sfams.Gaussian(link=slinks.log)}
		F_KEY = {'TW': sfams.Tweedie(), 'PO': sfams.Poisson(), 'NB': sfams.NegativeBinomial(), 'GA': sfams.Gamma(), 'NO': sfams.Gaussian()}

		
		if self.dist.upper()  == 'OLS':  		self.reg, self.execute  = sm.OLS, self.execute_REG
		elif self.dist.upper()[0] == 'G':
			self.reg, self.execute,  self.family =  scm.ZeroInflatedNegativeBinomialP,  self.execute_GIN, F_KEY['NB']

		elif self.dist.upper()[0] != 'Z': 		self.execute, self.permute, self.family = self.execute_GLM , self.permute_GLM, F_KEY[self.dist.upper()[0:2]]
		else:
				 
			if self.dist.upper()[0:3] in ['ZIP', 'ZPO']: 	self.reg, self.execute, self.family  = CUSTOM_ZPO, self.execute_ZIN, F_KEY['PO']
			elif self.dist.upper()[0:3] in ['ZIN', 'ZNB']:    self.reg, self.execute, self.family  = CUSTOM_ZNB, self.execute_ZIN, F_KEY['NB']
			elif self.dist.upper()[0:3] in ['ZGP', 'ZGP']:    self.reg, self.execute, self.family  = CUSTOM_ZGP, self.execute_ZIN, F_KEY['GP']



	def test(self,y): 

		self.y, self.yLen = y, len(y) 
		self.execute() 
		
#		self.resids = RegResiduals(self.X,self.y,self.dist).extract(self.model,self.zero_prob)  
#		self.process() 
		if self.valid: 
			self.output = [(p,b,x,i in self.X.predictor_idx) for i,(p,b,x) in enumerate(zip(self.model.pvalues,self.model.params,self.X.names))]

			try : self.pwr = {a: smp.FTestPower().solve_power(effect_size=np.sqrt(self.rsq/(1-self.rsq)),df_num=self.dfn,df_denom=self.dfd,alpha=a) for a in self.alphas}
			except: self.pwr = {a: 0.5 for a in self.alphas}
		else:
			self.output = [(0.99,0,x,i in self.X.predictor_idx) for i,x in enumerate(self.X.names)]
			self.pwr = {a: 0.5 for a in self.alphas}
			
		

		return self

	def execute_REG(self): 
		self.model, self.null_model = self.reg(self.y,self.X.array).fit(disp=0), self.reg(self.y,self.X.null).fit(disp=0) 

		self.rsq, self.rsa, self.bic, self.valid  = round(self.model.rsquared,5), round(self.model.rsquared_adj,3), round(self.model.bic,3), True



	def execute_ZIN(self): 
		#self.model, self.null_model = self.reg(self.y,self.X.array).fit(disp=0,maxiter=1000),  self.reg(self.y,self.X.null).fit(disp=0) # maxiter=100,method='cg') 
		#self.model, self.null_model = self.reg(self.y,self.X.array).fit(disp=0,maxiter=1000,method='cg'),  self.reg(self.y,self.X.null).fit(disp=0) # maxiter=100,method='cg') 
		#self.model, self.null_model = self.reg(self.y,self.X.array).fit(disp=0,maxiter=2000),  self.reg(self.y,self.X.null).fit(disp=0) # maxiter=100,method='cg') 
		self.model, self.null_model = self.reg(self.y,self.X.array,self.X.zp).fit(disp=0),  self.reg(self.y,self.X.null,self.X.zp).fit(disp=0) # maxiter=100,method='cg') 

		if any([np.isnan(pv) for pv in self.model.pvalues]): self.execute_GLM() 
		else:
			self.zero_prob = 1 / (1 + np.exp(self.model.params[-1]))
			self.bic, self.valid =  self.model.bic, True 
			self.rsq =  1 - (self.model.llf/self.null_model.llf)  
			self.rsa =  1- (((1-self.rsq)*(self.yLen-1)) / (self.yLen-self.X.len-1))

	def execute_GIN(self):

		self.y = np.array([y if i%5 != 0 else y+1 for i,y in enumerate(self.y)])
		self.y = [int(math.log(y+1,2)) for i,y in enumerate(self.y)]

		print self.y[0:20] 


#		self.model = self.reg(np.array(self.y), self.X.array).fit_regularized(disp=0)
#		self.model		 = sm.GLM(np.array(self.y), self.X.array, family= self.family).fit(disp=0)
#		print np.log(np.mean(self.y))
#		print self.model.summary()
#		print self.model.params
#		print self.model.pvalues
		
#		print  1.0 / (1 + np.exp(self.model.params[0]))
#		print len([y for y in self.y if y == 0]), len(self.y) 			
		
		self.model = scm.ZeroInflatedGeneralizedPoisson(np.array(self.y), self.X.array).fit_regularized(disp=0) 
		print self.model.summary()
		print self.model.params
		print self.model.pvalues
#		self.model = CUSTOM_ZNB(self.y, self.X.array).fit(disp=0) 
#		print self.model.summary()
#		print self.model.params
#		print self.model.pvalues

#		print  1.0 / (1 + np.exp(self.model.params[-1]))
		print  1.0 - (1.0 / (1 + np.exp(self.model.params[0])))
		print len([y for y in self.y if y == 0]), len(self.y) 			
		print len([y for y in self.y if y == 0]) / float(len(self.y))		


		sys.exit() 
	 

		
		

	def execute_GLM(self): 
		try: 
			self.model		 = sm.GLM(self.y, self.X.array, family= self.family).fit(disp=0)
			self.null_model 	 = sm.GLM(self.y, self.X.null, family = self.family).fit(disp=0)
			self.bic, self.valid =  self.model.bic, True
			self.rsq =  1 - (self.model.llf/self.null_model.llf)  
			self.rsa =  1- (((1-self.rsq)*(self.yLen-1)) / (self.yLen-self.X.len-1))
		except ValueError: 
			self.bic, self.rsq, self.rsa, self.valid =  100, 0.0, 0.0, False 

		


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







class RegResiduals: 
	def __init__(self,X,y,distribution):

		self.X, self.y, self.yLen, self.distribution = X, y, len(y), distribution  

		if self.distribution.upper()[0:2] == 'OLS': self.extract = self.extract_OLS
		else:					    self.extract = self.extract_GLM

		self.yLog, self.yVar = [np.log(z+1.0) for z in self.y]	, np.var(self.y) 
		self.predicted  =  dd(list) 
		self.subtracted =  dd(list) 
		self.vars       = {} 

	def extract_OLS(self, res, zp = 0):


		self.resids   =   [round(self.y[j]-sum([res.params[k]*self.X.array[j][k] for k in range(len(self.X.array[j]))]),3) for j in range(self.yLen)]
		self.n_resids =   [round(self.y[j]-sum([nes.params[k]*self.X.null[j][k]   for k in range(len(self.X.null[j]))]),3) for j in range(self.yLen)]
		self.c_resids =   [round(self.y[j]-sum([res.params[k]*self.X.array[j][k] for k in self.X.covariate_idx]),3) for j in range(self.yLen)]


	def extract_GLM(self, res, zp = 0):
 
		if   self.distribution.upper()[0:2] == "ZP" and len(res.params) ==  1+len(self.X.names): zp = 1 / (1 + np.exp(res.params[-1]))
		elif self.distribution.upper()[0:2] == "ZN" and len(res.params) ==  2+len(self.X.names): zp = 1 / (1 + np.exp(res.params[-1]))
		else:										 	 zp = 0 
		

		for i in range(self.yLen): 
			model_beta = sum([res.params[k]*self.X.array[i][k] for k in range(len(self.X.names))])
			covar_beta = sum([res.params[k]*self.X.array[i][k] for k in range(len(self.X.names)) if k not in self.X.predictor_idx])
			model_pred = np.exp(model_beta + np.log(1-zp)) 
			covar_pred = np.exp(covar_beta + np.log(1-zp))
			
			model_log_pred = model_beta * (1-zp)  
			covar_log_pred = covar_beta * (1-zp) 

			self.predicted['FULL_MODEL'].append(model_pred) 
			self.predicted['COVAR_MODEL'].append(covar_pred) 

			self.predicted['FULL_LOG_MODEL'].append(model_beta * (1-zp)) 
			self.predicted['COVAR_LOG_MODEL'].append(covar_beta * (1-zp)) 

			self.predicted['LOG_FULL_MODEL'].append(np.log(model_pred+1))
			self.predicted['LOG_COVAR_MODEL'].append(np.log(covar_pred+1))

		MIN_F = min(min(self.predicted['FULL_LOG_MODEL']), min(self.predicted['COVAR_LOG_MODEL']))

		FM_resids = [self.y[j] - self.predicted['FULL_MODEL'][j] for j in range(self.yLen)]
		CO_resids = [self.y[j] - self.predicted['COVAR_MODEL'][j] for j in range(self.yLen)]	
		var_FM, var_CO = np.var(FM_resids), np.var(CO_resids) 

		self.subtracted['FULL_MODEL'], self.subtracted['COVAR_MODEL'] = FM_resids, CO_resids
		self.vars['MODEL'] = [self.yVar, np.var(FM_resids),np.var(CO_resids)]


		FM_LOG_RESIDS = [np.log(self.y[j]+np.exp(MIN_F))-self.predicted['FULL_LOG_MODEL'][j] for j in range(self.yLen)]
		CO_LOG_RESIDS = [np.log(self.y[j]+np.exp(MIN_F))-self.predicted['COVAR_LOG_MODEL'][j] for j in range(self.yLen)]
		var_yLog, var_FMlog, var_COlog = np.var([np.log(self.y[j]+np.exp(MIN_F)) for j in range(self.yLen)]), np.var(FM_LOG_RESIDS), np.var(CO_LOG_RESIDS) 
		
		self.subtracted['FULL_LOG_MODEL'], self.subtracted['COVAR_LOG_MODEL'] = FM_LOG_RESIDS, CO_LOG_RESIDS 
		self.vars['MODEL_LOG'] = [var_yLog, var_FMlog, var_COlog]
		
		
		LOG_FM_RESIDS = [np.log(self.y[j]+1) - self.predicted['LOG_FULL_MODEL'][j] for j in range(self.yLen)]
		LOG_CO_RESIDS = [np.log(self.y[j]+1) - self.predicted['LOG_COVAR_MODEL'][j] for j in range(self.yLen)]
		var_LogY, var_LogFM, var_LogCO =  np.var(self.yLog), np.var(LOG_FM_RESIDS), np.var(LOG_CO_RESIDS) 

		self.subtracted['LOG_FULL_MODEL'], self.subtracted['LOG_COVAR_MODEL'] = LOG_FM_RESIDS, LOG_CO_RESIDS 
		self.vars['LOG_MODEL'] = [var_LogY, var_LogFM, var_LogCO]
		
		return self

























































































class CUSTOM_NB(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(CUSTOM_NB, self).__init__(endog, exog, **kwds)
        
    def nloglikeobs(self, params):
        alph = params[-1]
        beta = params[:-1]
	nY,nX = self.endog, self.exog
    	mu = np.exp(np.dot(nX, beta))
    	size = 1/alph
    	prob = size/(size+mu)
        nloglik        =      -  nbinom.logpmf(nY, size, prob)
	return nloglik
    
    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        # we have one additional parameter and we need to add it for summary
        self.exog_names.append('alpha')
        if start_params == None:
            # Reasonable starting values
            start_params = np.append(np.zeros(self.exog.shape[1]), 1.0)
	    start_params[-1] = 0.5 
 
            # intercept
            start_params[0] = np.log(self.endog.mean())
        return super(CUSTOM_NB, self).fit(start_params=start_params, maxiter=maxiter, maxfun=maxfun, **kwds) 





class CUSTOM_ZPO_back(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, missing='none', **kwds):
        # let them be none in case user wants to use inheritance

        super(CUSTOM_ZPO, self).__init__(endog, exog, missing=missing, **kwds)

        if exog is None:        self.exog = np.ones((self.nobs,1))
        self.nparams = self.exog.shape[1]

        obs_zp = len([e for e in self.endog if e == 0])/float(len(self.endog))
        pred_zp = poisson.pmf(0,self.endog.mean())
        additional_zp = obs_zp - pred_zp 

	if additional_zp > 0.05: 
       		init_z_val = np.log((1.0/additional_zp)-1) 
   		self.start_params = np.append(np.zeros(self.nparams), init_z_val)
                self.exog_names.append('zi')
		self.nloglikeobs = self.nloglikeobs_wzp
	else:
   		self.start_params = np.zeros(self.nparams)
		self.nloglikeobs = self.nloglikeobs_woz 
	
        self.start_params[0] = np.log(self.endog.mean())
        self.cloneattr = ['start_params']

    def nloglikeobs_woz(self, params):
	beta = params[::]
	endog = self.endog 
	XB = np.dot(self.exog, beta)
	nY,nX = self.endog, self.exog
   	nloglik = poisson.logpmf(nY, np.exp(XB))
        return -nloglik

    def nloglikeobs_wzp(self, params):
        beta = params[:-1]
	endog = self.endog 
        gamma = 1 / (1 + np.exp(params[-1]))  #check this
        XB = np.dot(self.exog, beta)
	nY,nX = self.endog, self.exog
	nloglik =            - np.log(1-gamma) - poisson.logpmf(nY,np.exp(XB))
	nloglik[endog==0] =  - np.log(gamma + np.exp(-nloglik[endog==0]))
        return nloglik

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params == None:
            start_params = self.start_params
        return super(CUSTOM_ZPO, self).fit(start_params=start_params, maxiter=maxiter, maxfun=maxfun,**kwds)









class CUSTOM_ZPO(GenericLikelihoodModel):
    def __init__(self, endog, exog, zp = 'none', missing='none', **kwds):
        # let them be none in case user wants to use inheritance

        super(CUSTOM_ZPO, self).__init__(endog, exog, missing=missing, **kwds)

	self.zp = np.array(zp) 
	self.cnter = 0 
        self.nparams = self.exog.shape[1]


        obs_zp = len([e for e in self.endog if e == 0])/float(len(self.endog))
        pred_zp = poisson.pmf(0,self.endog.mean())
        additional_zp = obs_zp - pred_zp 

	if additional_zp > 0.05: 
       		init_z_val = np.log((1.0/additional_zp)-1) 
   		self.start_params = np.append(np.zeros(self.nparams), init_z_val)
                self.exog_names.append('zi')
		self.nloglikeobs = self.nloglikeobs_wzp
	else:
   		self.start_params = np.zeros(self.nparams)
		self.nloglikeobs = self.nloglikeobs_woz 
	
        self.start_params[0] = np.log(self.endog.mean())
        self.cloneattr = ['start_params']

    def nloglikeobs_woz(self, params):
	beta = params[::]
	endog = self.endog 
	XB = np.dot(self.exog, beta)
	nY,nX = self.endog, self.exog
   	nloglik = poisson.logpmf(nY, np.exp(XB))
        return -nloglik

    def nloglikeobs_wzp(self, params):
        beta = params[:-1]
	endog = self.endog 
        gamma = 1 / (1 + np.exp(params[-1]))  #check this
        XB = np.dot(self.exog, beta)
	nY,nX = self.endog, self.exog
#	nloglik =            - np.log(1-gamma) - poisson.logpmf(nY,np.exp(XB))

#	print gamma

	nloglik = -np.log((1-gamma)*self.zp) - poisson.logpmf(nY,np.exp(XB))
	nloglik[endog==0] =  - np.log(gamma*self.zp[endog==0] + np.exp(-nloglik[endog==0]))

#	print nloglik, gamma
#	self.cnter += 1 
#	if self.cnter > 10: sys.exit() 
#	nloglik[endog==0] =  - np.log(gamma + np.exp(-nloglik[endog==0]))
        return nloglik

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params == None:
            start_params = self.start_params
        return super(CUSTOM_ZPO, self).fit(start_params=start_params, maxiter=maxiter, maxfun=maxfun,**kwds)








































class CUSTOM_ZNB(GenericLikelihoodModel):
    def __init__(self, endog, exog, zp='none', missing='none', **kwds):
        super(CUSTOM_ZNB, self).__init__(endog, exog, missing=missing, **kwds)

	self.zp = np.array(zp) 
        self.nparams = self.exog.shape[1]

        obs_zp = len([e for e in self.endog if e == 0])/float(len(self.endog))
	pred_zp = nbinom.pmf(0, 2, 2.0/(2.0+self.endog.mean())) 
        pred_zp = poisson.pmf(0,self.endog.mean())
        additional_zp = obs_zp - pred_zp

 
        self.exog_names.append('alpha')

	if additional_zp > 0.25: 
       		init_z_val = np.log((1.0/additional_zp)-1) 
        	self.start_params = np.hstack((np.zeros(self.nparams), 0.5,init_z_val))
                self.exog_names.append('zi')
		self.nloglikeobs = self.nloglikeobs_wzp
	else:
   		self.start_params = np.zeros(self.nparams)
        	self.start_params = np.hstack((np.zeros(self.nparams), 0.5))
		self.nloglikeobs = self.nloglikeobs_woz 
	
        self.start_params[0] = np.log(self.endog.mean())
        self.cloneattr = ['start_params']

    def nloglikeobs_woz(self, params):

	nY,nX,endog = self.endog, self.exog, self.endog 
        beta, alpha = params[:-1], params[-1]
        mu = np.exp(np.dot(nX, beta))
        size = 1/alpha
        prob = size/(size+mu)
        nloglik        =      -  nbinom.logpmf(nY, size, prob)
	return nloglik 






    def nloglikeobs_wzp(self, params):
	nY,nX,endog = self.endog, self.exog, self.endog 
        beta, alpha = params[:-2], params[-2]
        gamma = 1 / (1 + np.exp(params[-1]))  #check this
        mu = np.exp(np.dot(nX, beta))
        size = 1/alpha
        prob = size/(size+mu)
        nloglik        =      -np.log((1-gamma)*self.zp) -  nbinom.logpmf(nY, size, prob)
        nloglik[nY==0] =      - np.log(gamma*self.zp[nY==0] + np.exp(-nloglik[nY==0]))
        return nloglik

#	nloglik = -np.log((1-gamma)*self.zp) - poisson.logpmf(nY,np.exp(XB))
#	nloglik[endog==0] =  - np.log(gamma*self.zp[endog==0] + np.exp(-nloglik[endog==0]))

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params == None:
            start_params = self.start_params
            start_params[0] = np.log(self.endog.mean())
        return super(CUSTOM_ZNB, self).fit(start_params=start_params, maxiter=maxiter, maxfun=maxfun,**kwds)












class CUSTOM_ZNB_back(GenericLikelihoodModel):
    def __init__(self, endog, exog=None,  missing='none', **kwds):
        super(CUSTOM_ZNB, self).__init__(endog, exog, missing=missing, **kwds)

        if exog is None:        self.exog = np.ones((self.nobs,1))
        self.nparams = self.exog.shape[1]

        obs_zp = len([e for e in self.endog if e == 0])/float(len(self.endog))
	pred_zp = nbinom.pmf(0, 2, 2.0/(2.0+self.endog.mean())) 
        pred_zp = poisson.pmf(0,self.endog.mean())
        additional_zp = obs_zp - pred_zp

 
        self.exog_names.append('alpha')

	if additional_zp > 0.25: 
       		init_z_val = np.log((1.0/additional_zp)-1) 
        	self.start_params = np.hstack((np.zeros(self.nparams), 0.5,init_z_val))
                self.exog_names.append('zi')
		self.nloglikeobs = self.nloglikeobs_wzp
	else:
   		self.start_params = np.zeros(self.nparams)
        	self.start_params = np.hstack((np.zeros(self.nparams), 0.5))
		self.nloglikeobs = self.nloglikeobs_woz 
	
        self.start_params[0] = np.log(self.endog.mean())
        self.cloneattr = ['start_params']

    def nloglikeobs_woz(self, params):

	nY,nX,endog = self.endog, self.exog, self.endog 
        beta, alpha = params[:-1], params[-1]
        mu = np.exp(np.dot(nX, beta))
        size = 1/alpha
        prob = size/(size+mu)
        nloglik        =      -  nbinom.logpmf(nY, size, prob)
	return nloglik 

    def nloglikeobs_wzp(self, params):
	nY,nX,endog = self.endog, self.exog, self.endog 
        beta, alpha = params[:-2], params[-2]
        gamma = 1 / (1 + np.exp(params[-1]))  #check this
        mu = np.exp(np.dot(nX, beta))
        size = 1/alpha
        prob = size/(size+mu)
        nloglik        =      -np.log(1-gamma) -  nbinom.logpmf(nY, size, prob)
        nloglik[nY==0] =      - np.log(gamma + np.exp(-nloglik[nY==0]))
        return nloglik


    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params == None:
            start_params = self.start_params
            start_params[0] = np.log(self.endog.mean())
        return super(CUSTOM_ZNB, self).fit(start_params=start_params, maxiter=maxiter, maxfun=maxfun,**kwds)


















































    
if __name__=="__main__":
    
	#    M = pm.MCMC(locals())
	#    M.sample(10000, 5000, verbose=2)

	

	y1 = [1,2,3,5,0,7,8,12,11,0]
	y2 = [13,21,0,19,33,2,0,0,1,0]
	y3 = [36,0,0,29,33,2,0,0,1,0]
	y4 = [1,2,3,3,3,1,9,3,11,1]

	y5 = [19,0,0,33,0,20,15,0,13,0]
	y5 = [7,13,0,33,0,20,0,11,17,0]
	y6 = [1,2,3,0,3,1,1,0,0,1]
	


	y = y5+y6 
	X = X_data(y)


	ols = RegTest(X).test(y) 
	print "OLS",ols.bic
	for a,b,c,d in ols.output: print c,round(b,3),a
	print "" 




	ps2 = RegTest(X,dist='POISSON').test(y) 
	print "POISSON",ps2.bic 
	for a,b,c,d in ps2.output: print c,round(b,3),a
	print "" 



	nb2 = RegTest(X,dist='NB').test(y) 
	print "NB",nb2.bic
	for a,b,c,d in nb2.output: print c,round(b,3),a
	print "" 

	zp = RegTest(X,dist='ZIP').test(y) 
	print "ZIP",zp.bic
	for a,b,c,d in zp.output: print c,round(b,3),a
	print "" 


	gp = RegTest(X,dist='NOR').test(y) 
	print "NORM",gp.bic
	for a,b,c,d in gp.output: print c,round(b,3),a
	
	print gp.model.summary() 
	print gp.model.scale, 'uh'
	print "" 


	sys.exit() 
#	sys.exit()  



	












