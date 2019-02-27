#!/usr/bin/env python


import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
#import scipy.stats as stats
#from scipy.stats import variation as coVar 
#import scipy.spatial.distance as sds 
#from random import random
#import numpy as np
#import itertools
import random
#from math import fabs
#from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
#from scipy.stats import poisson as PSN 
#from scipy.stats import spearmanr as spearmanr
import seaborn
from math import log
import math
from random import shuffle
#from sklearn.cluster import KMeans	
#from sklearn.cluster import KMeans
#from sklearn.neighbors import KernelDensity
#from sklearn.preprocessing import MinMaxScaler
import matplotlib
import matplotlib.pyplot as plt

#import statsmodels.api as sm
#from statsmodels.stats import power as smp 
#from statsmodels.stats.outliers_influence import variance_inflation_factor as vif 
#from operator import mul

#from scipy.signal import savgol_filter as svgf 
#from math import exp
#from math import factorial 

#from scipy.stats import chisquare
#from scipy.stats import ttest_ind 
from Rage_IO import rage_outputs
from Rage_Plots import rage_regression_plots
from Rage_Transforms import rage_KDE
from Rage_Transforms import rage_DR
from Rage_Regression import rage_dextests as rt  
from Rage_Filters    import rage_filters 



#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 
import warnings
warnings.filterwarnings("ignore")


def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))




def regression_error(msg):
        sys.stderr.write('RageRegressionError: '+msg+'\n')
        sys.exit()



def compare_principal_components(pca1,pca2):

	print 'yo' 




def iterate_t_tests(c_groups):

	t_out = {}
	for i,c1 in enumerate(c_groups): 
		cdiff = [a for b in [c_groups[k] for k in [k for k in c_groups.keys() if k != c1]] for a in b]
		tv = ttest_ind(c_groups[c1],cdiff)[1] 
		m1,m2 = np.mean(c_groups[c1]),np.mean(cdiff)
		
		if m1 > m2: fc = m1 / (m2+0.001) 
		else:       fc = (-1*m2) / (m1+0.001)
		t_out[c1] = [tv,round(m1,3),round(fc,3)]
	return t_out 



def set_result_counter(model_result,model_type = 'full'):



	if model_type == 'full': 
		pvs, my_rs =  sorted([model_result['params'][i]['predictors'][0][0] for  i in range(len(model_result['params']))]), model_result['rs']
	else:
		pvs,my_rs = model_result['PV'][0], model_result['RS'][0]

        steps,rsk, maxR, pv5 = 5, [0.01,0.02], round(max(my_rs),2), np.percentile(pvs,10)
        if pv5 < 0.001: my_key = [0.01]
        elif pv5<0.01:  my_key = [0.01,0.005]
        else:           my_key  = [0.05,0.01,0.005]
        p,pv = pvs[0],0.001
        while True:
        	if p < pv: my_key.append(pv)
                else:      break
                pv /= 10.0
                if pv < 0.0000001: break
	my_counts = [0 for m in my_key]
	for i,m in enumerate(my_key):	my_counts[i] = len([p for p in pvs if p < m])

        step = round((maxR-rsk[-1])/steps,2)
        bar_key = sorted(list(set(rsk+[round(rsk[-1]+((1+i)*step),2) for i in range(steps+1)])))
	rs_counts = [0 for m in bar_key]
	for i,m in enumerate(bar_key):	rs_counts[i] = len([p for p in my_rs if p > m])

	return my_counts,rs_counts, my_key, bar_key





































class RegressionResult:
        def __init__(self):


		self.features = [] 
		self.result, self.sims = {} , {} 
		self.sim_vars, self.sim_rs, self.sim_pv = [],[] , [] 
		self.sim_vars_cnt, self.sim_rs_cnt, self.sim_pv_cnt = [],[] , [] 




	def add_full_result(self,f_dict): 
		self.result = f_dict 
		self.p1 = int(0.999+len(self.result.keys()) * 0.01)
		return self


	def summarize(self,X,steps=8):




		M, self.vif  = self.result.values() , rt.vif_test(X) 			
		self.rs, self.ars, self.bic  =   sorted([m['ars'] for m in M],reverse=True),sorted([m['rs'] for m in M],reverse=True), [m['bic'] for m in M]
		self.resid, self.pwr_05, self.pwr_001    =   [m['resids'] for m in M],   sorted([m['pwr-05'] for m in M],reverse=True), sorted([m['pwr-001'] for m in M], reverse=True) 

		if 'covariates' in m:  self.c_resid = [m['covariates']['resids'] for m in M]
		else: 		       self.c_resid = self.resid
		self.min_pvs  = sorted([sorted([m['params'][k][0][0] for k in m['params'].keys() if not X.COVARIATE[k]])[0] for m in M])
	
		if max(self.rs) < 0.1: self.rs_key = [0.01,0.02,0.03,0.04,0.05] 
		else: 		       self.rs_key = [0.01,0.05,0.10,0.25,0.50] 
		if min(self.min_pvs) < 0.0000001: self.pv_key = [0.01, 0.001, 0.0001, 0.00001,0.0000001] 
		else:				  self.pv_key = [0.05, 0.01, 0.001, 0.0001, 0.00001]
		
		self.pv_cnt = [len([p for p in self.min_pvs if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
		self.rs_cnt = [len([p for p in self.rs if p > self.rs_key[j]]) for j in range(len(self.rs_key))]

		self.stats, self.predictor_stats = {} , {} 

		for k,K in zip(['bic','rs','ars','pv','pwr-05','pwr-001'],[self.bic,self.rs,self.ars,self.min_pvs,self.pwr_05,self.pwr_001]): 
			self.stats[k] = {'mean': np.mean(K), 'std': np.std(K), 'p10': np.mean(K[0:10*self.p1]), 'p1': np.mean(K[0:self.p1])}
	
		self.pvs, self.predictor_pvs = dd(list), dd(list)   
		for m in M:
			for k,K in m['params'].items():
				for kPv,kBeta,kName in K: 
					if len(set([k,kName]))>1:	self.pvs[k+'='+kName].append(kPv) 
					else: 	       		   	self.pvs[k].append(kPv) 
			for k,K in m['predictors']['params'].items():
				for kPv,kBeta,kName in K:
					if len(set([k,kName]))>1: 	self.predictor_pvs[k+'='+kName].append(kPv) 
					else: 	       			self.predictor_pvs[k].append(kPv) 
		for k,K in self.pvs.items():		self.stats[k] = {'mean': np.mean(K), 'std': np.std(K), 'p10': np.mean(K[0:10*self.p1]), 'p1': np.mean(K[0:self.p1])}
		for k,K in self.predictor_pvs.items():	self.predictor_stats[k] = {'mean': np.mean(K), 'std': np.std(K), 'p10': np.mean(K[0:10*self.p1]), 'p1': np.mean(K[0:self.p1])}
		

		return self


	def add_simulated_result(self,sims,X):


		self.sims.append(sims) 

		dim_red = rage_DR.DR(None,False,len(sims))	
		self.sim_rs.append([sim['rs'] for sim in sims][0:5]) 		
		self.sim_pv.append(sorted([sorted([sim['params'][k][0][0] for k in sim['params'].keys() if not X.COVARIATE[k]])[0] for sim in sims]))
		self.sim_vars.append(dim_red.run_pca(np.matrix([sim['resids'] for sim in sims]).getT(),req='brief')['var_rates'])

		self.sim_pv_cnt.append([len([p for p in self.sim_pv[-1] if p < self.pv_key[j]]) for j in range(len(self.pv_key))])
		self.sim_rs_cnt.append([len([p for p in self.sim_rs[-1]  if p > self.rs_key[j]]) for j in range(len(self.rs_key))])
		


	def summarize_sims(self):




		self.var = np.mean(self.sim_vars) 
		r_out= [np.mean([rs[j] for rs in self.sim_rs_cnt]) for j in range(len(self.rs_cnt))]
		p_out= [np.mean([ps[j] for ps in self.sim_pv_cnt]) for j in range(len(self.pv_key))]

		return p_out,r_out,self.var




class Regression:
        def __init__(self,rage):

		self.rage = rage 
		self.progress = rage.progress
		self.options = rage.args 

		self.results = [] 

		if self.options.model == 'OLS':				self.regress = rt.regress_ols	
		elif self.options.model.upper() in ['GLM-NB']:		self.regress = self.regress_glmnb
		elif self.options.model.upper() in ['NBINOM','NB']:	self.regress = self.regress_nb 
		elif self.options.model.upper() in ['ZIP']: 		self.regress = self.regress_zip

	def run(self): 

		self.D = self.rage.data.filter_samples_by_attributes(self.options.predictors,self.options.covariates).normalize() 

		self.V = self.D.set_sample_variables(combine=False)

		#set_sample_variables(combine=False)

		if self.options.command == 'eval-model': 
    			self.progress.start_major('Running Regression Evaluation')
			self.evaluate_model() 	

		elif self.options.command == 'dex':
    			self.progress.start_major('Running Regression')
			self.run_model()

		elif self.options.command == 'eval-covariates':
    			self.progress.start_major('Running Regression Covariate Analysis')	
			self.eval_covariates() 
			
		else: print self.options.command		


	def run_model(self):
		#seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'pink'})
                self.progress.start_minor('Running Model Regression',len(self.D.features),True)
		X = self.V.select_variables()  

		model_result = self.regression_result(X,True,mtype = 'full') 

		sys.exit() 
		rage_outputs.regression_result(self.rage.args,model_result['dex']).write(self.D.predictor,self.D.covariates)
		self.progress.end() 
		sys.exit() 




	def run_regression_model(self,X,req = []):



		regResult = RegressionResult().add_full_result({f : self.regress([s.cnts[f.idx] for s in self.D.samples],X,req)}).summarize() 
		
		regResult.summarize(X) 
		return regResult 

		 

		for f in self.D.features:
			self.progress.mark() 
			res[f] = self.regress([s.cnts[f.idx] for s in self.D.samples],X,req)
		
		vif  = rt.vif_test(X) 


		print res.keys() 

		

	




	def evaluate_model(self,simulations=3):
		X = self.V.select_variables()  
		Y = [[s.cnts[f.idx] for s in self.D.samples] for f in self.D.features]


                self.progress.start_minor('Running Model Regressions',len(self.D.features),False)
		self.dim_red = rage_DR.DR(self.options,False,len(self.D.samples)) 
		pca_init = self.dim_red.run_pca(self.D.matrix(),req='brief')


		req = ['full','pwr','resids','predictors-only','covariates-only']

		M = RegressionResult().add_full_result({f : self.regress([s.cnts[f.idx] for s in self.D.samples],X,req) for f in self.D.features}).summarize(X) 
		rage_outputs.regression_result(self.rage.args).write(M,X) 
		pca_c_resid, pca_resid =    self.dim_red.run_pca(np.matrix(M.c_resid).getT(),req='brief'), self.dim_red.run_pca(np.matrix(M.resid).getT(),req='brief')

		for n in range(simulations): 
			Xs = self.V.select_variables(shuffle_items = [self.options.predictors[0]])
			M.add_simulated_result([self.regress(y,Xs,req=['resids']) for y in Y],Xs)
		sim_pvs,sim_rs,sim_vars = M.summarize_sims() 
 

                if self.options.featureKey == None:     mplot.add_model_table(M,X,{'sim_pvs': sim_pvs}).update()
                else:                                   mplot.add_model_table(M,X,{'ss':   self.calculate_ss(self.D.features)}).update()

 
                self.progress.start_minor('Plotting Results  ',100,False)
		skrees = [pca_init['var_rates'],pca_c_resid['var_rates'],pca_resid['var_rates'],sim_vars] 
		mplot.add_predictor_table(M,X,skrees).update()
	

		mplot.add_rs_bars(M.rs_cnt,M.sim_rs_cnt).update({'title': '$'+"\_".join(self.options.predictors)+'$ '+'$\  R^2\ Values$'})
		mplot.add_pv_bars(M.pv_cnt,M.sim_pv_cnt).update({'title': '$'+"\_".join(self.options.predictors)+'$ '+'$\  P\ \ Values$'})
		mplot.add_pca_pts(pca_init['pts'],self.D.samples,self.options.predictors[0],{'colspan':2}).update({'title': 'PCA Initial Values','yadd': 2,'colspan':2})
		mplot.add_pca_pts(pca_c_resid['pts'],self.D.samples,self.options.predictors[0],{'colspan':2}).update({'title': 'PCA Covariate Residuals','yadd': 2,'colspan':2})
		mplot.add_pca_pts(pca_resid['pts'],self.D.samples,self.options.predictors[0],{'colspan':2}).update({'title': 'PCA Model Residuals','yadd': 2,'colspan':2})
		mplot.save_mfig(self.options.model,self.options.predictors,self.options.covariates)

		self.progress.end() 















	def calculate_balance(self,cTypes,predictor_ids,cShort,cSpec):



		if cTypes == ['binary','binary']:	
			cLen,c_obs,chi_pv,chi_over =self.calculate_chi_enrichment([predictor_ids[j] for j in range(len(predictor_ids)) if self.D.samples[j].attributes[cShort] == cSpec])
			return chi_pv,chi_over

		elif cTypes == ['continuous','binary']:	
			try: 
				btop = sorted(iterate_t_tests({k: [self.D.samples[j].attributes[cShort] for j in G] for k,G in self.prc_seg[0].items()}).items(),key=lambda x: x[1][0])[0]
				return btop[1][0],btop[0] 
			except TypeError:
				return 0.01,'NA'


		elif cTypes == ['continuous','continuous']:
			balance  = pearsonr(predictor_ids,[s.attributes[cShort] for s in self.D.samples])
			return balance[1],str(round(balance[0],2))
		else:
			pMatch = [predictor_ids[j] for j in range(len(predictor_ids)) if self.D.samples[j].attributes[cShort] == cSpec]

			pDiff = [predictor_ids[j] for j in range(len(predictor_ids)) if self.D.samples[j].attributes[cShort] != cSpec]
			balance= [ttest_ind(pMatch,pDiff)[1],round(np.mean(pMatch),3),round(np.mean(pDiff),3)] 				
			return balance[0],'Div'


		return balance

































	def eval_covariates(self,number_of_sims=1):

		if len(self.options.covariates) == 0: regression_error('Eval Covariates requires covariates!') 		

		self.D = self.rage.data.filter_samples_by_attributes(self.options.predictors,self.options.covariates).set_sample_variables()
		X_names,  X = self.D.create_regression_array() 
		self.baseType, S, DVC, fLen,c_tests = 'continuoues',self.D.samples, self.D.variable_class, len(self.D.features), {} 


                self.progress.start_minor('Testing Base Predictor Model',len(self.D.features))


		Xi_sims, Xi,Xi_names = [], [[x[i] for i in range(len(x)) if self.D.variable_class[X_names[i]] != 'covariate'] for x in X], [n for n in X_names if self.D.variable_class[n] != 'covariate']
		Xi_result = self.regression_result(Xi_names,Xi,True, mtype = 'brief')
		self.vifs = self.score_vif(X,X_names)	
		predictor_ids = [s.attributes[self.D.predictor] for s in S]
		if S.attribute_class[self.D.predictor] == 'binary':
			self.prc_rates, self.prc_seg, self.baseType = {pid: pc / float(len(predictor_ids)) for pid,pc in cc(predictor_ids).items()}, S.segregate(self.D.predictor), 'binary'
	
		c_fin, cIdxs  = {}, [i for i,n in enumerate(X_names) if self.D.variable_class[n] == 'covariate']
		for i in cIdxs:
			cName,cVif,cType  = X_names[i],  self.vifs[X_names[i]], S.attribute_class[X_names[i].split('=')[0]]
			cShort,cSpec,cTypes = cName.split('=')[0],cName.split('=')[-1], [cType,S.attribute_class[self.D.predictor]]
                	self.progress.start_minor('Testing covariate: '+X_names[i]+' type='+cType+'..',len(self.D.features))


			
			
			if cType == 'binary': cSize = len([s for s in S if s.attributes[cShort] == cSpec])
			else:		      cSize = len([s for s in S if s.attributes[cShort] != 'NA'])



			c_save = {'vif': cVif,'cType': cType, 'cSize': cSize} 
			n_sims, c_sims, c_discs = self.run_shuffle_simulations(number_of_sims,sim_type='covariate-test',covar_name=cName)
			Xi_sims.extend(n_sims) 
			c_save['balance'] = self.calculate_balance(cTypes,predictor_ids,cShort,cSpec)

			Xc,Xc_names =  [[x[i],1.0] for x in X],[cName,'intercept']
			c_result = self.regression_result(Xc_names,Xc,True, mtype = 'brief')

			Xm,Xm_names = [[x[j] for j in range(len(x)) if (j == i) or (DVC[X_names[j]] != 'covariate')] for x in X], [X_names[j] for j in range(len(X_names)) if (j == i) or (DVC[X_names[j]] != 'covariate')]
			cm_result = self.regression_result(Xm_names,Xm,True, mtype = 'brief')

			c_save['stats'] = {'PV': c_result['PV'][1], 'RS': c_result['RS'][1], 'ARS': cm_result['RS'][1], 'APV': cm_result['PV'][1]} 

			c_discovery = {0.0001: [0,0], 0.001: [0,0], 0.05: [0,0]}
			for npv,cpv in zip(Xi_result['PV'][0],cm_result['PV'][0]):
				if npv != cpv:
					if min(cpv,npv) < 0.0001:	c_discovery[0.0001][npv<cpv] +=1
					elif min(cpv,npv) < 0.001:	c_discovery[0.001][npv<cpv] +=1
					elif min(cpv,npv) < 0.05:	c_discovery[0.05][npv<cpv] +=1
			sim_discovery = {k: [np.mean([cd[k][0] for cd in c_discs]),np.mean([cd[k][1] for cd in c_discs])] for k in c_discovery.keys()}
			c_save['discovery'] = (c_discovery,sim_discovery) 
			Xi_sims.extend(n_sims) 
			c_fin[cName] = c_save

		pv_counts, rs_counts, pv_key, rs_key = set_result_counter(Xi_result,model_type='brief') 
		eplot = rage_regression_plots.eplot(self.options,2+len(c_fin.keys()),{'r_key': rs_key, 'p_key': pv_key}) 
		eplot.add_base_model(self.baseType,self.D.predictors,[Xi_result,Xi_names,pv_counts,rs_counts],Xi_sims)
		eplot.add_covariate_data(c_fin)
		eplot.save_efig(self.options.model,[self.D.predictor],self.D.covariates)






































































































































	def run_shuffle_simulations(self,sims,sim_type,covar_name=None,compare_key={}):
		sim_vars,sim_rs, sim_pvs = [], [], [] 
		n_sim, c_sim, c_discs = [],[], []  



		if sim_type == 'full':
			for iX,(S_names,Xs) in enumerate([self.D.create_regression_array(shuffle_variables=True) for jX in range(sims)]):

				sim_result = self.regression_result(S_names,Xs,True)
				sim_vars.append(self.dim_red.run_pca(np.matrix(sim_result['resids']).getT(),req='brief')['var_rates'])
				sim_pvs.append(sorted([sim_result['params'][i]['predictors'][0][0] for  i in range(len(sim_result['params']))]))
				sim_rs.append(sim_result['rs'])

			return sim_vars,sim_rs,sim_pvs

		elif sim_type == 'compare' and len(compare_key) !=0:

			for iX,(S_names,Xs) in enumerate([self.D.create_regression_array(shuffle_variables=True) for jX in range(sims)]):

                		self.progress.start_minor('Starting Simulation '+str(iX+1),5000,False)
				sim_result = self.regression_result(S_names,Xs,True)
				sim_vars.append(self.dim_red.run_pca(np.matrix(sim_result['resids']).getT(),req='brief')['var_rates'])
				sim_pv = ([sim_result['params'][i]['predictors'][0][0] for  i in range(len(sim_result['params']))])
				sim_r = sim_result['rs'] 
				cnt_p,cnt_r = [0 for p in compare_key['pv']], [0 for p in compare_key['rs']] 
				for j,m in enumerate(compare_key['pv']): cnt_p[j] += len([p for p in sim_pv if p < m]) 
				for j,m in enumerate(compare_key['rs']): cnt_r[j] += len([p for p in sim_r if p > m]) 
				sim_rs.append(cnt_r); sim_pvs.append(cnt_p) 	

			return np.mean(sim_vars),sim_rs,sim_pvs
			r_out= [np.mean([rs[j] for rs in sim_rs]) for j in range(len(compare_key['rs']))]
			p_out= [np.mean([ps[j] for ps in sim_pvs]) for j in range(len(compare_key['pv']))]

			return np.mean(sim_vars), r_out, p_out
			sys.exit() 



		elif sim_type == 'covariate-test' and covar_name != None:		

			for iXc,(S_names,Xs) in enumerate([self.D.create_regression_array(shuffle_variables=True) for jX in range(sims)]):
                		self.progress.start_minor('Starting Simulation '+str(iXc+1),5000,False)
				
				Xic = [[x[i] for i in range(len(x)) if (self.D.variable_class[S_names[i]] != 'covariate') or (S_names[i] == covar_name)] for x in Xs]
				Xic_names = [n for n in S_names if (self.D.variable_class[n] != 'covariate') or (n == covar_name)]	
				c_result = self.regression_result(Xic_names,Xic,True,mtype='brief')

				Xis = [[x[i] for i in range(len(x)) if self.D.variable_class[S_names[i]] != 'covariate'] for x in Xs]
				Xis_names = [n for n in S_names if self.D.variable_class[n] != 'covariate']
				p_result = self.regression_result(Xis_names,Xis,True,mtype='brief')
				c_discovery = {0.0001: [0,0], 0.001: [0,0], 0.05: [0,0]}
				for npv,cpv in zip(p_result['PV'][0],c_result['PV'][0]):
					if npv != cpv:
						if min(cpv,npv) < 0.0001:	c_discovery[0.0001][npv<cpv] +=1
						elif min(cpv,npv) < 0.001:	c_discovery[0.001][npv<cpv] +=1
						elif min(cpv,npv) < 0.05:	c_discovery[0.05][npv<cpv] +=1
						else: continue 
				n_sim.append(p_result)
				c_sim.append(c_result) 
				c_discs.append(c_discovery) 

			return n_sim,c_sim,c_discs













































	def calculate_chi_enrichment(self,pc_ids,pc_ids2 = []):
		if len(pc_ids2) == 0:		
			cLen = len(pc_ids) 
			c_cc = cc(pc_ids) 
			c_exp = [cLen*self.prc_rates[k] for k in self.prc_rates]
			c_obs = [c_cc[k] if k in c_cc else 0 for k in self.prc_rates]

			chi_over = sorted([ (co-ce,k) for co,ce,k in zip(c_obs,c_exp,self.prc_rates)])[-1][1]
			

			chi_pv = chisquare(c_obs,f_exp=c_exp)[1]
			return cLen,c_obs,chi_pv,chi_over









	def enrichment_and_fold_change(self,seg_dict,min_valid=15):



		seg_lens = {k: len(V) for k,V in seg_dict.items()}
		sum_lens = float(sum(seg_lens.values()))
		seg_cnts =  sorted([a for b in [[(c,k) for c in seg_dict[k]] for k in seg_dict.keys()] for a in b])
		seg_means = sorted([(k, np.mean(V)) for k,V in seg_dict.items()], key = lambda x: x[1])
		seg_obs   = sorted([(k, len([v for v in V if v>0])/float(len(V))) for k,V in seg_dict.items()])
		seg_min,seg_max = seg_means[0][0], seg_means[-1][0] 

		seg_valid = len([x[1] for x in seg_cnts if x[0] > 0]) 


		if seg_valid < min_valid  or seg_valid < (seg_lens[seg_max] / 5.0): return {a: b for a,b in seg_means},{a: b for a,b in seg_obs},(seg_min,seg_max),(1.0,1.0)

		seg_means = sorted([(k, np.mean(V)) for k,V in seg_dict.items()], key = lambda x: x[1])
		seg_obs   = sorted([(k, len([v for v in V if v>0])/float(len(V))) for k,V in seg_dict.items()])
		seg_min,seg_max = seg_means[0][0], seg_means[-1][0] 
		min_len,max_len = seg_lens[seg_min], seg_lens[seg_max]
		min_seg = seg_cnts[0:min_len]
		i = len(min_seg) 
		while min_seg[-1][0] == seg_cnts[i][0]: 
			min_seg.append(seg_cnts[i])
			i+=1
			if i == len(seg_cnts): return {a: b for a,b in seg_means},{a: b for a,b in seg_obs},(seg_min,seg_max),(1.0,1.0)
		if max_len > len(seg_cnts) - i: 	max_seg = seg_cnts[i::] 
		else:
			seg_rev = seg_cnts[-1::-1]
			max_seg = seg_rev[0:max_len]
			i = len(max_seg) 
			while max_seg[-1][0] == seg_rev[i][0]: 
				max_seg.append(seg_rev[i])
				i+=1
				if i == len(seg_cnts): 	return {a: b for a,b in seg_means},{a: b for a,b in seg_obs},(seg_min,seg_max),(1.0,1.0)
		min_len,max_len = len(min_seg),len(max_seg) 
		AAexp,ABexp  = min_len * (seg_lens[seg_min] / sum_lens), min_len * (seg_lens[seg_max] / sum_lens)
		BAexp,BBexp  = max_len * (seg_lens[seg_min] / sum_lens), max_len * (seg_lens[seg_max] / sum_lens)
		AAobs,ABobs = len([x for x in min_seg if x[1] == seg_min]), len([x for x in min_seg if x[1] == seg_max])
		BAobs,BBobs = len([x for x in max_seg if x[1] == seg_min]), len([x for x in max_seg if x[1] == seg_max])
		chi_low = chisquare([AAobs,ABobs],f_exp=[AAexp,ABexp])[1]
		chi_hi  = chisquare([BAobs,BBobs],f_exp=[BAexp,BBexp])[1]


		return {a: b for a,b in seg_means},{a: b for a,b in seg_obs},(seg_min,seg_max),(chi_low,chi_hi)






	def dex_score(self,res,Y,Xi,Xi_names):


		i_res = self.regress(Y,Xi,Xi_names)

		dex_key = {} 

		for c,(p,f) in res['params']['covariates'].items():	dex_key[c] = [(p,f)] 
		for pv,fc,n in res['params']['predictors']: 		dex_key[n] = [(pv,fc)]
		for pv,fc,n in i_res['params']['predictors']: 		dex_key[n].append((pv,fc))
	

		y = [log(y+1.0) for y in Y]
		seg_dict =  {k: [y[i] for i in V] for k,V in self.seg.items()}
		F_list = self.enrichment_and_fold_change(seg_dict)

		return {'params': dex_key, 'fcs': F_list}




		

























	def regress_ols(self,Y,X,X_names,key={},log_transform=True):

		req  = [] 
		if 'req' in key:   	    req    = key['req']



		p_out, z_out, r_out = dd(lambda: {}) ,{} , [] 
		#x_out, r_out, p_out, z_out = {}, [], dd(lambda: {}), {}
		model = sm.OLS([log(y+1.0) for y in Y],np.array(X)).fit() 
		for pv,bw,n in zip(model.pvalues,model.params,X_names):
			if (self.D.variable_class[n] == 'covariate') or (n == 'intercept'): 	p_out['covariates'][n] = (pv,bw)
			else:									r_out.append((pv,bw,n))
		p_out['predictors'] = sorted(r_out) 

	 	x_out = {'params': p_out, 'rs': model.rsquared, 'ars': model.rsquared_adj, 'bic': model.bic}



		if 'resids' in req:	x_out['resids'] = model.resid

		if 'pwr' in req:
			if model.rsquared < 0: 
				x_out['pwr-05'],x_out['pwr-001'] = 0.5,0.1
			else:
				f_2 =  model.rsquared / (1-model.rsquared) 
				df_de, df_num = len(X[0]) -1 , len(Y) - len(X[0]) 
				for alp,srn in zip([0.05,0.001],['pwr-05','pwr-001']): 
					x_out[srn] = smp.FTestPower().solve_power(effect_size=np.sqrt(f_2), df_num=df_num, df_denom=df_de, alpha=alp)
		return x_out



	def regress_glmnb(self,Y,X,interest=None):


		r_out, p_out, alp= {}, dd(lambda: {}), 0.05
		null = sm.GLM(Y, [np.array(1) for x in X], family=sm.families.NegativeBinomial()).fit()
		model = sm.GLM(Y, X, family=sm.families.NegativeBinomial()).fit()

		for p in self.D.inferred_predictors: 					p_out[p.split('=')[0]][p.split('=')[1]] = (1,0) 
		for pv,bw,c in zip(model.pvalues,model.params,self.D.predictors):	p_out[c.split('=')[0]][c.split('=')[-1]] = (pv,bw)
		for a,b in p_out.items():	r_out[a] = sorted(b.items(),key=lambda loc: loc[1][0])

		x_out = {'rs': 1 - (model.llf / null.llf), 'ars':  1 - ((model.llf-len(X[0])) / null.llf), 'bic': model.bic}

		f_2 =  x_out['rs'] / (1-x_out['rs'])
		df_de, df_num = len(X[0]) -1 , len(Y) - len(X[0]) 
		pwr = smp.FTestPower().solve_power(effect_size=np.sqrt(f_2), df_num=df_num, df_denom=df_de, alpha=alp)

		x_out['pwr'] = pwr 
		x_out['resids'] = [log(x+1.0) for x in model.resid_pearson]
		x_out['params'] = r_out
		return x_out

	

	def regress_zip(self,Y,X,interest=None):

		r_out, p_out, alp= {}, dd(lambda: {}), 0.05
		Y = np.array([np.array(log(y+1.0)) for y in Y]) 

		null = msc.PoissonZiGMLE(Y,np.array([1 for x in X])).fit(disp=0)
		model = msc.PoissonZiGMLE(Y,np.array(X)).fit(disp=0)
		params = model.params
		try: pvals = model.pvalues
		except ValueError: pvals = [0.99 for p in params]

		for p in self.D.inferred_predictors: 					p_out[p.split('=')[0]][p.split('=')[1]] = (1,0) 
		for pv,bw,c in zip(pvals,params,self.D.predictors):			p_out[c.split('=')[0]][c.split('=')[-1]] = (pv,bw)
		for a,b in p_out.items():	r_out[a] = sorted(b.items(),key=lambda loc: loc[1][0])

		x_out = {'rs': 1 - (model.llf / null.llf), 'ars':  1 - ((model.llf-len(X[0])) / null.llf), 'bic': model.bic}
		f_2 =  x_out['rs'] / (1-x_out['rs'])
		df_de, df_num = len(X[0]) -1 , len(Y) - len(X[0]) 
		pwr = smp.FTestPower().solve_power(effect_size=np.sqrt(f_2), df_num=df_num, df_denom=df_de, alpha=alp)
		x_out['resids'] = Y
		x_out['params'] = r_out 
		return x_out 










	def regress_poisson(self,Y,X,interest=None):
				## FIRST POISSON ## 

		model = sm.Poisson(Y,X).fit(disp=0)
#        	poisson_mod = sm.Poisson(my_vals, [1 for v in my_vals])
 #               poisson_res = poisson_mod.fit(method="newton",disp=0)
#		poisson_pv =  poisson_res.pvalues[0]
#		pAIC,pBIC = poisson_res.aic, poisson_res.bic




	def regress_nb(self,Y,X,interest=None):


		print 'uh'

		#sm.GLM(data.endog, data.exog, family=sm.families.Gamma())
		foo = sm.GLM(Y, X, family=sm.families.NegativeBinomial(),variance=10).fit()
		foo = sm.GLM(Y, X, family=sm.families.NegativeBinomial()).fit()



		print foo.summary()


		model = sm.NegativeBinomial([log(y+1.0) for y in Y],X).fit()
#		print model.summary()
		sys.exit() 
	
		print len(model.pvalues)
		print len(model.params)
	
		print model.bic
#		print model.rsquared
#		print model.rsquared_adj

		for v in vars(model._results): print v 
#		print model.summary() 

		for v in vars(model.model):
			print v
		

		
#		p_out = {'params': r_out, 'bic': model.bic, 'rs': model.rsquared, 'ars': model.rsquared_adj, 'resids': model.resid, 'pwr': pwr}
		


		mPV,aPV = res_nbin.pvalues
		nbM,nbA = exp(res_nbin.params[0]),res_nbin.params[1] 
		estX,estP = convert_nb(nbM,nbA)
		my_comps = stats.nbinom.rvs(estX, estP, size=len(my_vals))
		chiT,chiP = self.bin_chi(my_vals,my_comps,min(binInt,int(len(my_vals)*binRate)))			
		self.tests['nbin'] = (chiT,chiP) 

		nbAIC,nbBIC = res_nbin.aic, res_nbin.bic

		print m.name,len(vals),len(dZ),val_type,'neg-binom',chiT,chiP,"|",mPV,aPV,'|', nbAIC,nbBIC


		sys.exit() 	






































	def calculate_ss(self,features,params,interest):
		f_key = dd(lambda: dd(float))
		headers = self.options.featureKey.readline().split() 
		for line in self.options.featureKey:
			line = line.split()
			for i in range(1,len(line)): 
				if float(line[i]) > 0: 
					if headers[i].split('=')[0] == interest: f_key[line[0]][headers[i].split('=')[-1]] = float(line[i])
		pvs = [0.05, 0.005, 0.0005, 0.00005] 
		ss_key = dd(lambda: dd(int))
		for i,f in enumerate(features):
			for pv in pvs:
				tp,fp,tn,fn = 0,0,0,0
				for (p,(a,b)) in params[i][interest]:
					if a < pv and f_key[f.name][p] > 0.0:	 tp+=1
					elif a < pv and f_key[f.name][p] == 0.0: fp+=1
					elif a > pv and f_key[f.name][p] != 0.0: fn+=1 
					else:					 tn+=1 
				if tp > 0: ss_key[pv]['TP'] += 1
				elif fp > 0: ss_key[pv]['FP'] += 1
				elif fn > 0: ss_key[pv]['FN'] += 1 
				else:        ss_key[pv]['TN'] += 1
		for pv in ss_key:
			ss_key[pv]['se'] =  float(ss_key[pv]['TP']+0.001) / (0.0001+ss_key[pv]['TP']+ss_key[pv]['FN'])
			ss_key[pv]['sp'] =  float(ss_key[pv]['TN']+0.001) / (0.00001+ss_key[pv]['TN']+ss_key[pv]['FP'])
		return ss_key 	



































































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
from scipy.stats import spearmanr as spearmanr
from scipy.stats import poisson as PSN 
#from scipy.stats import spearmanr as spearmanr
import seaborn
from math import log
import math
from random import shuffle
from sklearn.cluster import KMeans	
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import MinMaxScaler
import matplotlib
import matplotlib.pyplot as plt

from operator import mul
import scipy
import statsmodels.api as sm

from statsmodels.miscmodels import count  as msc
import pandas as pd 
from scipy.signal import savgol_filter as svgf 
from math import exp
from math import factorial 





def scale_vals(vals):
        scaler = MinMaxScaler()
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))

def poisson_approximate(s1,s2,s_both,sLen):

	if s1 >= s2: 
		np = s2*(s1/float(sLen))
	else:
		np = s1*(s2/float(sLen))
	if s_both == 0.0: return round(np,3),PSN.pmf(0,np) 
	elif s_both < np: return round(np,3),PSN.cdf(s_both,np)
	else:		  return round(np,3),1-PSN.cdf(s_both,np)
#a3 =  PSN.pmf(0,0.5)
#print PSN.cdf(2,0.5),'cdf'


def convert_nb(mu,theta):

#	var = mu + ((mu*mu)*alp)
	var = mu + theta * mu ** 2
	p = (var - mu) / var
    	return 1.0/theta, 1 - p
	

#	p = (var - mu) / float(var) 

#	r = (mu*mu) / (var - mu) 

#	return p,r 



class Dists:
	
        def __init__(self,args,members,member_space,progress=None):


		self.args, self.members, self.space, self.progress = args, members, member_space, progress 



	def bin_chi(self,vals,r,maxSize): 

		both = sorted(cc(vals+r).items()) 					
					
		bins, span, sT = [],[0], 0   	
		for v,c in both: 
			span.append(v) 
			sT += c 
			if sT >= maxSize: 
				bins.append((span[0],span[-1]))
				span,sT = [v+1],0
		if span[0] != span[-1]: bins.append((span[0],span[-1]))
		
		n = 0 
		rC,rK,rCnts = sorted(cc(r).items()),0,[1 for b in bins]
		vC,vK,vCnts = sorted(cc(vals).items()),0,[0 for b in bins]

		while n < len(bins): 

			while rK < len(rC) and rC[rK][0] < bins[n][0]: rK+=1
			while rK < len(rC) and rC[rK][0] <= bins[n][1]: 
				rCnts[n] += rC[rK][1]	
				rK+=1 
			
			while vK < len(vC) and vC[vK][0] < bins[n][0]: vK+=1
			while vK < len(vC) and vC[vK][0] <= bins[n][1]: 
				vCnts[n] += vC[vK][1]	
				vK+=1 

			n+=1				
		chiT,chiP= stats.chisquare(vCnts, f_exp=rCnts)
		return round(chiT,4),chiP







	def fit_binary(self,minSize=20,binInt=5,binRate=0.2):


		bin_dists = ['poisson','nbinom']
		bin_dists = ['poisson']

		for m in self.members:  
			vals,logV, dZ = [int(x) for x in m.cnts.values()], [log(v+1.0) for v in m.cnts.values()],[0 for i in range(self.space - len(m.cnts.values()))]
			if len(vals) < minSize: continue 
				
			val_key = {'RAW-NZ': vals, 'RAW-WZ': vals+dZ,  'LOG-NZ': logV, 'LOG-WZ':  logV+dZ}
			
			for val_type,my_vals in val_key.items(): 
				self.tests = {} 	
				vLen, bR, vMean = len(my_vals), int(len(my_vals)*binRate)  , np.mean(my_vals) 
				if val_type.split('-')[0] == 'LOG': continue 

				## FIRST POISSON ## 

                             	poisson_mod = sm.Poisson(my_vals, [1 for v in my_vals])
                                poisson_res = poisson_mod.fit(method="newton",disp=0)
				poisson_pv =  poisson_res.pvalues[0]
				pAIC,pBIC = poisson_res.aic, poisson_res.bic

				poisson_sample     = stats.poisson.rvs(vMean, size=len(my_vals))
				chiT,chiP = self.bin_chi(my_vals,poisson_sample,min(binInt,int(len(vals)*binRate)))
				self.tests['poisson'] = (chiT,chiP) 
				print m.name,len(vals),len(dZ),val_type,'poisson',chiT,chiP,'|',poisson_pv, 'NA','|',pAIC,pBIC

				## NEGATIVE BINOMIAL ## 
				
				mod_nbin = sm.NegativeBinomial(my_vals, [1 for v in my_vals])
				res_nbin = mod_nbin.fit(disp=0)
				
				mPV,aPV = res_nbin.pvalues
				nbM,nbA = exp(res_nbin.params[0]),res_nbin.params[1] 
				estX,estP = convert_nb(nbM,nbA)
				my_comps = stats.nbinom.rvs(estX, estP, size=len(my_vals))
				chiT,chiP = self.bin_chi(my_vals,my_comps,min(binInt,int(len(my_vals)*binRate)))			
				self.tests['nbin'] = (chiT,chiP) 

				nbAIC,nbBIC = res_nbin.aic, res_nbin.bic

				print m.name,len(vals),len(dZ),val_type,'neg-binom',chiT,chiP,"|",mPV,aPV,'|', nbAIC,nbBIC

				## NOW ZERO P ###
				if val_type.split('-')[-1] == 'NZ': continue
				zp_nbin = msc.PoissonZiGMLE(my_vals, [1 for v in my_vals])
				res_zp = zp_nbin.fit(disp=0)
				zpAIC,zpBIC = res_zp.aic, res_zp.bic
				zpM = exp(res_zp.params[0])
				zpZ  =  1 - (np.mean(my_vals) / zpM)

				try: 
					cPV,zPV = res_zp.pvalues
				except ValueError:
					cPV,zPV = 'NA','NA'
					print 'hmmm' 
				my_comps = [x if random.random() > zpZ else 0 for x in stats.poisson.rvs(zpM, size=len(my_vals))]
				chiT,chiP = self.bin_chi(my_vals,my_comps,min(binInt,int(len(my_vals)*binRate)))			
				self.tests['zp'] = (chiT,chiP) 
				print m.name,len(vals),len(dZ),val_type,'zip-po',chiT,chiP,"|",cPV,zPV,'|', zpAIC,zpBIC

		sys.exit() 	






































	def fit_dists(self):
		CUTOFF = 50
		anderson_dists = ['norm','expon','logistic','gumbel','extreme1']
		print '---','obs','zeros','datatype','mean','cv','|','dist','test','ts','pv'	
 
		for f in self.input.features: 

			if len(f.cnts.values()) < CUTOFF: continue 
			my_len = len(f.cnts.values())
			z_len = [0 for i in range(self.input.samples.len-my_len)]
			my_vals = [(f.cnts.values(),'raw-nonzero')]
			my_vals.append((my_vals[0][0]+z_len,'raw-withzero'))
			my_vals.append(([log(v+1.0) for v in my_vals[0][0]],'log-nonzero'))
			my_vals.append((my_vals[2][0]+z_len,'log-withzero'))


			for vals,val_type in my_vals:

				for a in anderson_dists: 
					at,cv,sl =  stats.anderson(vals,a)		
					sig = sl[0] 
					for x,y in zip(cv[-1::-1],sl[-1::-1]):
						if at > x: 
							sig = y 
						break 
					print f.name,my_len,len(z_len),val_type,np.mean(vals),stats.variation(vals),'|',a,'ANDERSON',round(at,4),sig

				for a in ['norm','beta','gamma','wald','t','lognorm','halflogistic']:
        				dist = getattr(scipy.stats, a)
        				param = dist.fit(vals)
        				try: 
						gf = stats.kstest(vals, a, param)
						print f.name,my_len,len(z_len),val_type,np.mean(vals),stats.variation(vals),'|',a,'KS',round(gf[0],4),gf[1]
					except ValueError: 
						print f.name,my_len,len(z_len),val_type,np.mean(vals),stats.variation(vals),'|',a,'KS','FAIL','NA'

		sys.exit() 


	def summarize_relationships(self):
		CUTOFF = 50 
		CUTOFF = 50
		f_sets = {} 
		sLen = len(self.input.samples)

		f_logs = {f: {s.idx: 0 if f.idx not in s.cnts else log(1.0+s.cnts[f.idx]) for s in self.input.samples} for f in self.input.features}
 
		for i,f in enumerate(self.input.features):

			if len(f.cnts.keys()) < CUTOFF: continue 

			f_sets[f] =   set(f.cnts.keys())
			#  set([  s for s in self.input.feature_vals[f].keys()])

		f_keys = f_sets.keys()
 
		print '--- f2 kind len1 len2 | intersect lambda pv | pR rV rS sV'
		for i in range(len(f_keys)-1):
			f1 = f_keys[i] 
			f1s = f_sets[f1] 
			for j in range(i+1,len(f_keys)):
				f2 = f_keys[j]
				f2s = f_sets[f2] 
				s_inter =  set.intersection(f1s,f2s)
				s_both =  len(s_inter)
				ld,pf = poisson_approximate(len(f1s),len(f2s),s_both,self.input.samples.len)
				if pf < 0.05 and s_both < ld:
					z=5
					print f1.name,f2.name,'NEG',len(f1s),len(f2s),'|',s_both,ld,pf
				elif s_both < ld or s_both < 10: continue 

				elif pf < 0.05 or s_both > 1000:
					


					R,rP = pearsonr([f_logs[f1][s] for s in s_inter],[f_logs[f2][s] for s in s_inter])						
					S,sP =spearmanr([f_logs[f1][s] for s in s_inter],[f_logs[f2][s] for s in s_inter])						
					print f1.name,f2.name,'POS',len(f1s),len(f2s),'|',s_both,ld,pf,'|',round(R,3),rP,round(S,3),sP

				else:
					continue 































	def summarize_sample_stats(self):

		self.progress.start_subtopic('Calculating Summary Stats','',self.sLen)
		res = dd(lambda: {}) 
		subplot = rage_subplots.subplot(3,1,self.args)  
		for s in self.input.samples: 	
			self.progress.mark_subtopic() 
			ordered_logs = sorted([log(1.0+self.input.sample_vals[s][f]) for f in self.input.sample_vals[s]],reverse=True)
			halfE,iX,k = sum(ordered_logs)*0.5,0,-1
			res['#Observed_Genes'][s] = len(ordered_logs) 
			res['#Genes_Above_Mean'][s] = len([x for x in ordered_logs if x > np.mean(ordered_logs)])/float(len(ordered_logs))
			while iX < halfE:
				k+=1;	iX+=ordered_logs[k] 
			res['%Genes_Required_For_HalfDepth'][s] = k / float(len(ordered_logs))

		subplot.add_hist(res['#Observed_Genes'].values()).update({'xlab':'genes per sample','ylab': 'occurences','title': 'Library Complexity'})	
		subplot.add_hist(res['#Genes_Above_Mean'].values()).update({'xlab':'%','ylab': 'occurences','title': '% genes above mean'})
		subplot.add_hist(res['%Genes_Required_For_HalfDepth'].values()).update({'xlab':'%Obs Genes','ylab': 'occurences','title': '% Genes Required For 50% Read Depth (Log Space)'})
		subplot.save('sample_summary.png',self.args) 
		
		rage_outputs.column_stats(self.args).write(res,self.input.samples,{'suffix': 'sample_stats.out','width': 20})

		self.progress.finish_subtopic() 

	def create_label_key(self):
		

		self.color_key = {'EB': 'lime', 'T': 'red', 'ES': 'cyan', 'O': 'grey','U': 'purple','H': 'orange'}

		try: 
			self.color_labels = [self.color_key[s.name[0]] if s.name[0] != 'E' else self.color_key[s.name[0:2]] for s in self.input.samples] 
		except KeyError:
			self.color_labels = ['k' for s in self.input.samples]

	def make_pca_and_tsne_plots(self):

		self.progress.start_subtopic('Calculating PCA/TSNE','',0)


		data_matrix = self.input.data_matrix('log')
		dr = rage_DR.DR(self.args,self.progress).set_matrix(data_matrix)
		dr.run_pca().run_kca().run_tsne()
#		log_vals = [[log(x+1.0) for x in c] for c in self.input.feature_cnts]

#		dr = rage_DR.DR(self.args,self).run_pca().run_tsne() 
#		dr. = rage_DR.DR(self.args,self).run_pca().run_kca().run_tsne() 



		subplot = rage_subplots.subplot(1,2,self.args)
		subplot = rage_subplots.subplot(1,3,self.args)
		#color_key = {'EB': 'lime', 'T': 'red', 'ES': 'cyan', 'O': 'grey','U': 'purple','H': 'orange'}
		#labels = [color_key[s[0]] if s[0] != 'E' else color_key[s[0:2]] for s in self.input.samples] 
		subplot.add_pca_data(dr.pca_pts,{'vars': dr.pca_vars,'title': 'PCA','colors':self.color_labels}).update() 
		subplot.add_pca_data(dr.kca_pts,{'type': 'kca', 'title': 'KCA','colors':self.color_labels}).update() 
		subplot.add_pca_data(dr.tsne_pts,{'type': 'tsne','colors':self.color_labels}).update() 
		subplot.add_legend(self.color_key.keys(),self.color_key.values())
		subplot.save('dim_red.png',self.args) 
		self.progress.finish_subtopic() 


	def get_feature_order(self):

		feature_ranks = dd(list) 
		for s in self.input.samples: 
			for i,(b,a) in enumerate(sorted([(b,a) for (a,b) in s.cnts.items()])):
				if i == 0: match,rank,m_list = b,1,[a]
				elif b == match: m_list.append(a) 
				else:
					for m in m_list: feature_ranks[m].append(rank) 
					match,rank,m_list = b,rank+1,[a]
			for m in m_list: feature_ranks[m].append(rank) 
		self.feature_order = [x[1] for x in sorted([(sum(b),a) for (a,b) in feature_ranks.items()])]







	def summarize_sample_pts(self,pt_label='val'):

		self.progress.start_subtopic('Plotting All Pts','',self.sLen)

		subplot = rage_subplots.subplot(1,1,self.args)  
		HUMAN_COLORS=False 
#		HUMAN_COLORS=True	
		for s in self.input.samples:
			self.progress.mark_subtopic() 

			if s[0] in ['U','H']: continue 
			ordered_vals = [self.input.sample_vals[s][f] if f in self.input.sample_vals[s] else 0 for f in self.feature_order]
			ordered_logs = [log(x+1.0) for x in ordered_vals]
			scaled_logs = scale_vals(ordered_logs) 
			if self.args.organism == 'human' and HUMAN_COLORS:
				if s[0] == 'H': continue 
				elif s[0] == 'U': continue 
				elif s[0] == 'T':    s_color = 'red' 
				elif s[0:2] == 'EB': s_color = 'lime'  
				elif s[0:2] == 'ES':  s_color = 'blue' 
				elif s[0:2] == 'OB':  s_color = 'gray' 
				else:		      s_color = 'black'

			XY = [(x,scaled_logs[x]) for x in range(len(scaled_logs))]
			color_groups, group_colors = [[xy for xy in XY if xy[1] == 0]], [0]
			for (a,b,c) in [(d/20.0,(d+1)/20.0,(d+d+1.0)/40.0) for d in range(0,20)]:
				color_groups.append([xy for xy in XY if xy[1] > a and xy[1] <= b])
				group_colors.append(c) 
			diff_colors = get_colors(group_colors, plt.cm.jet) 
			for g,grp in enumerate(color_groups):
				if len(grp) == 0: 		continue 
				if HUMAN_COLORS: 		clr = s_color 
				elif grp[0][1] == 0.0: 		clr = 'k'
				else: 				clr = diff_colors[g] 

				if grp[0][1] == 0.0: subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'color': clr, 'size': 0.1, 'alpha': 0.3,'yjitter': True})  
				elif grp[0][1] < 0.25: subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'color': clr, 'size': 0.1, 'alpha': 0.3,'yjitter': True})  
				elif grp[0][1] < 0.50: subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'color': clr, 'size': 0.1, 'alpha': 0.3,'yjitter': True})  
				elif grp[0][1] < 0.80: subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'color': clr, 'size': 0.1, 'alpha': 0.3,'yjitter': True})  
				elif grp[0][1] < 0.95: subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'color': clr, 'size': 0.1, 'alpha': 0.3,'yjitter': True})
  				else:	             subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'color': clr,'size': 0.1,'alpha':0.3,'yjitter': True})  
		edge = int(len(XY)*0.05)
		subplot.ax.set_xlim(0-edge,len(XY)+edge)
		subplot.ax.text(len(XY)/2.5,-0.1,'Ordered Genes',fontweight='bold',fontsize=15)
		#subplot.ax.text(len(XY)/2.6,1.1,'2k Mouse Cell Expression',fontweight='bold',fontsize=15)
		subplot.ax.text(len(XY)/2.6,1.1,'2k Single Cell Expression',fontweight='bold',fontsize=15)

		if HUMAN_COLORS: 
			subplot.add_legend(['EB','ES','OB','T'],['lime','blue','gray','red'])

		#plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.95,wspace=0.01,hspace=0.0)
		subplot.save('all_pts_dist.png',{'axis_off': True}) 
		self.progress.finish_subtopic() 
		


		


	def summarize_sample_dists(self):

		self.progress.start_subtopic('Plotting Sample Densities','',self.sLen)
		kde = rage_KDE.samples(0.3) 
		subplot,f_num = rage_subplots.subplot(10,10,self.args), 1 
		LOG=True

		s_id = ''
		my_class = 'CLASS_1' 
		my_class = 'CLASS_2' 
		my_class = 'CLASS_3' 
		for s in self.input.samples:
			self.progress.mark_subtopic() 	
			if len(self.input.sample_key['CLASS'].keys()) > 0: 
				if self.input.sample_key['CLASS'][s] != my_class: continue 

			sample_features = self.input.sample_vals[s].keys()  
			sample_items = self.input.sample_vals[s].items() 

			if LOG: 
				non_zeros = [log(x+1.0) for x in sorted([b for (a,b) in sample_items])]
				all_vals = [0 for x in range(self.fLen - len(sample_items))] + non_zeros
			else:
				non_zeros = [x for x in sorted([b for (a,b) in sample_items])]
				all_vals = [0 for x in range(self.fLen - len(sample_items))] + non_zeros

			x1,y1 = kde.run(all_vals)
			x2,y2 = kde.run(non_zeros)

			subplot.add_lines(x1,y1,None,None,'black')
			subplot.add_lines(x2,y2,None,None,'green')
			subplot.change_limits({'x0': -0.5,'x1': 8, 'y0': -0.1,'y1':1.5}) 
			subplot.ax.text(1.4,0.91,s+' ( '+str(len(non_zeros))+' )',color='blue')
			subplot.ax.set_xticklabels([]) 
			subplot.ax.set_yticklabels([]) 
			
			if not subplot.update({'clear_axes': True}): 
				plt.suptitle('Dual Dists: '+my_class) 
				plt.subplots_adjust(left=0.04, bottom=0.01, right=0.96, top=0.95,wspace=0.03,hspace=0.04)
				#fig.savefig('fig_'+my_class+"_"+str(f_num)+'.png', dpi=100)	
				subplot.save('fig_'+my_class+"_"+str(f_num)+'.png',{'title': 'Dual Dists: '+my_class})
				f_num += 1
				subplot = rage_subplots.subplot(10,10,self.args)  
				#break 
		plt.subplots_adjust(left=0.02, bottom=0.01, right=0.98, top=0.95,wspace=0.03,hspace=0.03)
		subplot.save('fig_'+my_class+"_"+str(f_num)+'.png',{'title': 'Dual Dists: '+my_class})

		self.progress.finish_subtopic() 

		



	def summarize_sample_pairs(self):

		from modules.Rage_Plots import rage_subplots

		xLen,yLen = 5,5
		subplot = rage_subplots.subplot(xLen,yLen,True)  
		total_features = len(self.input.features) 
		f_num = 1
		LOG=True
		feature_sample_ranks = dd(lambda: dd(float))
		for s in self.input.samples:
			for i,(b,a) in enumerate(sorted([(b,a) for (a,b) in self.input.sample_vals[s].items()])):
				if i == 0: match,rank,m_list = b,1,[a]
				elif b == match: m_list.append(a) 
				else:
					for m in m_list: feature_sample_ranks[s][m] = rank
					match,rank,m_list = b,rank+1,[a]
			feature_sample_ranks
			for m in m_list: feature_sample_ranks[s][m] =  rank
		f_num = 1 
		fig = matplotlib.pyplot.gcf()
		fig.set_size_inches(18.5, 9.5)
		s_id = ''
		for i in range(len(self.input.samples)):
			for j in range(i+1,len(self.input.samples)):

				s1,s2 = self.input.samples[i],self.input.samples[j]
				fr1,fr2 = feature_sample_ranks[s1],feature_sample_ranks[s2]
				fkeys = list(set(fr1.keys()+fr2.keys()))
				f_order = [x[1] for x in sorted([(fr1[f]+fr2[f],f) for f in fkeys])]
				x_range = range(len(f_order))
				v1  = [log(1.0+self.input.sample_vals[s1][f]) if f in self.input.sample_vals[s1] else 0 for f in f_order]
				v2  = [log(1.0+self.input.sample_vals[s2][f]) if f in self.input.sample_vals[s2] else 0 for f in f_order]

				vs1 = scale_vals(v1) 
				vs2 = scale_vals(v2) 
 				sv1 = svgf(vs1, 61, 2, mode='nearest')
 				sv2 = svgf(vs2, 61, 2, mode='nearest')
				subplot.add_line(x_range,sv1,{'lw': 0.2})
				subplot.add_line(x_range,sv2,{'lw': 0.2})
				sv_mix = [(sv1[z]+sv2[z])/2.0 for z in range(len(sv1))]

				step1,step2 = 50,100 
				subplot.add_line(x_range,sv_mix,{'lw': 0.5,'color':'k'}) 
				z_diffs,z_steps = [], [] 
				for z in range(step2,len(sv_mix),step1):
					z1 = sv1[z-step2:z+step2]
					z2 = sv2[z-step2:z+step2]
					z_diffs.append(sum([(z1[x]-z2[x])*(z1[x]-z2[x]) for x in range(len(z1))]))
					z_steps.append((z-step2,z+step2))

				
					#subplot.add_line(x_range[z-step2:z+step2],sv_mix[z-step2:z+step2],{'color': 'purple','alpha':0.4})
				diff_colors = get_colors(z_diffs, plt.cm.jet)
				for z in range(len(z_steps)):
					zA,zB = z_steps[z]		
					subplot.add_line(x_range[zA:zB],sv_mix[zA:zB],{'color': diff_colors[z],'alpha':0.5,'lw': 1})
 
		
				#subplot.change_limits({'x1': int(len(x_range)*1.08), 'y0': -0.05,'y1': 0.93}) 
				subplot.ax.text(int(len(x_range)*0.03),0.72,s1+'  '+s_id+' '+s2+' '+s_id,color='red')
				#subplot.ax.plot([0,len(x_range)],[0,0],color='k',linewidth=1,zorder=2) 
				if not subplot.update(): 
					plt.suptitle('Pair Comparison') 
					plt.subplots_adjust(left=0.04, bottom=0.01, right=0.96, top=0.95,wspace=0.03,hspace=0.03)
					fig.savefig('pairs_out'+str(f_num)+'.png', dpi=100)	
					f_num += 1
					if f_num > 10: sys.exit() 

				
		sys.exit()














































