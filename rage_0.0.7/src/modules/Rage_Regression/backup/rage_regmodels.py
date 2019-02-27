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
	def __init__(self,options,X,FULL = False):

		self.options, self.X, self.vif  = options, X, vif_test(X) 
		self.alp1, self.alp2 = 0.05, 0.001	
		
		if options.model == 'OLS':   self.regress = RegOLS(X,alphas= [self.alp1, self.alp2]) 
	

		self.out = dd(list) 	
		self.rs_key = [0.01,  0.03,   0.05,    0.10,     0.25]
		self.pv_key = [0.05, 0.01, 0.001, 0.0001, 0.00001,0.0000001]
		# self.pv_key = [0.01, 0.001, 0.0001, 0.00001,0.0000001]
		self.resid = [] 

	def run(self, Y):
		for y in Y: 
			t = self.regress.test(y) 
			self.out['params'].append(t.output) 
			if len(self.X.predictor_idx)>0:	 self.out['pvs'].append([p[0] for p in sorted(t.output) if p[-1]][0]) 
			for r,k in zip([t.rsq, t.rsa, t.bic, t.resids, t.c_resids,t.pwr[self.alp1],t.pwr[self.alp2]],['rsq','rsa','bic','resids','c_resids','pwr1','pwr2']): 
				self.out[k].append(r) 		

		return self	


	def aggregate(self,PARAMS=False): 

		if len(self.X.predictor_idx) > 0:	pv_list = [sorted([p[0] for p in P if p[-1]])[0] for P in self.out['params']]
		else:					pv_list = [] 

		self.pv_cnt = [len([p for p in pv_list if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
		self.rs_cnt = [len([p for p in self.out['rsq'] if p > self.rs_key[j]]) for j in range(len(self.rs_key))]
		
		if PARAMS: 	
			self.pv_dict = {self.X.names[i] : [p[i][0] for p in self.out['params']] for i in range(len(self.X.names))} 

		return self 		



	def run_permutation_tests(self,Y,V,minT=199,midT=1000,maxT=10000): 

		minCut, midCut  = 0.05 * minT, 0.01 * (midT+minT)


#		Xp = [ V.select_variables(V.variables, permute=V.predictors) for j in range(minT)]

#		self.permutations = {n: [1.0 for z in range(len(Y))] for i,n in enumerate(self.X.names) if i in self.X.predictor_idx}
		self.permutations = dd(list) 
#		return
		# self.pv_dict 
		for i,y in enumerate(Y):


			print "RUNNING GENE",i+1, 
			p_cands,self.p_key = [[a[2],a[0]] for a in [r for r in self.out['params'][i] if r[-1]]],{} 

			#print "" 	
			#print "" 	
			#print "INITIAL PREDS", p_cands

			for c,pv in p_cands: 			
				if pv > 0.10: self.permutations[c].append(round(pv,2)) 
				else: 	      self.p_key[c] = [pv,0,0]
		
			
			
			if len(self.p_key) > 0: 

				print "ROUND 1 TRIGGERED",

				self.regress.permute(y,[ V.select_variables(V.variables, permute=V.predictors) for j in range(minT)],self.p_key) 

				for c,(pv,L,G) in self.p_key.items(): 
					#if L > 3: iself.permutations[c].append((self.p_key.pop(c)[1])/float(minT))
					if L > minCut:    self.permutations[c].append((self.p_key.pop(c)[1])/float(minT))



				if len(self.p_key) > 0: 	
					print "ROUND 2 TRIGGERED",
					self.regress.permute(y,[ V.select_variables(V.variables, permute=V.predictors) for j in range(midT)],self.p_key) 
					for c,(pv,L,G) in self.p_key.items(): 
						if L> midCut:    self.permutations[c].append((self.p_key.pop(c)[1])/float(minT+midT))

					if len(self.p_key) > 0: 	
						print "ROUND 3 TRIGGERED",
						self.regress.permute(y,[ V.select_variables(V.variables, permute=V.predictors) for j in range(maxT)],self.p_key) 
						for c,(pv,L,G) in self.p_key.items(): 
							self.permutations[c].append((self.p_key.pop(c)[1])/float(minT+midT+maxT))

			print ""
			print "RESULTS:",
			for c,cI in self.permutations.items(): print c,cI[-1],
			print ""
			if i > 10: sys.exit()  



	def pv_discovery2(self,pvs,disc_pvs = [0.05,0.01,0.001,0.0001]):

		if len(self.pv_mins) != len(pvs):
			print 'diff lens wtf'
			sys.exit() 

#		disc_dict = dd(lambda: [0,0])
		disc_dict = {dpv: 0 for dpv in disc_pvs} 
		for p,s in zip(self.pv_mins,pvs):
			if p > disc_pvs[0] and s > disc_pvs[0]: continue 
			for i in range(len(disc_pvs)):
				if len([x for x in [p,s] if x<disc_pvs[i]]) ==1:
					if p < disc_pvs[0]: disc_dict[disc_pvs[i]] -=1
					else:		    disc_dict[disc_pvs[i]] +=1



				

		return disc_dict


	### EACH TIME IT JUST GETS A TWO DIMENSIONAL X ###


	### BUT WAIT HOW ABOUT THINGS LIKE BALANCE ???  ### 

	### BALANCE AND VIF ARE DONE X vs X ### 

	### CHI SQUARE FOR SEG IS NOT ### --- THAT IS DONE ON THE GENE LEVEL ### 





class RegOLS:
	def __init__(self,X,alphas=[0.05,0.001]):

		self.X, self.alphas = X, alphas 
		
		self.dfd  =  len(self.X.names) -1 
		self.dfn =   len(self.X.array) - len(self.X.names)  


	
	def test(self,y):


		model = sm.OLS(y,self.X.array).fit()

		self.rsq, self.rsa, self.bic  = round(model.rsquared,5), round(model.rsquared_adj,3), round(model.bic,3)
		self.output = [(p,b,x,i in self.X.predictor_idx) for i,(p,b,x) in enumerate(zip(model.pvalues,model.params,self.X.names))]


		try : self.pwr = {a: smp.FTestPower().solve_power(effect_size=np.sqrt(self.rsq/(1-self.rsq)),df_num=self.dfn,df_denom=self.dfd,alpha=a) for a in self.alphas}
		except: self.pwr = {a: 0.5 for a in self.alphas}
		self.resids, self.c_resids = model.resid, [sum([x[j] * model.params[j] for j in self.X.covariate_idx]) + y[i] for i,x in enumerate(self.X.array)]
		

		return self


	def permute(self,y,Xp,key):  # y = counts for a feature, Xp = list of permutated variable arrays, key = G/L counts for permuatations  #  		
		# key == dictionary -> predictor :  [ model-predictor-pvalue , L = number of observations less than, G = number of observations greater than ] 
		for P in Xp: 
			model = sm.OLS(y,P.array).fit()
			for i,(p,b,n) in enumerate(zip(model.pvalues,model.params,self.X.names)):
				if n in key:
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1




def regress_ols2(y,X,req=[]):

	model,f_out = sm.OLS(y,X.array).fit(), dd(list)


	## PARAMS ## 
	for pv,bw,n in zip(model.pvalues,model.params,X.names):	f_out[X.parent[n]].append((pv,bw,n)) 
	x_out = {'params': {n: sorted(p) for n,p in f_out.items()}, 'rs': model.rsquared, 'ars': model.rsquared_adj, 'bic': model.bic, 'pwr-05': 0.5, 'pwr-001': 0.01} 

	if 'resids' in req:	x_out['resids'] = model.resid
	if 'pwr' in req and model.rsquared > 0:
		df_de, df_num, f_2 = len(X.names) -1 , len(y) - len(X.names) , np.sqrt(model.rsquared / (1-model.rsquared) )
		x_out['pwr-05'],x_out['pwr-001'] = smp.FTestPower().solve_power(effect_size=f_2,df_num=df_num,df_denom=df_de,alpha=0.05),smp.FTestPower().solve_power(effect_size=f_2, df_num=df_num, df_denom=df_de, alpha=0.001)

	if 'predictors-only' in req and len(X.p_names)>1:

		p_model,p_out = sm.OLS(y,X.p_array).fit(), dd(list) 
		
		for pv,bw,n in zip(p_model.pvalues,p_model.params,X.p_names):	p_out[X.parent[n]].append((pv,bw,n)) 
		P_out = {'params': {n: sorted(p) for n,p in p_out.items()}, 'rs': p_model.rsquared, 'ars': p_model.rsquared_adj, 'bic': p_model.bic, 'resids': model.resid}
		x_out['predictors'] = P_out
		
		



	if 'covariates-only' in req and len(X.c_names)>1:

		c_model,c_out = sm.OLS(y,X.c_array).fit(), dd(list) 
		for pv,bw,n in zip(c_model.pvalues,c_model.params,X.c_names):	c_out[X.parent[n]].append((pv,bw,n)) 
		C_out = {'params': {n: sorted(p) for n,p in c_out.items()}, 'rs': c_model.rsquared, 'ars': c_model.rsquared_adj, 'bic': c_model.bic, 'resids': model.resid}
		x_out['covariates'] = C_out 
		
		
	return x_out 


























def list_stats(L):

	xSum,xLen,xObs = sum(L) , float(len(L)), len([x for x in L if x>0])
	avg,obs = xSum/xLen, xObs/xLen
	std = np.std(L) 
	cv = coVar(L) 
	return int(xLen),xObs,round(avg,3),round(obs,3),round(std,3),round(cv,3) 


def intersect(LOL):
	return [x for (x,y) in cc([a for b in LOL for a in b]).items() if y == len(LOL)]

def offset_list(my_list):
	if min(my_list) > 0: return my_list
	else:		     return [m-min(my_list) for m in my_list]


def sva_resids(y,X_interest,X_sva,TYPE='OLS'):	
	yResid =   offset_list(sm.OLS(y, X_interest).fit().resid) 
	yNew   =   offset_list(sm.OLS(yResid, X_sva).fit().resid)
	ySVA   =   offset_list([y[j] + (yNew[j] - yResid[j]) for j in range(len(y))])	
	return ySVA







def ttest(opt_cnts,opts):
	summary_data,tests,summary_key = [], [],{}
	for i in range(len(opts)): 
		A = opt_cnts[i]	
		if len(opt_cnts[i]) > 3:  		
			summary_data.append([opts[i],opt_cnts[i],sum(opt_cnts[i])/float(len(opt_cnts[i])), len([a for a in opt_cnts[i] if a>0]),len(opt_cnts[i])])
	for i in range(len(summary_data)):
		A_name,A_cnts ,A_mean, A_obs,A_len = summary_data[i] 
		summary_key[A_name] = [A_mean,A_obs,A_len,A_cnts]

		for j in range(i+1,len(summary_data)):
			B_name,B_cnts ,B_mean, B_obs,B_len = summary_data[j] 
			pv = stats.ttest_ind(A_cnts,B_cnts)[-1]
			if pv > 0.001: pv = round(pv,4) 
			if A_mean > B_mean: FC = round(A_mean/(B_mean+0.001),3)
			else:		    FC = -1*round(B_mean/(A_mean+0.001),3)		
			tests.append([(A_name,B_name),pv,FC])
	return [tests,summary_key]
			











#def regress(y,X,INTEREST=None,TYPE='OLS'): 
def regress(y,X,TYPE='OLS'):



	if TYPE == 'OLS': model = sm.OLS(y,X).fit()
#	pvs    = {model.pvalues._index[i]: model.pvalues[i] for i in range(len(model.pvalues))}
#	params = {model.params._index[i]: model.params[i]   for i in range(len(model.params))}


	params = {model.params._index[i]: (model.params[i],model.pvalues[i]) for i in range(len(model.pvalues))}
#	print model.pvalues._index[1]
#	print model.params._index[1]

#	print pvs
#	print params

#	RS =  model.rsquared
#	FP =  model.f_pvalue
	RESIDS = [r for r in model.resid]

#		iKeys = [k for k in pvs.keys() if k.split('=')[0] == INTEREST]
#		return [[(k,pvs[k],round(params[k],3)) for k in iKeys],(RS,FP,pvs,params),RESIDS,model]
	return [(model.rsquared,params),RESIDS,model]
		


def logitR(y,X):


	logit = LogisticRegression(penalty='l1')
	logit.fit(X,y)

	y_pred  = logit.predict(X)
	y_probs = logit.predict_proba(X)
	t0,t1,f0,f1 = 0,0,0,0


	
	hi, low  = [] , [] 

	for yp,pr,yr in zip(y_pred,y_probs,y):
		if yp == yr and yr == '1': t1 +=1.0
		elif yp == yr: 		 t0 +=1.0
		elif yr == '1': 		 f1 +=1.0
		else:			 f0 +=1.0
		

		if yp != '1': 
			if yp == yr: low.append([1.0-pr[0],'Y'])
			else: 	     low.append([1.0-pr[0],'N'])

		else:
			if yp == yr: hi.append([pr[0],'Y'])
			else: 	     hi.append([pr[0],'N'])


	scrs = (t1+t0)/float(len(y)),t1/(t1+f1+0.1),t0/(t0+f0+0.1)
	low.sort()
	hi.sort() 

	try: 
		hp = [x for x in hi if x[1] != 'Y'][0][0]
		xp = [x for x in hi if x[0] < hp]
		hist,hps = hp - xp[0][0], len(xp) 
	except IndexError: 
		hist,hps = 0.1,5 	

	try: 
		sp = [x for x in low if x[1] != 'Y'][0][0]
		lp = [x for x in low if x[0] < sp]
		dist,lps = sp - lp[0][0], len(lp) 
	except IndexError: 
		dist,lps = 0.1,5 



	return [(y_pred[i],y_probs[i]) for i in range(len(y_pred))],scrs,(lps,dist),(hps,hist)









#!/usr/bin/env python




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











class Regression:
	def __init__(self,rage):

		self.rage = rage 
		self.progress = rage.progress
		self.options = rage.args 


		self.D = self.rage.data 



		if self.options.model == 'OLS':				self.regress = self.regress_ols	
		elif self.options.model.upper() in ['GLM-NB']:		self.regress = self.regress_glmnb
		elif self.options.model.upper() in ['NBINOM','NB']:	self.regress = self.regress_nb 
		elif self.options.model.upper() in ['ZIP']: 		self.regress = self.regress_zip

















	def run(self): 

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

