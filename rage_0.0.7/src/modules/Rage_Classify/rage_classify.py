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
import rage_classify_outputs
#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 
import warnings
warnings.filterwarnings("ignore")


def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]










class BinaryClassifier:

        def __init__(self,options, data):
		self.options = options

		self.names = [f.name for f in data.features] 
		self.sample_names = [s.name for s in data.samples] 
		self.marker_stats = dd(lambda: dd(bool)) 


	def findMarkers(self,target,target_names,Y,F,OPT_MIN=0.10,MFC=1.2,MFF=1.3):

	
		for i,y in enumerate(Y):
			t_total = float(len([t for t in target if type(t) == int]))
			t_vals = [y[j] for j in range(len(y)) if type(target[j]) == int]				
			
			t_dict = {target_names[k]:  [y[j] for j in range(len(y)) if target[j] == k] for k in range(len(target_names))}
			t_stats = [[tk,float(len(tv))/t_total, len(tv),np.mean(tv),np.percentile(tv,90),len([v for v in tv if v>0])/float(len(tv))] for tk,tv in t_dict.items()]



			t_first = cc([sorted(t_stats,key=lambda X: X[j],reverse=True)[0][0] for j in range(3,6)])
			t_win = [t for t in t_first if t_first[t] >= 2] 
			if len(t_win) == 0: continue 

			t_opt,t_name = t_win[0], F[i] 
			avgM,avgQ,avgR=np.mean([t[3] for t in t_stats if t[0] != t_opt]),np.mean([t[4] for t in t_stats if t[0] != t_opt]),np.mean([t[5] for t in t_stats if t[0] != t_opt])
			opt_exp, opt_len, opt_mean, opt_qt, opt_rate = sorted(t_stats, key = lambda X: X[4], reverse = True)[0][1::]

			if opt_rate < OPT_MIN: continue 			
			res = {'obs': opt_rate,'opt': t_opt,'ofc':min(5.0,opt_rate/float(avgR+0.001)),'fc':min(5.0,opt_mean / float(avgM+0.001)),'qfc': min(5.0,opt_qt / float(avgQ+0.001))}
			t_cnts,n,m = sorted([(y[j],target_names[target[j]]) for j in range(len(y)) if type(target[j]) == int],reverse=True),3,max(3,len([t for t in t_vals if t == 0]))
			cnt_a, cnt_b, b_nz  = [t for t in t_cnts if t[1] == t_opt], [t for t in t_cnts if t[1] != t_opt], len([t for t in t_cnts if t[0] != 0 and t[1] != t_opt])
			gt,lt,enhi, a,b,A,B, minB = [], [],[],0,0, len(cnt_a), len(cnt_b) , min([bc[0] for bc in cnt_b])
			while a < len(cnt_a): 

				if cnt_a[a][0] <= minB:
					break 
				while b < len(cnt_b) and cnt_a[a][0] <= cnt_b[b][0] :	b += 1
				gt.append((B-b)/float(B))
				lt.append(b / float(B)) 
				a+=1

			res['comps'] = np.mean(gt) / np.mean(lt)
			while n + m < len(t_cnts): 
				t_n, t_m = t_cnts[0:n], t_cnts[0-m::]
				opt_hi, opt_lo = float(len([tn for tn in t_n if tn[1] == t_opt])), float(len([tm for tm in t_m if tm[1] == t_opt]))
				eH = ( opt_hi / (n+0.00001) ) / (opt_exp+0.0000001) 
				enhi.append(eH) 
				n+=1 
				m+=1

			res['enrich'],res['MARKER'] = round(np.mean(enhi),2), False

			if res['ofc'] > MFC and res['obs'] > 0.25 and res['fc'] > MFC and res['qfc'] > MFC and res['obs'] > 0.18 and res['comps']  > MFF and res['enrich'] > MFF: res['MARKER'] = True 	
			self.marker_stats[t_name] = res

		out = rage_classify_outputs.Marker_Output(self.options) 

		for f in F: 
			out.add_marker(f,self.marker_stats[f]) 

#		for y,f in zip(Y,F): 






	def logitR(self,target,target_names,Y,F,left_out=5):

		logit = LogisticRegression(penalty='l1')	
		id_Ys, coef_key = [[y[i] for y in Y] for i in range(len(target))], dd(lambda: dd(list)) 
		out = rage_classify_outputs.Classifier_Output(self.options) 

		TRAIN_IDXS, CAND_IDXS, TEST_IDXS, s_key, p_key = [], [], [] , {}, dd(list) 
		sample_grades, gene_grades = {}, dd(lambda: dd(list)) 
		

		for i in range(len(target)): 
			sample_grades[self.sample_names[i]] =  dd(list) 
			if type(target[i]) == int: 
				s_key[i] = [self.sample_names[i],True,target_names[target[i]]]
				TRAIN_IDXS.append(i) 
				CAND_IDXS.append(i) 
			else:
				s_key[i] = [self.sample_names[i],False,target[i]]
				TEST_IDXS.append(i) 



		while len(CAND_IDXS) > 0: 
			if len(CAND_IDXS) > left_out: test_set   = list([x for x in np.random.choice(CAND_IDXS,left_out,replace=False)])
			else: test_set   			 = CAND_IDXS
			train_set  = [i for i in TRAIN_IDXS if i not in test_set]
			train_opts = [target_names[target[t]] for t in train_set] 

		 	if len(list(set(train_opts))) != len(target_names):
				print 'uh'
				left_out -= 1 
				continue 

			Yj, Tj = [id_Ys[i] for i in train_set], [target[i] for i in train_set]

			
			logit.fit(Yj,Tj)
			p_coeffs = logit.coef_ 
			preds    = logit.predict(id_Ys) 
    			probs =    logit.predict_proba(id_Ys)
			

			for i in range(len(p_coeffs)): 
				for k in range(len(p_coeffs[i])):	
					coef_key[k][i].append(p_coeffs[i][k])




			for j in TEST_IDXS + test_set: 
				pr = sorted(probs[j],reverse=True) 
				p_key[j].append((target_names[preds[j]] , pr[0] - pr[1])) 
				sName,sBool,sTrue = s_key[j] 
				sPredict = target_names[preds[j]]
				if len(p_key[j]) > 0:


					for k,F_name in enumerate(F): 
						jVal = Y[k][j]  

						jMults = sorted({Xk: jVal*Xv[-1] for Xk,Xv in coef_key[k].items()}.items(),key=lambda XX: XX[1]) 
						for cIdx,cVals in enumerate(p_coeffs): 
							if cVals[k] != 0: 
								cMult = cVals[k] * jVal 

						
						if jMults[-1][1] > 0: 
							mTarget = target_names[jMults[-1][0]]
							sample_grades[sName][mTarget].append([F_name,-1*jMults[-1][1]])	
						if jMults[0][1] < 0: 
							mTarget = target_names[jMults[0][0]]
							sample_grades[sName][mTarget].append([F_name,-1*jMults[-1][1]])	
			


					sOutcome = (sTrue,sPredict) 
					if sTrue == sPredict: sBoolOutcome = 'YES'
					else: 		      sBoolOutcome = 'NO'
					
					for t in target_names: 

						tPos = sorted([(f_val,f_name) for f_name,f_val in sample_grades[sName][t] if f_val > 0],reverse=True)
						tNeg = sorted([(f_val,f_name) for f_name,f_val in sample_grades[sName][t] if f_val < 0],reverse=False)


						if sTrue != t: 
							for fi,f_data in enumerate(tPos):  
								gene_grades[f_data[-1]]['FALSE_POS'].append((fi+1,sBoolOutcome,sOutcome))	
							for fi,f_data in enumerate(tNeg): 
								gene_grades[f_data[-1]]['TRUE_NEG'].append((fi+1,sBoolOutcome,sOutcome))	
 						if sTrue == t:
							for fi,f_data in enumerate(tPos):  
								gene_grades[f_data[-1]]['TRUE_POS'].append((fi+1,sBoolOutcome,sOutcome))	
							for fi,f_data in enumerate(tNeg): 
								gene_grades[f_data[-1]]['FALSE_NEG'].append((fi+1,sBoolOutcome,sOutcome))	
 



					out.add_score(s_key[j],p_key[j]) 
					if j in CAND_IDXS: CAND_IDXS.remove(j) 
					elif j in TEST_IDXS: TEST_IDXS.remove(j) 
			  
		



		for g in gene_grades:
			GK = {} 
			for k in gene_grades[g].keys():
				k_data = gene_grades[g][k] 
				rank_all = [kd[0] for kd in k_data] 			
				rank_yes = [kd[0] for kd in k_data if kd[1] == 'YES'] 
				rank_no = [kd[0] for kd in k_data if kd[1] == 'NO'] 

				ram,raL = np.mean(rank_all), len(rank_all) 
				rym,ryL = np.mean(rank_yes), len(rank_yes) 
				rnm,rnL = np.mean(rank_no), len(rank_no) 

				GK[k] = [ram,raL,rym,ryL,rnm,rnL] 
			
			out.add_gene_grades(g,GK) 
	
			
		for k,f_name in enumerate(F): 		
			out.add_coefs(f_name,coef_key[k],target_names) 






























	def logitU(self,target,target_names,Y,F,left_out=5):

		logit = LogisticRegression(penalty='l1')	

		iterations = 5
		known_key, unk_key = {},{} 
		known_idxs,known_vals, unk_idxs,unk_vals,valid_names,valid_target,valid_key = [],[],[],[],[],[],{}
		for i  in range(len(target)): 
			if target_names[target[i]].split('~')[-1].upper()[0:3] == 'UNK': 
				unk_idxs.append(i) 
				unk_vals.append([y[i] for y in Y])
				unk_key[len(unk_idxs)-1] = [self.sample_names[i],target_names[target[i]]]
			else: 							         
				known_idxs.append(i) 
				known_vals.append([y[i] for y in Y])
				known_key[len(known_idxs)-1] = [self.sample_names[i],target_names[target[i]]]
				if target[i] not in valid_key:
					valid_names.append(target_names[target[i]])
					valid_key[target[i]] = len(valid_names)-1
				valid_target.append(valid_key[target[i]])

		


		novel_key = dd(list) 
		iter_key = dd(list) 

		left_idxs = np.random.choice(range(len(valid_target)),left_out,replace=False)
		while True:
			left_vals,left_target = [known_vals[i] for i in left_idxs],[valid_target[i] for i in left_idxs]

			iter_vals,iter_target = [v for i,v in enumerate(known_vals) if i not in left_idxs],[v for i,v in enumerate(valid_target) if i not in left_idxs] 

			if len(list(set(iter_target))) < len(valid_names):
				left_out -=1
				continue 


			logit.fit(iter_vals,iter_target)

			p_coeffs = logit.coef_ 
			pred_unk    = logit.predict(unk_vals) 
    			prob_unk    = logit.predict_proba(unk_vals)
			pred_val    = logit.predict(left_vals) 
			prob_val    = logit.predict_proba(left_vals) 
			 

			for i in range(len(unk_idxs)): 
				novel_key[i].append((valid_names[pred_unk[i]],prob_unk[i][pred_unk[i]]))

			for i,j in enumerate(left_idxs):
				iter_key[j].append((valid_names[pred_val[i]],prob_val[i][pred_val[i]]))


			left_cands = [i for i in range(len(valid_target)) if len(iter_key[i])<4]

			if len(left_cands) == 0: break 
			elif len(left_cands) <= left_out: left_idxs = left_cands
			else:	left_idxs = np.random.choice(left_cands,left_out,replace=False)


		out = rage_classify_outputs.Classifier_Unknown_Output(self.options) 
		

		for i,dubs in novel_key.items():
			votes = dd(float) 
			name,orig_id = unk_key[i] 
			for a,b in dubs: votes[a]+=b 
			scrs = sorted(votes.items(),key=lambda X: X[1],reverse=True)
			out.add_pred(name,orig_id,scrs[0][0],scrs[0][1]/sum([sc[1] for sc in scrs]))

		for i,dubs in iter_key.items():
			votes = dd(float) 
			name,orig_id = known_key[i] 
			for a,b in dubs: votes[a]+=b 
			scrs = sorted(votes.items(),key=lambda X: X[1],reverse=True) 
			out.add_pred(name,orig_id,scrs[0][0],scrs[0][1]/sum([sc[1] for sc in scrs]))


		
		













































































































































































































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
        def __init__(self,options):

		if options.model == 'OLS':   self.regress = RegOLS
	
		self.stats = dd(lambda: {}) 	
                self.rs_key = [0.01,  0.03,   0.05,    0.10,     0.25]
                self.pv_key = [0.01, 0.001, 0.0001, 0.00001,0.0000001]
		self.resid = [] 

	def run(self, Y, X,req=[]):

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
			if 'resids' in req:
				self.resid.append(r.model.resid) 



                self.pv_cnt = [len([p for p in self.pv_mins if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
                self.rs_cnt = [len([p for p in rs if p > self.rs_key[j]]) for j in range(len(self.rs_key))]

                for k,K in zip(['rs','pv'],[rs,sorted(self.pv_mins)]): self.stats[k] = ListStats(K) 
		return self

	




	def summarize(self):
		bic =sorted([r.model.bic for r in self.reg],reverse=True)
		rs, ars = sorted([r.model.rsquared for r in self.reg],reverse=True), sorted([r.model.rsquared_adj for r in self.reg],reverse=True)

		pwr_05,pwr_001 = sorted([r.pwr[0.05] for r in self.reg],reverse=True), sorted([r.pwr[0.001] for r in self.reg],reverse=True)



		try:		   
			self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.predictor_idx]) for r in self.reg]
		except ValueError: 
			try: self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.covariate_idx]) for r in self.reg]
			except: self.pv_mins = [1 for r in self.reg]

                self.pv_cnt = [len([p for p in self.pv_mins if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
                self.rs_cnt = [len([p for p in rs if p > self.rs_key[j]]) for j in range(len(self.rs_key))]




                for k,K in zip(['bic','rs','ars','pv','pwr_05','pwr_001'],[bic,rs,ars,sorted(self.pv_mins),pwr_05,pwr_001]): self.stats[k] = ListStats(K) 
		for i,n in enumerate(self.X.names): 	self.stats[n] =  ListStats(sorted([r.model.pvalues[i] for r in self.reg]))
		return self

	def out(self):
		rs, ars = sorted([r.model.rsquared for r in self.reg],reverse=True), sorted([r.model.rsquared_adj for r in self.reg],reverse=True)
		try:		   self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.predictor_idx]) for r in self.reg]
		except ValueError: self.pv_mins = [min([r.model.pvalues[idx] for idx in self.X.covariate_idx]) for r in self.reg]
                self.pv_cnt = [len([p for p in self.pv_mins if p < self.pv_key[j]]) for j in range(len(self.pv_key))]
                self.rs_cnt = [len([p for p in rs if p > self.rs_key[j]]) for j in range(len(self.rs_key))]
                for k,K in zip(['rs','pv'],[rs,sorted(self.pv_mins)]): self.stats[k] = ListStats(K) 
		return self


	def resids(self):
		if len(self.resid) > 0: return self.resid
		return [r.model.resid for r in self.reg]


	def pv_discovery(self,pvs,disc_pvs = [0.05,0.01,0.001,0.0001]):

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
        def __init__(self,y,X):


		print X.array 

		sys.exit() 


		self.model = sm.OLS(y,X.array).fit() 
		self.y,self.X = y, X 
		self.pwr = {0.05: 0.50, 0.001: 0.1}



	def test_pwr(self,alphas=[0.05,0.001]):
		if self.model.rsquared > 0: 
			df_de  =  len(self.X.names) -1 
			df_num =  len(self.y) - len(self.X.names)  
			f_2 = np.sqrt(self.model.rsquared / (1-self.model.rsquared) )
	
			for a in alphas:	self.pwr[a] = smp.FTestPower().solve_power(effect_size=f_2,df_num=df_num,df_denom=df_de,alpha=a)

		return self	




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

