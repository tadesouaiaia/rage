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

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection as fdr

# import formula api as alias smf
import statsmodels.formula.api as smf
# formula: response ~ predictor + predictor

from scipy.stats import pearsonr as pearsonr


import random
from math import fabs
#from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
import pickle
from math import log
import pylab 
from matplotlib.patches import Rectangle as Rect
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib._png import read_png
from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import Circle as Circ
from operator import mul
from sklearn.linear_model import LogisticRegression




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







class model:
	def __init__(self,progress=None):
		self.progress  = progress 
		self.variables = {} 
		self.types = {} 
		self.minGroupSize = 10 

	def load_dataset(self,dataset):

		self.samples = dataset.samples  
		self.features = dataset.features
		self.feature_vals     = dataset.feature_cnts
		self.feature_fracs   = dataset.feature_fracs 
		return self

	def add_covariate_key(self,key):

		for k in key:
			k_vals = [key[k][s] if s in key[k] else 'NA' for s in self.samples]
			try: 
				self.variables[k] = [float(x) if x != 'NA' else 'NA' for x in k_vals]
				self.types[k] = 'continuous'
			except ValueError:
				self.variables[k] = k_vals
				self.types[k] = 'binary' 

	

	def prepare_variables(self,variables):
		#sample_idxs = self.prepare_idxs(variables)
		values, labels, idxs = [],[], [i for i in range(len(self.samples)) if 'NA' not in [self.variables[v][i] for v in variables]]
		for v in variables:
			if self.types[v] == 'binary':
				kVals,kCount = [self.variables[v][i] for i in idxs], cc([self.variables[v][i] for i in idxs])
				kPass = [kn for (kn,kv) in kCount.items() if kv >= min(self.minGroupSize,max(kCount.values()))]
				if len(kPass) == len(kCount.keys()):	kPass = kPass[1::]
				kV,kL = [[1 if v == opt else 0 for v in kVals] for opt in kPass],[v+'='+opt for opt in kPass]
			else:
				kV,kL = [[self.variables[v][i] for i in idxs]],[v]
			values.extend(kV)
			labels.extend(kL) 
		return values,labels,idxs

#  		X_data =  pd.DataFrame(np.array([np.array([values[j][i] for j in range(len(values))]) for i in range(len(idxs))]), columns = labels)
#		Y_data =  pd.DataFrame(np.array([np.array([self.feature_vals[j][i] for j in range(len(self.feature_vals))]) for i in idxs]),columns = self.features)
#		Y_data =  pd.DataFrame(np.array([np.array([log(1.0+self.feature_vals[j][i]) for j in range(len(self.feature_vals))]) for i in idxs]),columns = self.features)
#		Y_data =  pd.DataFrame(np.array([np.array([counts[i] for counts in D.counts]) for i in idxs]), columns = D.features)
#		return idxs,[self.samples[i] for i in idxs],X_data,Y_data

	def run_analysis(self,variables,interest=None):

		values,labels,idxs = self.prepare_variables(variables) 

		samples = [self.samples[i] for i in idxs]
  		X_data  =  pd.DataFrame(np.array([np.array([values[j][i] for j in range(len(values))]) for i in range(len(idxs))]), columns = labels)
		Y_data  =  pd.DataFrame(np.array([np.array([log(1.0+self.feature_vals[j][i]) for j in range(len(self.feature_vals))]) for i in idxs]),columns = self.features)
		#X_full, X_interest, X_cov = sm.add_constant(X_data[X_cols]) , sm.add_constant(X_data[X_INTS]), sm.add_constant(X_data[X_COVS]) 
		#X_full, X_interest, X_cov = sm.add_constant(X_data[X_cols]) , sm.add_constant(X_data[X_INTS]), sm.add_constant(X_data[X_COVS]) 
		X_full = sm.add_constant(X_data[X_data.columns]) 
	#	xRun = X_full
		#xBoth =  X_data['TOTAL_FEATURES']*X_data['TOTAL_READS']
#		xBoth =  X_data['TOTAL_READS']/X_data['TOTAL_FEATURES']
		xBoth =  X_data['TOTAL_READS']*X_data['TOTAL_FEATURES']
		xRun = X_data.assign(xBoth=xBoth)
	#	xRun = sm.add_constant(xRun)
		z=0 
		zR= 0 
		for feature in Y_data.columns:
#			pd.concat(X_full,Y_data[feature])
#			xBoth =  X_data['TOTAL_FEATURES']*X_data['TOTAL_READS']
#			df1 = X_data.assign(feature=Y_data[feature])
#			test = smf.ols(formula='feature ~ TOTAL_FEATURES * TOTAL_READS',data = df1).fit()
			model  = sm.OLS(Y_data[feature],xRun).fit()
			RS =  model.rsquared
			FP =  model.f_pvalue
			zR += RS 
			z += 1
			if z %50 ==0: print zR/z
			continue

	def run_observation_analysis(self,variables):

		values,labels,idxs = self.prepare_variables(variables) 

		f_res = {} 
  		X_data  =  pd.DataFrame(np.array([np.array([values[j][i] for j in range(len(values))]) for i in range(len(idxs))]), columns = labels)
#		xBoth =  X_data['TOTAL_READS']*X_data['TOTAL_FEATURES']
		xBoth =  X_data['TOTAL_FEATURES']*X_data['LOG_TOTAL_READS']
	#	xRun = X_data.assign(LOG_TOTAL_READSxTOTAL_FEATURES=xBoth) 
	#	xRun = sm.add_constant(xRun) 
		xRun = sm.add_constant(X_data) 
		Y_obs   =  pd.DataFrame(np.array([np.array(['1' if self.feature_vals[j][i]>0 else '0' for j in range(len(self.feature_vals))]) for i in idxs]),columns = self.features)
		Y_data  =  pd.DataFrame(np.array([np.array([log(1.0+self.feature_vals[j][i]) for j in range(len(self.feature_vals))]) for i in idxs]),columns = self.features)
	#	for feature,log_vals in zip(self.features,[[log(self.feature_vals[j][i]+1.0) for i in idxs] for j in range(len(self.feature_vals))]):



		for feature,log_vals in zip(self.features,[[self.feature_fracs[j][i] for i in idxs] for j in range(len(self.feature_vals))]):


			if feature.split(';')[1] in ['boo']: continue 

			y_all = pd.DataFrame(np.array(np.array(log_vals)), columns = [feature])[feature]
			reg_all = regress(y_all,xRun)
			RS_all,params_all = reg_all[0] 


			logitObs,logitScrs,logLow,logHi = logitR(Y_obs[feature],xRun)


			if logitScrs[0] < 0.4 or logitScrs[1] < 0.08 or logitScrs[2] < 0.33: continue 	
			if len([l for l in logitScrs if l < 0.66]) > 1: continue 
			if logitScrs[0] < 0.66 and logitScrs[2] < 0.66 and logitScrs[1] < 0.66: continue 
	
			if sum(logitScrs) < 1.3: continue 


			if logLow[0] < 10 or logLow[1] < 0.03 or logHi[0] < 10 or logHi[1] < 0.03: continue 
			if logLow[1] + logHi[1] < 0.11: continue 

			obs_tuples = [(k,v) for k,v in enumerate(log_vals) if v >0]
			obs_idxs,obs_vals = [x[0] for x in obs_tuples],[x[1] for x in obs_tuples]
			y_obs = pd.DataFrame(np.array(np.array(obs_vals)), columns = [feature])[feature]
  			X_obs  =  pd.DataFrame(np.array([np.array([values[j][m] for j in range(len(values))]) for m in obs_idxs]), columns = labels)
#			xBoth =  X_obs['TOTAL_READS']/X_data['TOTAL_FEATURES']
#			xBoth =  X_obs['LOG_TOTAL_READS']*X_data['TOTAL_FEATURES']
#			xObs = X_obs.assign(LOG_TOTAL_READSxTOTAL_FEATURES=xBoth) 
			#xObs = sm.add_constant(xObs)

			xObs = sm.add_constant(X_obs)
			reg_obs = regress(y_obs,xObs) 
			RS,params= reg_obs[0] 

			dG,dL = 0, 0 
			s_vals, t_preds, v_preds, t_diffs, v_diffs = [],[],[],[],[] 
			preds = []

			t_px = []  
			for i,idx in enumerate(idxs):
				s_id, s_key, s_val, v_pred, t_pred = self.samples[idx], {v: self.variables[v][idx] for v in variables}, log_vals[i], 0, 0
				for p,(v,sr) in params.items():
					if p == 'const': 	v_pred += v 	
					elif p in s_key:	v_pred += v*s_key[p]  
					else:			v_pred += v*reduce(mul,[s_key[ps] for ps in p.split('x')])
				for p,(v,sr) in params_all.items():
					if p == 'const': 	t_pred += v 	
					elif p in s_key:	t_pred += v*s_key[p]  
					else:			t_pred += v*reduce(mul,[s_key[ps] for ps in p.split('x')])


				if v_pred < -1: v_pred = -1
				if t_pred < -1: t_pred = -1 
				preds.append([v_pred,t_pred,log_vals[i]])	

				if s_val > 0:
					if v_pred < log(2.0): v_pred = log(2.0)
					s_vals.append(s_val)
					t_preds.append(t_pred) 
					v_preds.append(v_pred) 
					v_diffs.append((s_val-v_pred)*(s_val-v_pred))
					t_diffs.append((s_val-t_pred)*(s_val-t_pred))
					#if v_pred > s_val: dG +=1
					#else:		   dL +=1 	

			obR,obP = pearsonr(v_preds,s_vals)
			allR,allP = pearsonr([x[1] for x in preds],[x[2] for x in preds])

			if allR > 0.40 or obR < 0.5: continue 
		#	print feature, allR,allP, obR, obP , '|',logLow[0],logLow[1],logHi[0],logHi[1],feature
		#	if pearsonr(v_preds,s_vals)[0] < 0.5: continue 

 
			f_res[feature] = {'res': preds, 'bin': logitObs, 'obs': reg_obs, 'all': reg_all}
			if len(f_res.keys()) > 4: 
				break
		return f_res 











	def evaluate_model(self,outstr):
		idxs,samples,X_data,Y_data = self.prepare_data(self.members) 
		X_cols, Y_cols = [xc for xc in X_data.columns], [yc for yc in Y_data.columns]
   		out = open(outstr,'w')
		vDict = {xc: 'INT' for xc in X_cols if xc.split('=')[0] in self.interests}
		for feat in Y_cols: 			
			y,res = Y_data[feat],dd(bool) 
			for xCol in X_cols:
				if xCol not in vDict.keys(): vDict[xCol] = 'COVAR'
				X_eval = sm.add_constant(X_data[[xCol]])
				res[xCol] = reg_choice(y, X_eval)[1][0:2] 
	
			for i,aCol in enumerate(X_cols):
				for j in range(i+1,len(X_cols)):
					bCol = X_cols[j] 
					sR,sF = res[aCol]
					iR,iF = res[bCol]
					X_pair = sm.add_constant(X_data[[aCol,bCol]])
					pR,pF,pvalK,paramK = reg_choice(y, X_pair)[1] 
					Pscr,Rscr = (iF/pvalK[bCol])-1, (pR/(sR+iR))-1
					if   Rscr < -0.25 and Pscr < -0.25: verdict = 'COLINEAR'
					elif Rscr >  0.25 and Pscr  > 0.25: verdict = 'INDEP'
					else:			    verdict = 'UNCLEAR'
					aT,bT = vDict[aCol],vDict[bCol]
					out.write('%-30s %30s %6s %8.5f %8.3e %30s %6s %8.5f %8.3e | %10f %10f %10f %8s\n' % (feat,aCol,aT,sR,sF,bCol,bT,iR,iF,pR,Rscr,Pscr,verdict))
	
			



	def evaluate_vars(self,outstr):

		idxs,samples,X_data,Y_data = self.prepare_data(self.members) 
		X_cols, Y_cols = [xc for xc in X_data.columns], [yc for yc in Y_data.columns]
   		out = open(outstr,'w')
		out.write('%-35s %25s %15s %10s %10s  | %6s %6s %5s %5s %5s  \n' % ('---','VAR','TYPE','R_2','FV','MEAN','OBS','mR','oR','TOT'))
#		if self.dex.args.verbose: sys.stderr.write('Update: Printing Partial Rs\n') 
		for feat in Y_cols: 			
			y,res = Y_data[feat],dd(bool) 
			yMean = sum(y)/float(len(y))
			bin_dict = dd(list) 
			re_dict = {}
			for c,xCol in enumerate(X_cols):
				xGroup = xCol.split('=')[0] 
				X_eval = sm.add_constant(X_data[[xCol]])
				sR,sF = reg_choice(y, X_eval)[1][0:2] 
				if self.dex.types[xGroup] == 'binary': 
					xCnts = [b for (a,b) in zip(X_eval[xCol],y) if a != 0] 
					xMean,xObs = sum(xCnts)/float(len(xCnts)), len([a for a in xCnts if a>0])/float(len(xCnts))
					bin_dict[xGroup].append([xMean,xObs,sR,sF,xCol])
					if xGroup not in re_dict.keys(): re_dict[xGroup] = [1 for xc in range(len(X_eval[xCol]))]
					for j,xc in enumerate(X_eval[xCol]): 
						if xc > 0: re_dict[xGroup][j] = 0
				else:
					bin_dict[xCol] = [[sR,sF]] 
			for b,K in bin_dict.items():
				if len(K) == 1:	out.write('%-35s %25s %15s %10.5f %10.2e  | %6s %6s %5s %5s %5s  \n' % (feat,b,'CONT',K[0][0],K[0][1],'NA','NA','NA','NA','1'))
				else:
					missing = re_dict[b]
  					X_tmp =  pd.DataFrame(np.array([np.array([m])  for m in missing]), columns = ['MISSING'])
					X_miss = sm.add_constant(X_tmp[['MISSING']])
					sR,sF = reg_choice(y, X_miss)[1][0:2] 
					xCnts = [bb for (aa,bb) in zip(X_miss['MISSING'],y) if aa != 0] 
					xMean,xObs = sum(xCnts)/float(len(xCnts)), len([a for a in xCnts if a>0])/float(len(xCnts))
					bin_dict[b].append([xMean,xObs,sR,sF,b+'=MISS'])
					mRank,oRank = {},{} 
					mSrt = sorted(K,reverse=True) 
					members = len(mSrt) 
					for i in range(len(mSrt)):
						if i == 0:
							if mSrt[i][0] > mSrt[i+1][0]: mRank[mSrt[i][-1]]=1 
							else:			      mRank[mSrt[i][-1]]=2
						
						else:
							if mSrt[i][0] == 0: mRank[mSrt[i][-1]] = members 
							else: 		    mRank[mSrt[i][-1]] = i+1 
					mSrt = sorted(K,key = lambda x: x[1],reverse=True) 
					for i in range(len(mSrt)):
						if i == 0:
							if mSrt[i][1] > mSrt[i+1][1]: oRank[mSrt[i][-1]]=1 
							else:			      oRank[mSrt[i][-1]]=2
						else:
							if mSrt[i][1] == 0: oRank[mSrt[i][-1]] = members 
							else: 		    oRank[mSrt[i][-1]] = i+1 
					for kMean,kObs,kR,kP,kName in K: 
						out.write('%-35s %25s %15s %10.5f %10.2e  | %6.3f %6.3f %5d %5d %5d  \n' % (feat,b,kName,kR,kP,kMean,kObs,mRank[kName],oRank[kName],members))
					
		



	def check_model(self,interest):

		D,res,rstats = self.dex, {} ,{}
		idxs,samples,X_data,Y_data = self.prepare_data([interest]) 
		ids,opts =  [D.key[interest][s] for s in samples],list(set([D.key[interest][s] for s in samples]))
		X_cols, Y_cols = [xc for xc in X_data.columns], [yc for yc in Y_data.columns]
		X_INT, X_SVA = [xc for xc in X_cols if xc.split('=')[0] == interest], [xc for xc in X_cols if xc.split('=')[0] != interest]
		X_full,X_interest, X_sva = sm.add_constant(X_data[X_cols]) , sm.add_constant(X_data[X_INT]), sm.add_constant(X_data[X_SVA]) 

		INIT=False
		
		
		if D.types[interest] == 'binary':
			for feature in Y_cols: 
				y = Y_data[feature]
				opt_cnts = [[y[i] for i in range(len(y)) if D.key[interest][samples[i]] == opt] for opt in opts]
				opt_lists = [(list_stats(opt_cnts[x]),opts[x]) for x in range(len(opts))]
				rstats[feature] = list_stats(y) 
				res[feature] = opt_lists 			


		self.stats = rstats					

		return res













































































	def evaluate_markers(self,outstr):
		out = open(outstr,'w') 
	#	out = sys.stdout

		
		out.write('%-30s %25s %10s | %13s %13s %13s | %10s\n' % ('---','INTEREST','OBS','obsN','totN','enrich','stat-data'))
		for (feature,counts) in zip(self.dex.features,self.dex.counts):
			for interest in self.results.keys():
				if self.dex.types[interest] != 'binary': continue
				#for k in self.results[interest][feature]:
				#for k in self.results[interest][feature]:
				#for k in range(1):
				scores,cnts = self.results[interest][feature]['ttest']
				if min([s[1] for s in scores]) > 0.01: continue 
				if max([s[2] for s in scores]) < 2.0: continue 
				my_cnts = [] 
				for c in cnts.keys():	my_cnts.extend([(s,c) for s in cnts[c][-1]])
				my_cnts.sort(reverse=True) 
				my_cc,my_len = cc([m[1] for m in my_cnts]),len(my_cnts)
				total = sum(my_cc.values())
				
				tots,exp,obs,n,enrich = {x: y for x,y in my_cc.items()},{x: y/float(total) for  x,y in my_cc.items()},dd(float),0.0,dd(list)  

				for i,(mS,m) in enumerate(my_cnts):
					if mS == 0: break 
					n += 1.0
					obs[m] += 1.0 
					if i > 10 and i%3==0: 
						for a,aX in obs.items():
							ech = (aX/n) / exp[a]
							if ech > 1.25: 	
								if len(enrich[a])  < 3 and ech > 1.25: 				enrich[a].append([ech,aX,n]) 
								elif len(enrich[a])>=3 and ech > enrich[a][-1][0]: 		enrich[a].append([ech,aX,n]) 
					

				if len(enrich.keys()) == 0: continue  			
				for (e,eX) in enrich.items():
					if obs[e] > tots[e]/3.0: 
						my_scrs = [s for s in scores if e in s[0]]
						out.write('%-30s %25s %10.3f |' % (feature,e,obs[e]/float(tots[e])))
						out.write(' %13s %13s %13.3f |' % (eX[-1][1],eX[-1][2],eX[-1][0]))
						#print feature,e,round(obs[e] / float(tots[e]),3)
						#print '|',eX[-1][1],eX[-1][2],round(eX[-1][0],3),
						for (a,b),pv,fc in my_scrs: 
							
							out.write(' %10s %5.2e  %5.2f' % (a+','+b,pv,fc))
						out.write('\n')
						#print '|',eX[-1][1],eX[-1][2],round(eX[-1][0],3),










	def run_shuffled_sims(self):
		self.sim_results = {} 


		xLoc,yLoc = 0,0 
		out = open(self.args.output+'.reg_simulate.txt','w') 
		#out = sys.stdout
		out.write('%-15s %20s %13s %13s %13s %13s\n' % ('---','p-value','obs','sims','sim-obs','empirical-pv'))
		for interest in self.interests:		


			#params = 
			#res =  sorted([vv[2][0] for v in self.regression_results[interest].values()]) 
			pvs =  sorted([v[2][0] for v in self.regression_results[interest].values()]) 
			res =  sorted([(v[2][0],v[1][0]) for v in self.regression_results[interest].values()]) 
			simList,realCnt = dd(list), dd(int) 

			if pvs[0] < 0.000000001: pvBins = [0.000000001,0.0000001,0.00001,0.001,0.01]
			elif pvs[0] < 0.0000001: pvBins = [0.0000001,  0.000001, 0.00001,0.001,0.01]
			elif pvs[0] < 0.000001:  pvBins = [0.000001,   0.00001,  0.0001, 0.001,0.01]
			elif pvs[0] < 0.00001:   pvBins =  [0.00001,    0.0001,   0.001, 0.01,0.05]
			else: 	                 pvBins =  [0.0005,    0.001,   0.005, 0.01,0.05]
				
			for pv in pvs: 
				for v in pvBins: 
					if pv <= v: realCnt[v] +=1 
		
			for n in range(self.args.shuffle): 
				pvKey = dd(int)
				random_pvs = [s[2][0] for s in self.reg_model(interest,'SHUFFLE').values()]
				random_res = [(s[2][0],s[1][0]) for s in self.reg_model(interest,'SHUFFLE').values()]
				for r in random_pvs: 
					for p in pvBins: 
						if r <= p: pvKey[p] +=1

				for v in pvBins: simList[v].append(pvKey[v]) 			

			barData = [[v,realCnt[v],simList[v]] for v in sorted(pvBins)]
			if yLoc == self.yLen: 
				xLoc,yLoc = xLoc+1,0

			ax = plt.subplot2grid((self.xLen, self.yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
			bars = BarPlot(barData)  
			bars.draw(ax,interest) 	
			yLoc +=1


			for v,realCnt,sL in barData:
				s_obs = len([x for x in sL if x>=realCnt])
				s_pv = s_obs / float(len(sL))
				out.write('%-15s %20f %13d %13d %13d %13.5f\n' % (interest,v,realCnt,len(sL),s_obs,s_pv))

 
 		plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.90,wspace=0.3,hspace=0.4) 

		plt.savefig(self.args.output+'_'+self.covStr+'_'+'fig_sim.png')
		if not self.args.hide: plt.show() 












 




#!/usr/bin/env python


import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
import scipy.stats as stats
from scipy.stats import variation as coVar 
import scipy.spatial.distance as sds 
from random import random
import numpy as np
import itertools
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

import statsmodels.api as sm
from statsmodels.stats import power as smp 
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif 
from operator import mul

from scipy.signal import savgol_filter as svgf 
from math import exp
from math import factorial 

from scipy.stats import chisquare
from scipy.stats import ttest_ind 
from Rage_IO import rage_outputs
from Rage_Plots import rage_regression_plots
from Rage_Transforms import rage_KDE
from Rage_Transforms import rage_DR
from Rage_Regression import rage_regression 
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




	def run_model(self):
		#seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'pink'})

		print len(self.D.samples) 
		self.V = self.D.filter_samples_by_attributes(self.options.predictors,self.options.covariates).set_sample_variables(combine=False)
		
		X = self.V.select_variables() 


		
#	
#		if self.D.COMBINED_PREDICTOR:
#			self.seg,self.seg_lens = self.D.samples.segregate(self.D.predictor)	
#		else:
#			print 'hello'
#			print self.D.predictors
		


                self.progress.start_minor('Running Model Regression',len(self.D.features),True)


		model_result = self.regression_result(X,True, mtype = 'full')
		rage_outputs.regression_result(self.rage.args,model_result['dex']).write(self.D.predictor,self.D.covariates)

		print 'uh' 
		self.progress.end() 
		sys.exit() 





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




	


	def evaluate_model(self,simulations=3):
		#seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'pink'})
		self.D = self.rage.data.filter_samples_by_attributes(self.options.predictors,self.options.covariates).set_sample_variables()
		X_names,  X = self.D.create_regression_array() 



		self.seg,self.seg_lens = self.D.samples.segregate(self.D.predictor)	
                self.progress.start_minor('Running Model Regression',len(self.D.features),False)
		self.dim_red = rage_DR.DR(self.options,False,len(self.D.samples)) 
		pca_init = self.dim_red.run_pca(self.D.matrix('log'),req='brief')

		model_result = self.regression_result(X_names,X,True, mtype = 'full')
		rage_outputs.regression_result(self.rage.args,model_result['dex']).write(self.D.predictor,self.D.covariates)

		pca_c_resid, pca_resid =    self.dim_red.run_pca(np.matrix(model_result['c-resids']).getT(),req='brief'), self.dim_red.run_pca(np.matrix(model_result['resids']).getT(),req='brief')
		my_pvs =  sorted([model_result['params'][i]['predictors'][0][0] for  i in range(len(model_result['params']))])
		my_rs =  model_result['rs']

		pv_counts, rs_counts, pv_key, rs_key = set_result_counter(model_result) 
		sim_vars, sim_rs, sim_pvs = self.run_shuffle_simulations(simulations,sim_type='compare',compare_key = {'pv': pv_key,'rs': rs_key})

		mplot = rage_regression_plots.mplot(3,2,self.options,{'r_key': rs_key, 'p_key': pv_key}) 

		if self.options.featureKey == None:  mplot.add_model_table(model_result,{'interest': self.D.predictor, 'preds': self.D.predictors,'pvs': (pv_counts,sim_pvs)}).update()
		else:	mplot.add_model_table(model_result,{'interest':self.D.predictor,'preds':self.D.predictors,'ss':self.calculate_ss(self.D.features,model_result['params'],self.D.predictor)}).update()

 
                self.progress.start_minor('Plotting Results  ',100,False)
		skrees = [pca_init['var_rates'],pca_c_resid['var_rates'],pca_resid['var_rates'],sim_vars] 
		mplot.add_predictor_table(model_result,{'interest': self.D.predictor, 'preds': self.D.predictors,'skrees': skrees}).update()
		
		mplot.add_rs_bars(rs_counts,sim_rs).update({'title': '$'+"\_".join(self.D.predictor.split('_'))+'$ '+'$\  R^2\ Values$'})
		mplot.add_pv_bars(pv_counts,sim_pvs).update({'title': '$'+"\_".join(self.D.predictor.split('_'))+'$ '+'$\  P\ \ Values$'})
		mplot.add_pca_pts(pca_init['pts'],self.D.samples,self.D.predictor,{'colspan':2}).update({'title': 'PCA Initial Values','yadd': 2,'colspan':2})
		mplot.add_pca_pts(pca_c_resid['pts'],self.D.samples,self.D.predictor,{'colspan':2}).update({'title': 'PCA Covariate Residuals','yadd': 2,'colspan':2})
		mplot.add_pca_pts(pca_resid['pts'],self.D.samples,self.D.predictor,{'colspan':2}).update({'title': 'PCA Model Residuals','yadd': 2,'colspan':2})
		mplot.save_mfig(self.options.model,[self.D.predictor],self.D.covariates)

		self.progress.end() 














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





	def score_vif(self,X,X_names):


		print X_names		

		print self.D.sample_predictors, self.D.sample_covariates

		vd,vd_out = dd(list),{}

		sys.exit() 
		X_tmp = [[x[i] for i in range(len(x)) if X_names[i] != 'intercept'] for x in X]
		N_tmp = [n for n in X_names if n != 'intercept']
		N_iter = [n for n in range(len(N_tmp)) if self.D.variable_class[N_tmp[n]] != 'covariate']

		XX = np.array(X_tmp) 
		for j,n in enumerate(N_tmp):
			vd_out[n] = round(vif(XX,j),2) 
		
		return vd_out 
	

		if len(X_tmp[0]) > 2 and len(N_iter) > 1:

			for i in range(len(N_iter)):
				X_iter = [[x[j] for j in range(len(x)) if j != i] for x in X_tmp]
				X_iter_names = [N_tmp[j] for j in range(len(N_tmp)) if j != i]
				XX = np.array(X_iter)
				for j,n in enumerate(X_iter_names):
					vd[n].append(vif(XX,j))

			for n,v in vd.items():
				vd_out[n] = round(np.mean(v),2)
		else:
			XX = np.array(X_tmp) 
			for j,n in enumerate(N_tmp):
				vd_out[n] = round(vif(XX,j),2) 
		
		return vd_out 

	
	def regression_result(self,X_names,X,log_transform=True,mtype = 'pvs'):


		res = dd(list) 
		if mtype.upper() == 'FULL':
			res['dex'], res['vif'] = {}, self.score_vif(X,X_names) 
			Xi = [[x[i] for i in range(len(x)) if self.D.variable_class[X_names[i]] != 'covariate'] for x in X]
			Xi_names = [n for n in X_names if self.D.variable_class[n] != 'covariate']
			Xc = [[x[i] for i in range(len(x)) if (X_names[i] == 'intercept') or (self.D.variable_class[X_names[i]] == 'covariate')] for x in X]
			Xc_names = [n for n in X_names if (n=='intercept') or (self.D.variable_class[n] == 'covariate')]
				
		for f in self.D.features:

			if sum(f.cnts) == 0: continue 
			self.progress.mark() 
			res['features'].append(f.name) 
			y = [s.cnts[f.idx] for s in self.D.samples]

			if mtype.upper() == 'FULL':	reg_result = self.regress(y,X,X_names,{'req': ['pwr','resids']})
			elif mtype.upper() == 'BRIEF':	reg_result = self.regress(y,X,X_names,{'req': []})
			else:				reg_result = self.regress(y,X,X_names,{'req': ['resids']})
			
			for k,v in reg_result.items():	res[k].append(v)

			if mtype.upper() == 'FULL':
				res['dex'][f.name] = self.dex_score(reg_result,y,Xi,Xi_names)
				res['c-resids'].append(self.regress(y,Xc,Xc_names,{'req': ['resids']})['resids'])

		if mtype.upper() != 'BRIEF':
			return res

		else:
			fLen, b_out = len(res['features']),{}
			f_1,f_10  = int(fLen*0.01+0.999),int(fLen*0.1+0.999)
			res_rs  =       res['rs']

			if len(res['params'][0]['predictors']) > 0:
				res_pvs = 	[res['params'][j]['predictors'][0][0] for j in range(len(res['params']))]
			else:
				res_pvs =   [res['params'][j]['covariates'].values()[0][0] for j in range(len(res['params']))]
		
			for a,b in zip([res_rs,res_pvs],['RS','PV']):

				if b == 'RS': t_srt = sorted(a,reverse=True) 
				else:	      t_srt = sorted(a) 
				t_tops = [np.mean(t_srt),t_srt[f_10],t_srt[f_1]]
				b_out[b] = [a,t_tops]			
			return b_out 







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














































