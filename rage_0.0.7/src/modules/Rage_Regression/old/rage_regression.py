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












 




