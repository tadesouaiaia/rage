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

from scipy.signal import savgol_filter as svgf 
from math import exp
from math import factorial 


#from ..Rage_Transforms import rage_KDE

#from ..Rage_Transforms import rage_DR
from ..Rage_IO import rage_outputs 



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


def perc_diff(a,b): 
	delta = fabs(a-b) 
	return round(100 * (delta / ((a+b) / 2.0)),2) 



class features: 
        def __init__(self,rage):
		self.progress = rage.progress
		self.args = rage.args 
		self.D    = rage.data 

	def get_f_ratios(self,RSTAT='MED'): 


		self.r_key = {}
		self.mean_key = {}
		self.RSTAT = RSTAT 


		self.HOUSEKEEPING = dd(bool) 
		if self.args.ratios != None:
   			self.progress.start_minor('Reading DropOut Aware Feature Parameters From File',len(self.D.features)*(len(self.D.features)/2.0),True)
			f_idx = {f.name: f.idx for f in self.D.features}
			r_header = self.args.ratios.readline() 
			for line in self.args.ratios:
				self.progress.mark()
				line = line.split() 
#				f1,f2,fMatch,fRate = line[0],line[1],int(line[-2]),float(line[-1])
				f1,f2,fRate = line[0],line[1],float(line[-1])

				try: 
					self.HOUSEKEEPING[f_idx[f1]] = True
					self.HOUSEKEEPING[f_idx[f2]] = True
					self.r_key[(f_idx[f1],f_idx[f2])] = fRate
				except KeyError: continue
		else:
   			self.progress.start_minor('Selecting DropOut Aware Feature Parameters',len(self.D.features),True)
			wMed = open(self.args.prefix+'_medRatios.out','w') 
			wMean = open(self.args.prefix+'_meanRatios.out','w') 
			wMean.write('%-40s %40s %5s %5s | %10s %25s %25s\n' % ('---','---','obs1','obs2','matches','medRatio','meanRatio'))
			wMed.write('%-40s %40s %5s %5s | %10s %25s %25s\n' % ('---','---','obs1','obs2','matches','meanRatio','medRatio'))
			for fi in range(len(self.D.features)): 
				self.progress.mark() 
				f_key = dd(list) 
				for si,ci in self.D.features[fi].cnts.items():
					for fj,cj in [st for st in self.D.samples[si].cnts.items() if st[0] > fi]: 
						f_key[(fi,fj)].append(ci/float(cj))


				for fp,fl in f_key.items():
					fLen,fMed,fMean = len(fl), np.median(fl), np.mean(fl) 
#					self.r_key[fp] = fMed
					if RSTAT == 'MED':	self.r_key[fp] = fMed
					else:			self.r_key[fp] = fMean
#					self.mean_key[fp] = fMean
					c1 = len(self.D.features[fp[0]].cnts)
					c2 = len(self.D.features[fp[1]].cnts)
					self.HOUSEKEEPING[fp[0]] = True 
					self.HOUSEKEEPING[fp[1]] = True
					wMed.write('%-40s %40s %5d %5d | %10d %25f %25f\n' % (self.D.features[fp[0]].name,self.D.features[fp[1]].name,c1,c2,fLen,fMean,fMed))
					wMean.write('%-40s %40s %5d %5d | %10d %25f %25f\n' % (self.D.features[fp[0]].name,self.D.features[fp[1]].name,c1,c2,fLen,fMed,fMean))
		
		self.progress.end() 		

		return self






	def predict_missing_sample_values(self,STAT='MEAN'): 

		EXP_TOTALS = []
		self.expected_totals = [] 
		s_totals = [sum(s.cnts.values()) for s in self.D.samples] 		
		w = open(self.args.prefix+'_mR_effective.out','w')
                w.write('%-40s %10s %20s %20s\n' % ('---','idx','actual_total','inferred_total'))
		for si,s in enumerate(self.D.samples): 	
			
			s_missing      = [i for i in range(len(self.D.features)) if i not in s.cnts] 
			s_housekeeping = [i for i in s.cnts if i in self.HOUSEKEEPING] 
			s_pred = {i: 0 for i in range(len(self.D.features))} 
                        self.progress.start_minor('Inputing Missing Data For Sample '+s.name,1+len(s_missing),False)

			for m in s_missing: 
				self.progress.mark() 
				m_infer = [] 
				for i in s_housekeeping: 
					try: 
						if i < m: m_infer.append(s.cnts[i] /  self.r_key[(i,m)]) 
						else:     m_infer.append(s.cnts[i] *  self.r_key[(m,i)]) 
					except KeyError: continue 
				
			
				if len(m_infer) > 0:
					if STAT == 'MED': s_pred[m] = np.median(m_infer) 					
					else:		  s_pred[m] = np.mean(m_infer) 
					#s_pred[m] = np.median(m_infer) 					
					#s_pred[m] = np.mean(m_infer) 

			inferred_missing_total = sum(s_pred.values())
			expected_total = s_totals[si] + inferred_missing_total
			self.expected_totals.append([expected_total, s_totals[si], s.idx])
			w.write('%-40s %10d %20f %20f\n' % (self.D.samples[s.idx].name,s.idx,s_totals[si],expected_total))

			EXP_TOTALS.append(expected_total) 
		EXP_TOTALS = [e/(j+1) for j,e in enumerate(EXP_TOTALS)]

		#print "TOTAL_MEAN-SPAN-VAR-CV", self.RSTAT,STAT,np.mean(EXP_TOTALS), max(EXP_TOTALS) - min(EXP_TOTALS) , np.var(EXP_TOTALS), coVar(EXP_TOTALS) 

		return self









	def predict_known_ratio_values(self): 
		f_pred, s_pred = dd(lambda: dd(float)), dd(lambda: dd(float))
		f_totals = [sum(f.cnts.values()) for f in self.D.features] 
		f_log_totals = [log(ft) for ft in f_totals] 
		s_totals = [sum(s.cnts.values()) for s in self.D.samples] 
		s_log_totals = [log(st) for st in s_totals] 
		for si,s in enumerate(self.D.samples): 	
			s_contained    = [i for i in s.cnts] 
			s_housekeeping = [i for i in s.cnts if i in self.HOUSEKEEPING] 
			for m in s_contained: 
				m_infer = [] 
				for i in s_housekeeping: 
					if i == m: continue
					try:
						if i < m: m_infer.append(s.cnts[i] / self.r_key[(i,m)]) 
						if i > m: m_infer.append(s.cnts[i] * self.r_key[(m,i)]) 
					except KeyError: 
						continue 
				if len(m_infer) == 0: infer_val = 0 
				else:                 infer_val = np.mean(m_infer) 
				#else:                 infer_val = np.mean(m_infer) 
				f_pred[m][s.idx] = infer_val 
				s_pred[s.idx][m] = infer_val 
			
		wf = open(self.args.prefix+'_summarize_featureRatios.out','w') 
                wf.write('%-40s %15s %15s %6s %18s %10s %10s %10s %10s\n' % ('---','total_reads','total_obs','cv','predicted_total','perc_diff','R-depth','R-pred','R-log-pred'))
		for fi,f in enumerate(self.D.features): 

			f_key, f_name = f.cnts.keys(), f.name 
			f_true, f_log_true =        [f.cnts[k] for k in f_key] , [log(f.cnts[k]) for k in f_key] 
			f_predicted = [f_pred[fi][k] for k in f_key]
			f_predicted_total = sum(f_predicted) 
			f_cv = coVar(f_true) 
			p_diff = perc_diff(f_predicted_total, f_totals[fi]) 
			fs_log_totals =  [s_log_totals[k] for k in f_key]
			fs_totals =      [s_totals[k] for k in f_key]
			fRT =  stats.pearsonr(fs_log_totals,f_log_true)[0] 
			fRP = stats.pearsonr(f_predicted, f_true)[0]
			fRLP = stats.pearsonr([log(x) for x in f_predicted], f_log_true)[0] 
                	wf.write('%-40s %15d %15d %6.2f %18.1f %10.2f %10.2f %10.2f %10.2f\n' % (f.name,f_totals[fi],len(f.cnts),coVar(f_true),f_predicted_total,p_diff,fRT,fRP,fRLP))
		
		ws = open(self.args.prefix+'_summarize_sampleRatios.out','w') 
                ws.write('%-40s %15s %15s %6s %18s %10s %10s %10s %10s\n' % ('---','total_reads','total_obs','cv','predicted_total','perc_diff','R-depth','R-pred','R-log-pred'))
		for si,s in enumerate(self.D.samples): 

			s_key, s_name = s.cnts.keys(), s.name 
			s_true, s_log_true =        [s.cnts[k] for k in s_key] , [log(s.cnts[k]) for k in s_key] 
			s_predicted = [s_pred[si][k] for k in s_key]
			s_predicted_total = sum(s_predicted) 
			s_cv = coVar(s_true) 
			p_diff = perc_diff(s_predicted_total, s_totals[si]) 
			fs_log_totals =  [f_log_totals[k] for k in s_key]
			fs_totals =      [f_totals[k] for k in s_key]
			fRT, fRP, fRLP =  stats.pearsonr(fs_log_totals,s_log_true)[0] , stats.pearsonr(s_predicted, s_true)[0], stats.pearsonr([log(x) for x in s_predicted], s_log_true)[0] 
                	ws.write('%-40s %15d %15d %6.2f %18.1f %10.2f %10.2f %10.2f %10.2f\n' % (s.name,s_totals[si],len(s.cnts),coVar(s_true),s_predicted_total,p_diff,fRT,fRP,fRLP))
		self.progress.end() 
		return self



















