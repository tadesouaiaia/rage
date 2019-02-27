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
import warnings
warnings.filterwarnings("ignore")

#from ..Rage_Transforms import rage_KDE

#from ..Rage_Transforms import rage_DR
from ..Rage_IO import rage_outputs 
#from ..Rage_Comps import rage_comps
#from ..Rage_Comps import rage_compares as rage_comps
from ..Rage_Comps import rage_feature_comparisons as rage_comps

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








class Input_Norm:
	
        def __init__(self,rage):

		self.rage = rage 
		self.progress = rage.progress
		self.args = rage.args 
#		self.input = rage.input 
		self.D = rage.data 


	def rpkm(self): 

		
		F,S = self.D.features,self.D.samples 
	
		for s in S: 
			s.notes['total'] = sum(s.cnts.values()) 

		total_factor =  np.percentile(([s.notes['total'] for s in S if s.notes['total']>0]),10) 
		rpkm_output = [] 
		for f in F: 
			rpkm_output.append([0 if s.notes['total'] == 0 else round((f.cnts[i]*total_factor) / s.notes['total'],4) for i,s in enumerate(S)])
	
		rpkm_out = rage_outputs.count_file(self.args).write_row_col_data(F,S,rpkm_output,{'name': self.args.prefix+'_rpkm.cnts'}) 

	def quantile(self):

		F,S = self.D.features,self.D.samples 
		f_vals = dd(int) 
		for s in S: 
			s.notes['q_srt'] = sorted(s.cnts.items(), key = lambda X: X[1]) 
			s_offset = len(F) - len(s.notes['q_srt'])
			for i,(f_idx,f_cnt) in enumerate(s.notes['q_srt']): 
				f_vals[i+s_offset] += f_cnt
		f_norms =  [f_vals[fk]/float(len(S)) for fk,fv in f_vals.items()]
		for s in S: 
			s_len = len(s.notes['q_srt']) 
			s_quantiles = dd(float) 
			s_start_offset = 0 
			s_end_offset   = len(f_norms) - len(s.notes['q_srt']) 
			fI,fC = [s.notes['q_srt'][0][0]],[s.notes['q_srt'][0][1]]
			for (fi,fc) in s.notes['q_srt'][1::]: 
				if fc == fC[-1]: 
					fI.append(fi)
					fC.append(fc) 
				else:
					s_end_offset += len(fI) 
					f_norm_avg = np.mean(f_norms[s_start_offset:s_end_offset]) 
					for f_idx in fI: s_quantiles[f_idx] = f_norm_avg
					fI, fC = [fi], [fc] 
					s_start_offset = s_end_offset
			s_end_offset += len(fI) 
			f_norm_avg = np.mean(f_norms[s_start_offset:s_end_offset]) 
			for f_idx in fI: s_quantiles[f_idx] = f_norm_avg
			s.notes['quantile_values'] = s_quantiles 

		quantile_output   = [[round(s.notes['quantile_values'][f.idx],2) for s in S] for f in F]
		quantile_out = rage_outputs.count_file(self.args).write_row_col_data(F,S,quantile_output,{'name': self.args.prefix+'_quantile.cnts'}) 





        def size_factors(self):


		F,S = self.D.features,self.D.samples 
		F_avg = [sum(f.cnts.values())/float(len(S)) for f in F] 
#		F_avg = [np.power(np.product(f.cnts.values()),1/float(len(f.cnts))) for f in F] 

		s_factors = [] 
#		log_avgs = [math.log(sum(f.cnts)/float(len(S))) for f in F] 
		QT = 98 
		QT_vals = [np.percentile([0 for j in range(len(F)-len(s.cnts))]+s.cnts.values() , QT) for s in S] 
		QT_avg = np.mean(QT_vals) 
#		QT_exp =  [exp(np.log(qt)-np.log(QT_avg)) for qt in QT_vals]
		QT_log =  [np.log(qt)-np.log(QT_avg) for qt in QT_vals]
		QT_offset = [qt-np.mean(QT_log) for qt in QT_log] 
		QT_fin =  [np.exp(qt) for qt in QT_offset]

		w = open(self.args.prefix+'_sizeFactors_qt_'+str(QT)+'.out','w') 
		w.write('%-40s %40s %30s\n' % ('---','s_idx','size_factor'))
		for i,s in enumerate(S): 
			w.write('%-40s %40s %30f\n' % (s.name,i,QT_fin[i]))

			
		size_counts =	[[round(f.cnts[i]*QT_fin[i],3) for i in range(len(S))] for f in F] 
		size_out = rage_outputs.count_file(self.args).write_row_col_data(F,S,size_counts,{'name': self.args.prefix+'_sizeFactors_qt_'+str(QT)+'.cnts'}) 
			







	def downsample(self):

		F,S = self.D.features,self.D.samples 
		if self.args.ds != None:
			self.progress.start_minor('Reading Downsampling Parameters',len(self.D.samples),False)

			if self.args.ds.name.split('.')[-1] == 'vals':

				ds_key = dd(lambda: {}) 
				vkey = {'Impute-Missing': 'impDrop', 'Impute-Relative': 'impRel','Infer-Missing': 'infDrop','Infer-Relative': 'infRel'}
				dkey = {'25000': '25k', '50000': '50k','100000': '100k'}
				for line in self.args.ds: 
					line = line.split() 
					if line[0] == '---': continue 
					
					sample,dstyle = line[0], vkey[line[4]]	
					for i in range(5,len(line),3):
						dtot = dkey[line[i]]
						ds_fail = line[i+1]  
						dval  = int(line[i+2]) 
						if ds_fail == 'True':	ds_key[dtot][sample] = False
						else: 			ds_key[dtot][sample] = dval

				for dtot in ds_key:
					pc = 'Rnorm-'+dstyle+'-'+dtot
					for s in self.D.samples:  				
						self.progress.mark() 
						s.notes[pc] = ds_key[dtot][s.name]
		
					missed = [s.name for s in self.D.samples if s.cnt_total < s.notes[pc]]
					self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled + ('+",".join(missed)+'\n')
					self.run_downsampler(method='notes', param = pc ) 

				
			else:
				ds = {}  
				for line in self.args.ds: 
					line = line.split() 
					if line[0] == '---': continue 
					else: ds[line[0]] = float(line[-1]) 
				for s in self.D.samples:
					self.progress.mark() 
					s.notes['precomp'] = ds[s.name] 	
				missed = [s.name for s in self.D.samples if s.cnt_total < s.notes['precomp']] 
				self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled + ('+",".join(missed)+'\n')
				self.run_downsampler(method='notes', param = 'precomp') 
		else:
			self.progress.start_minor('Selecing Downsampling Parameters',5000,False)
			if self.args.size == 0:
				ds_attempt =  np.percentile([s.cnt_total for s in self.D.samples],3) 
				if ds_attempt > 10000: 
					if ds_attempt < 15000:     ds_attempt    = 10000
					elif ds_attempt < 20000:   ds_attempt  = 15000
					elif ds_attempt < 30000:   ds_attempt  = 20000
					elif ds_attempt < 50000:   ds_attempt  = 25000
					elif ds_attempt < 75000:   ds_attempt  = 50000
					elif ds_attempt < 150000:  ds_attempt = 100000
					elif ds_attempt < 250000:  ds_attempt = 150000
					elif ds_attempt < 500000:  ds_attempt = 250000
					elif ds_attempt < 1000000: ds_attempt = 500000
					else:			   ds_attempt = 1000000
				self.args.size = ds_attempt
			self.progress.notate('Size Parameter Selected: '+str(self.args.size)+'\n')
			missed = [s.name for s in self.D.samples if s.cnt_total < self.args.size] 
			if len(missed) > 100: 	self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled\n')
			elif len(missed) > 0: 	
				self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled + ('+",".join(missed)+'\n')
			self.run_downsampler(method = 'value', param = int(self.args.size)) 


		self.progress.end() 














	def rsample(self,F_TOTALS=[('25k',25000),('50k',50000),('100k',100000)]): 

		F,S = self.D.features,self.D.samples 

		saved_cnts = {} 	
		feature_comps = rage_comps.features(self.rage).get_totals()

		## amplification_key

		mInf = int(feature_comps.minInfer)
		mImp = int(feature_comps.minImpute)

		INFER_TOTALS = [(str(mInf),mInf)]   + [x for i,x in enumerate(F_TOTALS) if (x[1] > mInf and (i==0 or x[1] < mInf*5))]
		IMPUTE_TOTALS = [(str(mImp),mImp)] + [x for i,x in enumerate(F_TOTALS) if (x[1] > mImp and (i==0 or x[1] < mInf*5))]



		for s in S:

			saved_cnts[s.name] = s.cnts.items() 

			for a,b in INFER_TOTALS: 	s.notes['fh_infer_'+a] = b*feature_comps.fraction_key[s.name][0] 
			for a,b in IMPUTE_TOTALS: 	s.notes['fh_impute_'+a] = b*feature_comps.fraction_key[s.name][1] 

			
			


		for a,b in INFER_TOTALS: 
			missed = [s.name for s in self.D.samples if s.cnt_total < s.notes['fh_infer_'+a]]
   			self.progress.start_minor('Preparing Inferred Downsampling for Total '+a,len(self.D.samples),False)
			if len(missed) < 20 and len(missed)>0: self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled + ('+",".join(missed)+'...')
			elif len(missed) >=20: 		       self.progress.notate('Warning: '+str(len(missed))+' samples and will be upsampled...')
			self.run_downsampler(method='notes', param = 'fh_infer_'+a)
		
		for a,b in IMPUTE_TOTALS: 
			missed = [s.name for s in self.D.samples if s.cnt_total < s.notes['fh_impute_'+a]]
   			self.progress.start_minor('Preparing Imputed Downsampling for Total '+a,len(self.D.samples),False)
			if len(missed) < 20 and len(missed)>0: self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled + ('+",".join(missed)+'...')
			elif len(missed) >=20: 		       self.progress.notate('Warning: '+str(len(missed))+' samples and will be upsampled...')
			self.run_downsampler(method='notes', param = 'fh_impute_'+a)

		self.progress.end() 


		




	
	def gad(self): 
		F,S = self.D.features,self.D.samples 
		mf1 = {f.name: sum([float(b) / self.D.samples[a].cnt_total for (a,b) in f.cnts.items()]) / len(f.cnts) for f in self.D.features}
		mf2 = {f.name: sum([float(b) / self.D.samples[a].cnt_total for (a,b) in f.cnts.items()]) / len(self.D.samples) for f in self.D.features}

		mt1 = sum(mf1.values())
		mt2 = sum(mf2.values())
		mt1_cands, mt2_cands = [], [] 
		for s in self.D.samples:
			frac1 = sum([mf1[self.D.features[fi].name] for fi in s.cnts.keys()]) / mt1 
			frac2 = sum([mf2[self.D.features[fi].name] for fi in s.cnts.keys()]) / mt2 
			s.notes['frac1'] = frac1
			s.notes['frac2'] = frac2
			mt1_cands.append((s.cnt_total / frac1,s.idx)) 
			mt2_cands.append((s.cnt_total / frac2,s.idx)) 
		
#		mt1_choice = sorted(mt1_cands)[25][0]
#		mt2_choice = sorted(mt2_cands)[25][0]
		mt1_choice = sorted(mt1_cands)[25][0]
		mt2_choice = sorted(mt2_cands)[25][0]
		for s in self.D.samples: 
			s.notes['gadC']  =  mt1_choice * s.notes['frac1'] 
			s.notes['gadA']  =  mt2_choice * s.notes['frac2'] 

		self.run_downsampler(method='notes', param = 'gadC') 
		self.run_downsampler(method='notes', param = 'gadA') 



	def run_downsampler(self,method, param ):

   		self.progress.start_minor('Running Downsampler..' ,len(self.D.features),False)
		F,S = self.D.features,self.D.samples 


		if method == 'notes': 	
			ds_vals = {s: s.notes[param] for s in self.D.samples if s.notes[param] != False}
			name = self.args.prefix+'_'+param+'.cnts'
	

		elif method == 'value':

			ds_vals = {s: param for s in self.D.samples}
			name = self.args.prefix+'_ds'+str(param)+'.cnts'
	

		for s,dsv in ds_vals.items():
			self.progress.mark() 
			f_idx, f_rand, f_rep, f_iter, f_key  =  [[(0,0),'INIT']], [], 0, 0, dd(int)
			for f,c in s.cnts.items(): f_idx.append([(f_idx[-1][0][1],f_idx[-1][0][1]+int(c)),f]) 
			while s.cnt_total < dsv: 
				dsv -= s.cnt_total
				f_rep         += 1 

			for r in sorted(random.sample(range(int(s.cnt_total)), int(dsv))):
				while r > f_idx[f_iter][0][1]: f_iter+=1 
				f_key[f_idx[f_iter][1]] += 1
		
			for f,c in s.cnts.items(): 
				fv = f_key[f] + (f_rep*c) 
				if fv > 0: s.norms['ds'][f] = fv 


		if method == 'value':	
			f_output = [[int(s.norms['ds'][f.idx]) for s in self.D.samples] for f in self.D.features]
			ds_out = rage_outputs.count_file(self.args).write_row_col_data(self.D.features,self.D.samples,f_output,{'name': name})
		else:
			S = [s for s in self.D.samples if s.notes[param] != False] 
			f_output = [[int(s.norms['ds'][f.idx]) for s in S] for f in self.D.features]
			ds_out = rage_outputs.count_file(self.args).write_row_col_data(self.D.features,S,f_output,{'name': name})

			#f_output = [[int(s.norms['ds'][f.idx]) for s in self.D.samples] for f in self.D.features]
			#ds_out = rage_outputs.count_file(self.args).write_row_col_data(self.D.features,self.D.samples,f_output,{'name': name})









	def rank_norm(self):

		F,S = self.D.features,self.D.samples 
		total = len(self.D.features) 
		f_ranks, f_data = [] , [] 
   		self.progress.start_minor('Running Rank Normalization..' ,len(S)+len(F),False)
		for s in S: 
			self.progress.mark() 
			f_srt = sorted(s.cnts.items(), key = lambda X: X[1],reverse=True) 
			f_ranks.append({f[0]: i+1 for i,f in enumerate(f_srt)})
			


		for i,f in enumerate(F):
			self.progress.mark() 
			f_data.append([fr[f.idx] if f.idx in fr else total for fr in f_ranks])

		rank_out = rage_outputs.count_file(self.args).write_row_col_data(F,S,f_data,{'name': self.args.prefix+'_rank.cnts'}) 
		self.progress.end()




















	def top_norm(self):

		if  self.args.size == 0: top_rank = 1000 
		else: 			 top_rank = self.args.size


		F,S = self.D.features,self.D.samples 
		total = len(self.D.features) 
		f_ranks, f_data = [] , [] 
   		self.progress.start_minor('Running Rank Normalization..' ,len(S)+len(F),False)

		s_lens = [len(s.cnts) for s in S] 
		q_key = dd(list) 
		for s in S: 
			self.progress.mark() 
			f_srt = sorted(s.cnts.items(), key = lambda X: X[1],reverse=True) 

			s_key = {}
			for i,f in enumerate(f_srt):
				q_key[i].append(f[1]) 
				s_key[f[0]] = i
				if i >= top_rank: break 
			f_ranks.append(s_key) 			
		s_len = len(S) 
		mean_key = {i: round(np.mean(q_key[i]),2) for i in range(len(q_key))}


		c_data, q_data = [],[] 

		for i,f in enumerate(F):
			s_c,s_q = [],[] 
			for j,fr in enumerate(f_ranks):
				if f.idx in fr: 
					s_q.append(round(mean_key[fr[f.idx]],2))
					s_c.append(int(S[j].cnts[f.idx]))

				else:
					s_c.append(0) 
					s_q.append(0) 
			c_data.append(s_c) 
			q_data.append(s_q) 

		q_out = rage_outputs.count_file(self.args).write_row_col_data(F,S,q_data,{'name': self.args.prefix+'_top_'+str(top_rank)+'_qt.cnts'}) 
		c_out = rage_outputs.count_file(self.args).write_row_col_data(F,S,c_data,{'name': self.args.prefix+'_top_'+str(top_rank)+'_raw.cnts'}) 


		self.progress.end()





