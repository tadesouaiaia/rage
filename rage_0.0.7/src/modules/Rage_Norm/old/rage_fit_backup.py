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








class Input_Norm:
	
        def __init__(self,rage):

		self.progress = rage.progress
		self.args = rage.args 
		self.input = rage.input 



	def rpkm(self): 

		d = self.input 
		
		F,S = self.input.features,self.input.samples 
	
		for s in S: 
			s.notes['total'] = sum(s.cnts.values()) 

		total_factor =  np.percentile(([s.notes['total'] for s in S if s.notes['total']>0]),10) 
		rpkm_output = [] 
		for f in F: 
			rpkm_output.append([0 if s.notes['total'] == 0 else round((f.cnts[i]*total_factor) / s.notes['total'],4) for i,s in enumerate(S)])
	
		rpkm_out = rage_outputs.count_file(self.args).write_row_col_data(F,S,rpkm_output,{'name': self.args.prefix+'_rpkm.cnts'}) 

	def quantile(self):

		d = self.input 
		F,S = self.input.features,self.input.samples 

		#pass_idx = [f.idx for f in F if 100 * (len(f.cnts) / float(len(S))) >= self.args.rankObservations]

		#pass_scores = dd(int) 
		f_vals = dd(int) 
		for s in S: 
			s.notes['q_srt'] = sorted(s.cnts.items(), key = lambda X: X[1]) 
			s_offset = len(F) - len(s.notes['q_srt'])
			for i,(f_idx,f_cnt) in enumerate(s.notes['q_srt']): 
				f_vals[i+s_offset] += f_cnt


			
#			s.notes['quantiles'] = sorted([(s.cnts[idx],idx) for idx in pass_idx])
#			for i in range(len(pass_idx)):	pass_scores[i] += s.notes['quantiles'][i][0]
#		pass_values = [log((x/float(len(S)))+1) for x in pass_scores.values()]
#		pass_values = [x/float(len(S)) for x in pass_scores.values()]
#		pass_values = [x for x in pass_scores.values()]
#		pass_total  = float(len(pass_values)) 
#		f_norms =  {fk: f_vals[fk]/float(len(S)) for fk,fv in f_vals.items()}
#		f_norms =  [f_vals[fk]/float(len(S)) for fk,fv in f_vals.items()]
#		f_norms =  [f_vals[fk]/float(len(S)) for fk,fv in f_vals.items()]
		f_norms =  [f_vals[fk]/float(len(S)) for fk,fv in f_vals.items()]
#		f_norms =  [f_vals[fk]/float(len(S)) for fk in range(len(F))] # in f_vals.items()]
		
#		print len(F), len(f_norms) 


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






	def downsample(self):


   		self.progress.start_minor('Selecing Downsampling Parameters',5000,False)

		d = self.input 
		if self.args.size == 0:
			ds_attempt =  np.percentile([s.cnt_total for s in d.samples],3) 
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


		missed = [s.name for s in d.samples if s.cnt_total < self.args.size] 
		if len(missed) > 100: 	self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled\n')

    		elif len(missed) > 0: 	
			self.progress.notate('Warning: '+str(len(missed))+' samples will be upsampled + ('+",".join(missed)+'\n')
			#sys.exit() 



		self.run_downsampler(method = 'value', param = int(self.args.size)) 


	def rsample(self): 

		d= self.input
		r_key = {} 
		expected_totals = [] 
   		self.progress.start_minor('Selecting DropOut Aware Feature Parameters',len(d.features)/10000.0,True)

		w = open(self.args.prefix+'_medRatios.out','w') 

	
		w.write('%-40s %40s %5s %5s | %10s %25s\n' % ('---','---','idx1','idx2','matches','ratio'))
		for fi in range(len(d.features)): 
			self.progress.mark(0.75) 
			f_key = dd(list) 
			for si,ci in d.features[fi].cnts.items():
				for fj,cj in [st for st in d.samples[si].cnts.items() if st[0] > fi]: 
					f_key[(fi,fj)].append(ci/float(cj))
			for fp,fl in f_key.items():
				fLen,fMed = len(fl), np.median(fl)  
				r_key[fp] = fMed
				w.write('%-40s %40s %5d %5d | %10d %25f\n' % (d.features[fp[0]].name,d.features[fp[1]].name,fp[0],fp[1],fLen,fMed))
			


#		sys.exit() 
#		print r_key[2]
#


#		for i in range(len(d.features)):
#			self.progress.mark()
#
#			for j in range(i+1,len(d.features)):
#				r_key[(i,j)] = np.median([s.cnts[i]/float(s.cnts[j]) for s in d.samples if i in s.cnts and j in s.cnts])
					



#		print  r_key[2] 

#		sys.exit() 




   		self.progress.start_minor('Inputing Missing Sample Data',len(d.samples),False)
		for s in d.samples: 
			self.progress.mark()
			missing = [i for i in range(len(d.features)) if i not in s.cnts]
			missing_cnts = [] 
			for m in missing: 
				m_preds = [s.cnts[i]/r_key[(i,m)] if i<m else s.cnts[i]*r_key[(m,i)] for i in s.cnts]
				missing_cnts.append(np.mean(m_preds)) 

			expected_total = sum(s.cnts.values())+sum(missing_cnts)
			expected_totals.append([expected_total,sum(s.cnts.values()),s.idx]) 

		t_srt = [[t[1]/float(t[0]),t[1],t[2]] for t in expected_totals]
		t_min =  min([t[1]/t[0] for t in t_srt]) 



   		self.progress.start_minor('Preparing Downsampler',len(d.samples),False)
		for t in t_srt:
			#aj = t[0]*t_min 
			#s = d.samples[t[-1]] 
			d.samples[t[-1]].notes['ratio'] = t[0]*t_min




		self.run_downsampler(method='notes', param = 'ratio')

		

	
	def gad(self): 
		d = self.input 
		mf1 = {f.name: sum([float(b) / d.samples[a].cnt_total for (a,b) in f.cnts.items()]) / len(f.cnts) for f in d.features}
		mf2 = {f.name: sum([float(b) / d.samples[a].cnt_total for (a,b) in f.cnts.items()]) / len(d.samples) for f in d.features}

		mt1 = sum(mf1.values())
		mt2 = sum(mf2.values())
		mt1_cands, mt2_cands = [], [] 
		for s in d.samples:
			frac1 = sum([mf1[d.features[fi].name] for fi in s.cnts.keys()]) / mt1 
			frac2 = sum([mf2[d.features[fi].name] for fi in s.cnts.keys()]) / mt2 
			s.notes['frac1'] = frac1
			s.notes['frac2'] = frac2
			mt1_cands.append((s.cnt_total / frac1,s.idx)) 
			mt2_cands.append((s.cnt_total / frac2,s.idx)) 
		
#		mt1_choice = sorted(mt1_cands)[25][0]
#		mt2_choice = sorted(mt2_cands)[25][0]
		mt1_choice = sorted(mt1_cands)[25][0]
		mt2_choice = sorted(mt2_cands)[25][0]
		for s in d.samples: 
			s.notes['gadC']  =  mt1_choice * s.notes['frac1'] 
			s.notes['gadA']  =  mt2_choice * s.notes['frac2'] 

		self.run_downsampler(method='notes', param = 'gadC') 
		self.run_downsampler(method='notes', param = 'gadA') 



	def run_downsampler(self,method, param ):

   		self.progress.start_minor('Running Downsampler..' ,len(self.input.features),False)
		d = self.input 

		if method == 'notes': 	
			ds_vals = {s: s.notes[param] for s in d.samples}
			name = self.args.prefix+'_'+param+'.cnts'
		elif method == 'value': 
			ds_vals = {s: param for s in d.samples}
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
			f_output = [[int(s.norms['ds'][f.idx]) for s in d.samples] for f in d.features]
		else:
			f_output = [[int(s.norms['ds'][f.idx]) for s in d.samples] for f in d.features]



		ds_out = rage_outputs.count_file(self.args).write_row_col_data(d.features,d.samples,f_output,{'name': name})
		self.progress.end() 





	def fit_binary_dists(self,CUTOFF=4):

		for f in self.input.features: 
			vals,logV, dZ = [int(x) for x in f.cnts.values()], [log(v+1.0) for v in f.cnts.values()],[0 for i in range(self.input.samples.len - len(f.cnts.values()))]
			if len(vals) < CUTOFF: continue 
			if len(vals) > 10: continue 
				
			val_key = {'RAW-NZ': vals, 'RAW-WZ': vals+dZ,  'LOG-NZ': logV, 'LOG-WZ':  logV+dZ}
			
			for val_type,vals in val_key.items(): 
					
				vLen, v10, vMean = len(vals), int(len(vals)*0.10), np.mean(vals) 
				if val_type.split('-')[0] == 'LOG': continue 
				else:
					r      = stats.poisson.rvs(vMean, size=len(vals))
					both = sorted(cc(vals+r).items()) 					
					
					bins, span, sT = [],[], 0   	
					for v,c in both: 
						span.append(v) 
						sT += c 
						if sT > v10:
							bins.append((span[0],span[-1]))
							span,sT = [v+1],0 
					
					print bins 	


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














































