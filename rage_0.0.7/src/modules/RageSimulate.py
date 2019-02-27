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

from operator import mul

from scipy.signal import savgol_filter as svgf 
from math import exp
from math import factorial 


from Rage_IO import rage_outputs
from Rage_Plots import rage_subplots
from Rage_Plots import  rage_dimplots as dplot 

from Rage_Transforms import rage_KDE
from Rage_Transforms import rage_DR
from Rage_Regression import rage_regression 
from Rage_Filters    import rage_filters 
#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 



def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))

def compare_principal_components(pca1,pca2):

	print 'yo' 


class Simulate:
        def __init__(self,rage):

		self.rage = rage 
		self.progress = rage.progress
		self.options = rage.args 



	def run(self): 


		if self.rage.data: 
			self.rage.data.filter_samples_by_attributes(self.options.predictors)
			if len(self.options.predictors) > 0: 
				self.simulate_dex()

		sys.exit() 
		if self.options.command == 'samples': 
			print 'ok' 
		elif self.options.command == 'dex':
			if not self.rage.data: 
				print 'no data'
				print 'need data for dex' 
				sys.exit()  

			self.simulate_dex() 
		else:
			print 'huh'	
			sys.exit()



	def get_betas(self,preds,inferred_preds): 
		
		if 'intercept' in preds: betas = ['B_o'] 
		else:		       	 betas = []
		for p in [P for P in  preds if P != 'intercept']:
			betas.append('B_'+str(len(betas))+'('+p+')')
		return betas	


	def simulate_dex(self):
		D = self.rage.data
#		print D.samples
#		print D.features
		 

		dex_key, dex_data, dex_set = dd(lambda: dd(bool)), [], set([]) 
		suffix = 'dexsim_'+"_".join(self.options.predictors)
		for p in self.options.predictors: 
			p_vals = [s.attributes[p] for s in D.samples] 
			D.samples.attribute_class['SIM:'+p] = D.samples.attribute_class[p] 
			if self.options.randomize: shuffle(p_vals)
			for j,s in enumerate(D.samples): 
				s.attributes['SIM:'+p] = p_vals[j] 
		

		
		self.options.predictors = ['SIM:'+p for p in self.options.predictors]
  
		X = D.create_sample_predictor_variables(self.options.predictors) 
		preds = [p for p in D.predictors if p != 'intercept'] 


		for f in D.features:
			cnts = [f.cnts[s.idx] for s in D.samples]
			log_cnts = [log(1.0+f.cnts[s.idx]) for s in D.samples]
			s_changes = [[] for s in D.samples]
			for p in preds: 
				if random.random() < self.options.dex_rate:
			#		effect_mu = 0.15
			#		effect_var = 1.0

					if len(s_changes[0])!=0: continue

					effect_mu, effect_var = random.choice([(0.25,0.05),(0.5,0.1),(2.0,0.5)])
			#		sys.exit() 
 
					dex_key[f.name][p] = (effect_mu,effect_var)
					dex_set.add(p) 
					if len(p.split('=')) == 1: 
						sv  = scale_vals([s.attributes[p] for s in D.samples])
						for j,s in enumerate(D.samples):
							s_changes[j].append(np.random.normal(effect_mu,effect_var,1)[0] * sv[j])
					else:
						for j,s in enumerate(D.samples):
							if s.attributes[p.split('=')[0]] == p.split('=')[-1]:
								s_changes[j].append(np.random.normal(effect_mu,effect_var,1)[0])

					


			if len(dex_key[f.name]) == 0: dex_data.append([f.cnts[s.idx] for s in D.samples])
			else:			      dex_data.append([int(math.exp(log(f.cnts[s.idx]+1.0) + sum(s_changes[j]))-1.0) for j,s in enumerate(D.samples)])


		sim_out = rage_outputs.sim_counts(self.options)
		sim_out.write_row_col_data(D.features,D.samples,dex_data,{'suffix': suffix+'.cnts'})
 		sim_out.write_feature_key(dex_key,sorted(list(dex_set)),{'suffix': suffix+'.features.key'})
		sim_out.write_sample_key(D.samples,{'suffix': suffix+'.samples.key'}) 



  


	def evaluate_predictor(self,interest,simulations=4):
		D = self.rage.data
		X = self.rage.data.create_sample_predictor_variables(self.options.predictors,intercept=True) 

		X_sim = [self.rage.data.simulate_sample_predictor_variables(self.options.predictors,interest) for i in range(simulations)]

		dim_red = rage_DR.DR(self.options,self.progress,len(D.samples)) 
		#seaborn.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'w'})
		#res, subplot = dd(lambda: {}), rage_subplots.subplot(3,2,options) 
		### INITIAL PCA ### 
		pca_init = dim_red.run_pca(D.matrix())['run']

		### MODEL RESULT ### 
		model_result = [self.regress(y,X) for y in [[log(s.cnts[f.idx]+1.0) for s in D.samples] for f in D.features]]
		pca_resid =    dim_red.run_pca(np.matrix([m[0].resid for m in  model_result]).getT())['run']
		model_dex_top = [m[1][interest][0] for m in model_result]

		betas = self.get_betas(D.predictors,D.inferred_predictors) 

		print D.predictors,D.inferred_predictors
		print 'initial variance pca',sum(pca_init.explained_variance_) 



		print 'model variance pca',sum(pca_resid.explained_variance_) 
		print 'model rsquared',    np.mean([m[0].rsquared for m in model_result]) 
		print 'model radj',        np.mean([m[0].rsquared_adj for m in model_result]) 

		print 'model_sig_features_05', len([1 for m in model_dex_top if m[1][0] < 0.05])

		for i,Xi in enumerate(X_sim):
			sim_result = [self.regress(y,Xi) for y in [[log(s.cnts[f.idx]+1.0) for s in D.samples] for f in D.features]]
			pca_sim =    dim_red.run_pca(np.matrix([m[0].resid for m in  sim_result]).getT())['run']
			sim_dex_top = [m[1][interest][0] for m in sim_result]
			
			print 'sim '+str(i+1) 
			print 'sim variance pca',sum(pca_sim.explained_variance_) 
			print 'sim rsquared',    np.mean([m[0].rsquared for m in sim_result]) 
			print 'sim radj',        np.mean([m[0].rsquared_adj for m in sim_result]) 
			print 'sim_sig_features_05', len([1 for m in sim_dex_top if m[1][0] < 0.05])
			print "" 

#		print pca_resid.noise_variance_
#		print pca_resid.explained_variance_
#		print pca_post.noise_variance_
#		print sum(pca_post.explained_variance_ratio_)
		sys.exit() 



	def regress_ols(self,y,X):



		model = sm.OLS(y,X).fit()
		p_out = dd(lambda: {}) 
		r_out = {} 
		for p in self.rage.data.inferred_predictors: p_out[p.split('=')[0]][p.split('=')[1]] = (1,0) 
		for pv,bw,c in zip(model.pvalues,model.params,self.rage.data.predictors):	p_out[c.split('=')[0]][c.split('=')[-1]] = (pv,bw)

		for a,b in p_out.items():	r_out[a] = sorted(b.items(),key=lambda loc: loc[1][0])
		return model,r_out  
		print model.fvalue
		print model.rsquared, model.f_pvalue, model.aic, model.bic, model.rsquared_adj, model.ssr
		print model.rsquared
		print model.aic
		print model.bic 
		print model.rsquared
		print model.ssr 

		print model.params
		print model.pvalues
		mrf =  model.df_model
		mdf  = model.df_resid 
		mf = np.sqrt(model.fvalue) 









	def rungo(self): 
		D = self.rage.data
		self.rage.data.create_sample_predictor_variables(self.options.predictors) 



		X = [np.array([s.preds[x] for x in D.predictors]) for s in D.samples]


#		print D.predictors

 

		for f in D.features: 
			y =  [log(s.cnts[f.idx]+1.0) for s in D.samples]
			model = sm.OLS(y,X).fit()
#			for v in vars(model): print v
#			print model.fvalue
#			print model.rsquared, model.f_pvalue, model.aic, model.bic, model.rsquared_adj, model.ssr
#			print D.predictors
#			print model.pvalues
#			print model.params
			#gamma_model = sm.GLM(data.endog, data.exog,family=sm.families.Gamma())
#			gamma_model = sm.GLM(y, X,family=sm.families.Gamma())
#			g_res = gamma_model.fit() 
#			print g_res.summary() 
#			print model.summary() 
			mrf =  model.df_model
			mdf  = model.df_resid 
#			print 0.1/0.9 
#			print np.sqrt(model.fvalue),'yo',0.9/0.1,9*0.1, 0.1/0.9
			mf = np.sqrt(model.fvalue) 
#			print mf, 0.1/0.9
#			print mf / (mf+1.0) 

			print f.name, model.rsquared, model.f_pvalue, model.params			


#			print model.rsquared / (1.0 - model.rsquared) , mf 
			#
			# 

#			print smp.FTestPower().solve_power(effect_size=mf, df_num=mrf, df_denom=mdf, alpha=0.0005),'yo'
#			print smp.FTestPower().solve_power(effect_size=np.sqrt(0.1/(1-0.1)), df_num=89, df_denom=5, alpha=0.05)






	def run2(self):
		R = self.rage 
		sys.exit() 
		self.prepare_predictors(predictors) 
		print R.data.samples.len 
		




		sys.exit() 


		if R.args.command == 'samples':
			R.progress.start_major('SampleSummary')

			if R.args.raw: 
				R.progress.start_minor('Calculating Summary Stats',R.data.samples.len)
				sample_stats = summary_hists(R.data.samples,R.data.features,R.args,R.progress)  
				rage_outputs.column_stats(R.args).write(sample_stats,R.data.samples,{'suffix': 'samplestats','width': 20})
				sys.exit() 


			else:
				
				#sample_dists = summary_dists(R.data.samples,R.data.features,R.args,R.progress) 
				R.progress.start_minor('Summarizing Dimensional Reduction',R.data.samples.len)
				pca, tsne, kca = make_dr_plots(R) 	
				rage_outputs.column_coefs(R.args).write(pca['coefs'],R.data.features,{'suffix': 'PCAcoeffs.features.out','width': 20})
				rage_outputs.dr_matrix(R.args).write(kca['run'].X_transformed_fit_,R.data.samples,{'suffix': 'kcaFit.samples.out'})
				for dr in [pca,tsne,kca]:	rage_outputs.dr_pts(R.args).write(dr['pts'],R.data.samples,{'suffix': dr['type']+'.pts.out'})


		elif R.args.command in ['genes','features']:
			R.progress.start_major('FeatureSummary')
			if R.args.raw: 
				R.progress.start_minor('Calculating Summary Stats',R.data.features.len)
				feature_stats = summary_hists(R.data.features,R.data.samples,R.args,R.progress)  
				rage_outputs.column_stats(R.args).write(feature_stats,R.data.features,{'suffix': 'featurestats.out','width': 20})

			else:
				R.progress.start_minor('Summarizing Dimensional Reduction',R.data.features.len)
				pca, tsne, kca = make_dr_plots(R,'features') 
				rage_outputs.column_coefs(R.args).write(pca['coefs'],R.data.samples,{'suffix': 'PCAcoeffs.samples.out','width': 20})
				rage_outputs.dr_matrix(R.args).write(kca['run'].X_transformed_fit_,R.data.features,{'suffix': 'kcaFit.features.out'})
				for dr in [pca,tsne,kca]:	rage_outputs.dr_pts(R.args).write(dr['pts'],R.data.features,{'suffix': dr['type']+'.featurePts.out'})


				 	 
					


		sys.exit() 



def prepare_predictor_data(X,Y,options,progress,X_NAME='SAMPLES'):


















			#if self.args.pca: sample_summary.make_pca_and_tsne_plots() 
			#sample_summary.summarize_sample_stats() 
	#		sample_summary.summarize_sample_pts()
		#	sample_summary.summarize_sample_dists()
#			sample_summary.summarize_sample_pairs() 

		sys.exit() 
		if self.args.command == 'genes': 
			#feature_summary = rage_summarize_features.Overview(self)
			


			for f in self.input.features: 
				fc = sorted(f.cnts.values()) 
				if len(fc) < 10: continue 
				q1,q3 = int(len(fc)*0.25),int(len(fc)*0.75)
				iqr = fc[q3]-fc[q1] 
				cM,cV = np.mean(fc),np.var(fc) 

				bO = [x for x in fc if x > fc[q3]+(2*iqr)]

				s3 = [x for x in fc if x > cM+(3*(cV**0.5))]


				print f.name,len(fc),cM,cV,max(fc),'basic/s3',len(bO),len(s3)

			sys.exit() 

			if self.args.fitdists: 	
				self.dists = rage_summarize_dists.Dists(self.args,self.input.features, self.input.samples.len) 
				self.dists.fit_binary() 
				feature_summary.fit_binary_dists() 
				feature_summary.fit_dists()
		sys.exit() 




def summary_hists(X,Y,options,progress,X_NAME='SAMPLES'):

	#seaborn.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'lightgray'})
	seaborn.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'w'})
	res, subplot = dd(lambda: {}), rage_subplots.subplot(3,2,options) 
	cMax = float(sum([y.len for y in Y]))
	for x in X: 
		progress.mark() 
		res['#ComparisonIndex'][x] = sum([Y[y].len for y in x.cnts.keys()]) / cMax
		ordered_logs = sorted([log(1.0+c) for c in x.cnts.values()],reverse=True)	
		res['#Total'][x], res['#Obs'][x],  halfE,iX,k = x.cnt_total, len(x.cnts),  sum(ordered_logs) * 0.5, 0, -1 		
		res['#Obs_AboveMean'][x] = len([l for l in ordered_logs if l > np.mean(ordered_logs)])/float(len(ordered_logs))
		while iX < halfE:
			k+=1;	iX+=ordered_logs[k] 
		res['%Obs_HalfDepth'][x] = k / float(len(x.cnts))
		res['CoeffVar'][x] = coVar(ordered_logs) 

	subplot.add_hist(res['#Total']).update({'xlab':'reads','ylab': 'occurences','title': 'Depth'})	
	subplot.add_hist(res['#Obs']).update({'xlab':'observations','ylab': 'occurences','title': 'Complexity'})	
	subplot.add_hist(res['#Obs_AboveMean']).update({'xlab':'%','ylab': 'occurences','title': '% obs above mean'})
	subplot.add_hist(res['CoeffVar']).update({'xlab':'CV','ylab': 'occurences','title': 'Coefficient of Variation (Log Space)'})
	subplot.add_hist(res['%Obs_HalfDepth']).update({'xlab':'%Obs','ylab': 'occurences','title': '% Obs For 50% Read Depth (Log Space)'})
	subplot.add_hist(res['#ComparisonIndex']).update({'xlab':'%Comparisons','ylab': 'occurences','title': 'Comparison Index'})
	plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.90,wspace=0.1,hspace=0.40)
	if X.label == 'samples': 	subplot.save(options.prefix+'_sample_summary.png',{'title': 'Sample Summary Histograms'}) 
	elif X.label == 'features':	subplot.save(options.prefix+'_feature_summary.png',{'title': 'Feature Summary Histograms'}) 
	return res 	





def make_dr_plots(R,choice='samples'):

	if choice == 'samples': 
		r_members = R.data.samples 
		r_matrix = R.data.matrix('log') 	
		out_name = R.args.prefix+'_samples_'
	else:
		r_members = R.data.features
		r_matrix = R.data.matrix('log').getT()
		out_name = R.args.prefix+'_features_'

	dr = rage_DR.DR(R.args,R.progress)
	pca_run = dr.run_pca(r_matrix) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,pca_run,{'title':'PCA','out': out_name+'pca.pdf'})
	tsne_run = dr.run_tsne() 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,tsne_run,{'title':'TSNE','out': out_name+'tsne.pdf'})  
	kca_run = dr.run_kca(r_matrix) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA','out': out_name+'kca.pdf','zoom': True})  
	return pca_run, tsne_run, kca_run 








def simsample_items(L,size=200):

	L_new,L_now = [], [int(round(10.0 * v,0)) for v in L]
	shuffle(L_now)
	while len(L_now) + len(L_new) < size: L_new.extend(L_now) 
	L_new.extend(random.sample(L_now,size-len(L_new)))
	Lc= cc(L_new)
	return [Lc[a] if a in Lc else 0 for a in range(0,11)]





def summary_dists(X,Y,options,progress,X_NAME='SAMPLES'):
	seaborn.set(rc={'axes.facecolor':'lightpink', 'figure.facecolor':'lightgray'})
	progress.start_major('Plotting Distribution Densities',len(X))
	kde = rage_KDE.samples(0.3) 
	f_num,subplot =  1, rage_subplots.subplot(6,2,options) 
	LOG=True

	dr = rage_DR.DR(options,progress) 
	y_vals = scale_vals([log((1.0+sum(y.cnts.values())) ) for y in Y])
	x1,y1 = kde.run(y_vals)
	subplot.add_lines(x1,y1,None,None,'black').update({'title': 'Global Distribution'})





	iter_data = [] 
	for x in X:
		progress.mark()
		non_zeros = scale_vals([log(v+1.0) for v in x.cnts.values()]+[0.0])
		all_vals = scale_vals([0 for v in range(Y.len-(1+len(non_zeros)))] + [log(v+1.0) for v in x.cnts.values()])
		nz = simsample_items(non_zeros) 
		sz = simsample_items(all_vals)
		x.notes['iter'] = [non_zeros,all_vals]


		iter_data.append([non_zeros,all_vals,nz,sz]) 
	r_matrix =  np.matrix([it[2] for it in iter_data])
	pca_run = dr.run_pca(r_matrix)
	kmean_run = dr.run_kmeans(r_matrix) 
	subplot.add_pca_data(pca_run['pts'],{'title': 'PCA on binned distribution values'})#.update({'clear_axes': True}) 
	for i in range(len(kmean_run['labels'][0])):
		X[i].notes['km'] = kmean_run['labels'][0][i] 
		if X[i].notes['km'] == 0:
			subplot.ax.scatter(pca_run['pts'][i][0],pca_run['pts'][i][1],color='yellow')
	subplot.update({'clear_axes': True})

	X0,X1 = [x for x in X if x.notes['km'] == 0], [x for x in X if x.notes['km'] == 1]

	for x1,x2 in zip(X0,X1):
	
		if x2.name in ['EB321','EB1015']: continue
		nz1,az1 = x1.notes['iter'] 
		nz2,az2 = x2.notes['iter'] 

		a1,b1 = kde.run(az1)
		a2,b2 = kde.run(nz1)
		
		subplot.add_lines(a1,b1,None,None,'black')
		subplot.add_lines(a2,b2,None,None,'cyan')
		subplot.update({'clear_axes': True,'title': x1.name}) 

		a1,b1 = kde.run(az2)
		a2,b2 = kde.run(nz2)
		
		subplot.add_lines(a1,b1,None,None,'black')
		subplot.add_lines(a2,b2,None,None,'cyan')
		subplot.update({'clear_axes': True,'title': x2.name}) 
		f_num += 1
		if not subplot.update or f_num > 15: 
			break
		

	plt.subplots_adjust(left=0.07, bottom=0.01, right=0.93, top=0.95,wspace=0.2,hspace=0.6)
	subplot.save(options.prefix+'fig_dists'+str(f_num)+'.png',{'title': 'Dual Dists: '})
	progress.end()

		




















		
def make_pca_and_tsne_plots(self):

	seaborn.set(rc={'axes.facecolor':'black', 'figure.facecolor':'cornflowerblue'})
	my_sizes = scale_vals([len(s.cnts.keys()) for s in self.input.samples],20,55)
	self.progress.start_subtopic('Calculating PCA/TSNE','',0)
	data_matrix = self.input.data_matrix('log')
	dr = rage_DR.DR(self.args,self.progress).set_matrix(data_matrix)
	dr.run_pca().run_kca().run_tsne().run_ica() 
	subplot = rage_subplots.subplot(2,2,self.args)
	subplot.add_legend(self.color_key.keys(),self.color_key.values())
	subplot.add_pca_data(dr.pca_pts,{'vars': dr.pca_vars,'title': 'PCA','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True}) 
	subplot.add_pca_data(dr.kca_pts,{'type': 'kca', 'title': 'KCA','colors':self.color_labels,'zoom': True,'sizes': my_sizes}).update({'clear_axes': True}) 
	subplot.add_pca_data(dr.ica_pts,{'type': 'ica', 'title': 'ICA','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True}) 
	subplot.add_pca_data(dr.tsne_pts,{'type': 'tsne','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True}) 
	#subplot.add_legend(self.color_key.keys(),self.color_key.values())
	subplot.save(self.args.prefix+'_dimred.png',{}) 
	self.progress.finish_subtopic() 





















































class SampleSummary:
	
        def __init__(self,summary):

		self.progress = summary.progress
		self.args = summary.args 
		self.input = summary.input 

		self.sLen = len(self.input.samples) 
		self.fLen = len(self.input.features) 	
		
		self.create_label_key() 

		self.get_feature_order()	


	def summarize_sample_stats(self):

		seaborn.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'lightgray'})
		self.progress.start_subtopic('Calculating Summary Stats','',self.sLen)
		res = dd(lambda: {}) 
		subplot = rage_subplots.subplot(3,2,self.args)  
		for s in self.input.samples: 	
			self.progress.mark_subtopic()
			ordered_logs = sorted([log(1.0+c) for c in s.cnts.values()],reverse=True)	
 			res['#Reads'][s] = s.cnt_total 
			halfE,iX,k = sum(ordered_logs)*0.5,0,-1
			res['#Observed_Genes'][s] = len(ordered_logs) 
			res['#Genes_Above_Mean'][s] = len([x for x in ordered_logs if x > np.mean(ordered_logs)])/float(len(ordered_logs))
			while iX < halfE:
				k+=1;	iX+=ordered_logs[k] 
			res['%Genes_Required_For_HalfDepth'][s] = k / float(len(ordered_logs))
			res['CoeffVar'][s] = coVar(ordered_logs) 
			res['#topVals'][s] = 0 	
	
		for f in self.input.features:
			for a,b in sorted([(b,a) for (a,b) in f.cnts.items()])[-5::]:
				res['#topVals'][self.input.samples[b]]+=1

		subplot.add_hist(res['#Reads'].values()).update({'xlab':'reads per sample','ylab': 'occurences','title': 'Depth'})	
		subplot.add_hist(res['#Observed_Genes'].values()).update({'xlab':'genes per sample','ylab': 'occurences','title': 'Library Complexity'})	
		subplot.add_hist(res['#Genes_Above_Mean'].values()).update({'xlab':'%','ylab': 'occurences','title': '% genes above mean'})
		subplot.add_hist(res['%Genes_Required_For_HalfDepth'].values()).update({'xlab':'%Obs Genes','ylab': 'occurences','title': '% Genes Required For 50% Read Depth (Log Space)'})
		subplot.add_hist(res['CoeffVar'].values()).update({'xlab':'CV','ylab': 'occurences','title': 'Coefficient of Variation Across Genes (Log Space)'})
		subplot.add_hist(res['#topVals'].values()).update({'xlab':'TopVals','ylab': 'occurences','title': 'Number of maximal values (top5)'})


		plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.90,wspace=0.1,hspace=0.40)
		subplot.save('sample_summary.png',{'title': 'Sample Summary Histograms'}) 
		
		rage_outputs.column_stats(self.args).write(res,self.input.samples,{'suffix': 'sample_stats.out','width': 20})

		self.progress.finish_subtopic() 



	def make_pca_and_tsne_plots(self):

		seaborn.set(rc={'axes.facecolor':'black', 'figure.facecolor':'cornflowerblue'})
		my_sizes = scale_vals([len(s.cnts.keys()) for s in self.input.samples],20,55)
		self.progress.start_subtopic('Calculating PCA/TSNE','',0)
		data_matrix = self.input.data_matrix('log')
		dr = rage_DR.DR(self.args,self.progress).set_matrix(data_matrix)
		dr.run_pca().run_kca().run_tsne().run_ica() 
		subplot = rage_subplots.subplot(2,2,self.args)
		subplot.add_legend(self.color_key.keys(),self.color_key.values())
		subplot.add_pca_data(dr.pca_pts,{'vars': dr.pca_vars,'title': 'PCA','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True}) 
		subplot.add_pca_data(dr.kca_pts,{'type': 'kca', 'title': 'KCA','colors':self.color_labels,'zoom': True,'sizes': my_sizes}).update({'clear_axes': True}) 
		subplot.add_pca_data(dr.ica_pts,{'type': 'ica', 'title': 'ICA','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True}) 
		subplot.add_pca_data(dr.tsne_pts,{'type': 'tsne','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True}) 
		#subplot.add_legend(self.color_key.keys(),self.color_key.values())
		subplot.save(self.args.prefix+'_dimred.png',{}) 
		self.progress.finish_subtopic() 





	def summarize_sample_pts(self,pt_label='val'):
		seaborn.set(rc={'axes.facecolor':'lightpink', 'figure.facecolor':'lightgray'})
		self.progress.start_subtopic('Plotting All Pts','',self.sLen)
		subplot = rage_subplots.subplot(1,1,self.args)  
		for si,s in enumerate(self.input.samples):
			self.progress.mark_subtopic() 
			if s.name[0] in ['U','H']: continue 
			#ordered_vals = [s.cnts[fi] for fi in range(len(self.input.features))]
			ordered_vals = [s.cnts[fi] for fi in self.feature_order] #range(len(self.input.features))]
			ordered_logs = [log(x+1.0) for x in ordered_vals]
			scaled_logs = scale_vals(ordered_logs) 
			s_color,s_mark = self.color_labels[si], self.mark_labels[si] 
			XY = [(x,scaled_logs[x]) for x in range(len(scaled_logs))]
			color_groups, group_colors = [[xy for xy in XY if xy[1] == 0]], [0]
			for (a,b,c) in [(d/20.0,(d+1)/20.0,(d+d+1.0)/40.0) for d in range(0,20)]:
				color_groups.append([xy for xy in XY if xy[1] > a and xy[1] <= b])
				group_colors.append(c) 
			diff_colors = get_colors(group_colors, plt.cm.jet) 
			for g,grp in enumerate(color_groups):
				if len(grp) == 0: 		continue
				elif grp[0][1] == 0.0: 
					clr, sz, alp = 'k',20,0.5 
					subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'mark': s_mark,'color': 'k', 'size': 0.1, 'alpha': 0.3,'yjitter': True})  
				else:
					clr, sz, alp  = diff_colors[g], (grp[0][1] + 1.2) * 2  , (0.1+grp[0][1])*0.6
					subplot.scatter_pts([gr[0] for gr in grp],[gr[1] for gr in grp],{'mark': s_mark, 'color': clr, 'size': sz, 'alpha': alp,'yjitter': True})  
		edge = int(len(XY)*0.05)
		subplot.ax.set_xlim(0-edge,len(XY)+edge)
		subplot.ax.text(len(XY)/2.5,-0.1,'Ordered Genes',fontweight='bold',fontsize=15)
		subplot.ax.text(len(XY)/2.6,1.1,'2k Single Cell Expression',fontweight='bold',fontsize=15)
		subplot.save(self.args.prefix+'_allpts.png',{'axis_off': True}) 
		self.progress.finish_subtopic() 
		


		


	def summarize_sample_dists(self):
		seaborn.set(rc={'axes.facecolor':'lightpink', 'figure.facecolor':'lightgray'})
		self.progress.start_subtopic('Plotting Sample Densities','',self.sLen)
		kde = rage_KDE.samples(0.3) 
		subplot,f_num = rage_subplots.subplot(10,10,self.args), 1 
		LOG=True
		for s in self.input.samples:
			self.progress.mark_subtopic() 	

			if LOG: non_zeros = [log(x+1.0) for x in s.cnts.values()]
			else:	non_zeros = [x for x in s.cnts.values()]
			all_vals = [0 for x in range(self.input.features.len-len(non_zeros))] + non_zeros

			x1,y1 = kde.run(all_vals)
			x2,y2 = kde.run(non_zeros)
			subplot.add_lines(x1,y1,None,None,'black')
			subplot.add_lines(x2,y2,None,None,'cyan')
			subplot.change_limits({'x0': -0.5,'x1': 8, 'y0': -0.1,'y1':1.5}) 
			subplot.ax.text(1.4,0.91,s.name+' ( '+str(len(non_zeros))+' )',color='blue')
			subplot.ax.set_xticklabels([]) 
			subplot.ax.set_yticklabels([]) 
			
			if not subplot.update({'clear_axes': True}): 
				plt.subplots_adjust(left=0.04, bottom=0.01, right=0.96, top=0.95,wspace=0.03,hspace=0.04)
				subplot.save(self.args.prefix+'fig_dists'+str(f_num)+'.png',{'title': 'Dual Dists: '})
				f_num += 1
				subplot = rage_subplots.subplot(10,10,self.args)  
		plt.subplots_adjust(left=0.02, bottom=0.01, right=0.98, top=0.95,wspace=0.03,hspace=0.03)
		subplot.save(self.args.prefix+'fig_dists'+str(f_num)+'.png',{'title': 'Dual Dists: '})
		self.progress.finish_subtopic() 

		

	def get_pairwise_dists(self):

		my_sizes = scale_vals([len(s.cnts.keys()) for s in self.input.samples],20,60)


		data_matrix = self.input.data_matrix('log')

		
#		print data_matrix.shape 

	
		### WEIGHS AS LONG AS FEATURES ###
		w=[100,3,9,122,5,99999,999999]
		
#		print len(w) 
		eud_dists = sds.pdist(data_matrix)
#		man_dists = sds.pdist(data_matrix,metric='mahalanobis')
		corr_dists = sds.pdist(data_matrix,metric='correlation')
#		mink_dists = sds.pdist(data_matrix,metric='wminkowski',w=w)
		mink_dists = sds.pdist(data_matrix,metric='minkowski')
#		print sds.squareform(corr_dists)


                mds_eud = rage_DR.DR(self.args,self.progress).run_mds(sds.squareform(eud_dists))
                mds_corr = rage_DR.DR(self.args,self.progress).run_mds(sds.squareform(corr_dists))
                mds_mink = rage_DR.DR(self.args,self.progress).run_mds(sds.squareform(mink_dists))
		
		seaborn.set(rc={'axes.facecolor':'black', 'figure.facecolor':'lightcyan'})


		subplot = rage_subplots.subplot(1,3,self.args)
                subplot.add_pca_data(mds_eud.mds_pts,{'type': 'mds-eud','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True})
                subplot.add_pca_data(mds_corr.mds_pts,{'type': 'mds-corr','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True})
                subplot.add_pca_data(mds_mink.mds_pts,{'type': 'mds-mink','colors':self.color_labels,'sizes': my_sizes}).update({'clear_axes': True})
                subplot.add_legend(self.color_key.keys(),self.color_key.values())
                subplot.save(self.args.prefix+'_mdsred.png',{})
		n=0
		print '--- s2 euclid corr mink' 
		for i in range(len(self.input.samples)): 
			s1 = self.input.samples[i] 
			for j in range(i+1,len(self.input.samples)):
				s2 = self.input.samples[j] 
				print s1.name, s2.name, eud_dists[n], corr_dists[n], mink_dists[n] 
				n+=1
		

		sys.exit() 

	def summarize_sample_pairs(self):


		self.get_pairwise_dists() 
		sys.exit() 
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



































	def create_label_key(self):

		circles,squares,triangles,others = ["o",".","8",'|'], ['s','p','D','_'],['>',"<","v","^"],['*','H','x',',']
		mark_list  = ['.']+[item for sublist in [[x[i] for x in [circles,triangles,squares,others]] for i in range(4)] for item in sublist]
		color_list = ['red','blue','orange','green','purple','magenta','lime','cyan','brown','crimson','darkred','r','g','y']
		color_list = ['orange','purple','magenta','brown','crimson','darkred','r','g','y']
		if len(self.input.samples.attributes) != 0: 
			j,k = 0,0 
			self.color_key = {'UHR': 'gray', 'HBR': 'silver','ADULT_TP': 'cyan','ADULT_HP_DG': 'blue', 'ADULT_HP_CA': 'dodgerblue','FETAL_ES':'yellow','FETAL_CTX':'green'}
			self.color_key['FETAL_CR'] = 'lime'
			self.color_key['FETAL_GERM'] = 'red'
			self.color_labels, self.mark_labels = [], [] 
			self.mark_key = {} 
			for s in self.input.samples: 
				s_att = s.attributes.keys()
				g_val = s.attributes['GRP']  
				s_val = s.attributes['LOC']
			#	s_val = s.attributes['GRP']
				if g_val == 'HP': 
					if s_val in ['DG','DG_SUB','SUB']:   s_val = 'ADULT_HP_DG'
					elif s_val in ['CA1','CA3']: s_val = 'ADULT_HP_CA' 
					else:			     s_val = 'NA'
				elif g_val == 'TP': s_val = 'ADULT_TP'
				elif g_val == 'EB': 	
					if s_val == 'SP_IZ_SVZ': s_val = 'NA' 
					elif s_val in ['IZ','IZ_SP','SP','CP','CP_SP']: s_val = 'FETAL_CTX'
					elif s_val in ['SVZ','IZ_SVZ']: s_val = 'FETAL_GERM'
					elif s_val == 'MZ': s_val = 'FETAL_CR'
				elif g_val == 'OB':	    s_val = 'FETAL_OB'
				elif g_val == 'ES':	    s_val = 'FETAL_ES'

				if s_val == 'NA' or s_val == 'UHR': 
					self.color_labels.append('white') 
					self.mark_labels.append('o') 
					continue

				if s_val not in self.color_key: 
					self.color_key[s_val] = color_list[k] 
					k+=1
				if s_val not in self.mark_key: 
					self.mark_key[s_val] = mark_list[j]
					j+=1
				self.mark_labels.append(self.mark_key[s_val])
				self.color_labels.append(self.color_key[s_val]) 
		else:

			self.color_key = {'EB': 'lime', 'T': 'red', 'ES': 'cyan', 'O': 'grey','U': 'purple','H': 'orange'}
			self.color_labels = [self.color_key[s.name[0]] if s.name[0] != 'E' else self.color_key[s.name[0:2]] for s in self.input.samples] 
			self.mark_labels = ['o' for s in self.input.samples]






