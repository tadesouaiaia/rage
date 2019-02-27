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

from operator import mul

from scipy.signal import savgol_filter as svgf 
from math import exp
from math import factorial 


from Rage_IO import rage_outputs

from Rage_Transforms import rage_DR
from Rage_Plots import rage_scatterplots
from Rage_Plots import rage_subplots
from Rage_Comps import rage_comps


#from Rage_Plots import rage_subplots
#from Rage_Plots import  rage_dimplots as dplot 

#from Rage_Transforms import rage_KDE
#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 



def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))




class Summary:
        def __init__(self,rage):

		self.rage = rage 
		self.options = rage.args

	def run(self):
		R = self.rage 

		# R.data.filter_samples_by_attributes().normalize()


		if R.args.command == 'samples':
			R.progress.start_major('SampleSummary')
		
			R.data.samples.create_plot_labels(R.args) 

			if R.args.pca or R.args.tsne:  
				R.progress.start_minor('Performing Dimensional Reduction',len(R.data.samples))
				dim = rage_DR.DR(R.args,R.progress).set_fit_matrix(R.data.matrix('log'))
				pca = dim.pca() 
				if R.args.tsne:	dim_plot = rage_scatterplots.DimR(R.args,R.progress,1,2).add_dim_run(pca,R.data.samples).add_dim_run(dim.tsne(),R.data.samples).save()
				else:		dim_plot = rage_scatterplots.DimR(R.args,R.progress).add_dim_run(pca,R.data.samples).save()

				rage_outputs.column_coefs(R.args).write(pca['coefs'],R.data.features,{'suffix': 'PCAcoeffs.features.out','width': 15})
				rage_outputs.dr_pts(R.args).write(pca['pts'],R.data.samples,{'suffix': 'pca.pts.out'})


			R.progress.start_minor('Calculating Summary Stats',len(R.data.samples))
			sample_stats = summary_hists(R.data.samples,R.data.features,R.args,R.progress)  
			rage_outputs.column_stats(R.args).write(sample_stats,R.data.samples,{'suffix': 'samplestats','width': 15})

			sample_trends = summary_trends(R.data.samples,R.data.features,R.args,R.progress)


		elif R.args.command == 'features': 
			R.progress.start_major('FeatureSummary')

			R.data.features.create_plot_labels(R.args) 
			if R.args.pca or R.args.tsne:  
				R.progress.start_minor('Performing Dimensional Reduction',len(R.data.features))
				dim = rage_DR.DR(R.args,R.progress).set_fit_matrix(R.data.matrix('log',TRANSPOSE=True))
				pca = dim.pca() 
				if R.args.tsne:	dim_plot = rage_scatterplots.DimR(R.args,R.progress,1,2).add_dim_run(pca,R.data.features).add_dim_run(dim.tsne(),R.data.features).save()
				else:		dim_plot = rage_scatterplots.DimR(R.args,R.progress).add_dim_run(pca,R.data.features).save()

				rage_outputs.column_coefs(R.args).write(pca['coefs'],R.data.features,{'suffix': 'PCAcoeffs.features.out','width': 15})
				rage_outputs.dr_pts(R.args).write(pca['pts'],R.data.samples,{'suffix': 'pca.pts.out'})


			R.progress.start_minor('Calculating Summary Stats',len(R.data.samples))

			feature_stats = summary_hists(R.data.features,R.data.samples,R.args,R.progress)  
			rage_outputs.column_stats(R.args).write(feature_stats,R.data.features,{'suffix': 'featurestats.out','width': 15})

			feature_trends = summary_trends(R.data.features,R.data.samples,R.args,R.progress)



		elif R.args.command == 'ratios': 

			feature_comps = rage_comps.features(self.rage).get_f_ratios()

                	HOUSEKEEPING, r_key = feature_comps.HOUSEKEEPING, feature_comps.r_key

                	feature_comps.predict_known_ratio_values()















def summary_trends(X,Y,options,progress,X_NAME='SAMPLES'):

	seaborn.set(rc={'axes.facecolor':'lightcyan', 'figure.facecolor':'whitesmoke'})
	res, subplot = dd(lambda: {}), rage_subplots.subplot(3,3,options,{'titlepos': [0.0,1.05]}) 
	qts, maxV, obsR, means, totals, log_totals, cnt_means = [],[], [] , [] , [] , [] , [] 
	trends = {} 
	stds = [] 
	cvs       = [] 
	vsx = [] 
	for x in X:
		zC = len(Y)-len(x.cnts) 
		x_all = x.cnts.values() + [0 for s in range(len(Y)-len(x.cnts))]
		
		x_mean = np.mean(x.cnts.values()) 	

		#x_log = [log(p+1.0) for p in x.cnts.values()] + [0 for s in range(len(Y)-len(x.cnts))]

		#qts.append(log(1.0+np.percentile(x_all,95)))
		#maxV.append(log(1.0+max(x_all)))
		#means.append(log(1.0+np.mean(x_all)))
			
		

		qts.append(np.percentile(x_all,95))
		maxV.append(log(max(x_all),2))
		means.append(log(np.mean(x_all),2))
		obsR.append(len(x.cnts)/float(len(Y)))
		
		totals.append(sum(x_all))
		log_totals.append(log(sum(x_all)+1.0,2))
		cnt_means.append(log(x_mean,2)) 

		stds.append(np.std(x_all)) 
		vsx.append(log(np.var(x_all),2)) 
		cvs.append(coVar(x_all)) 


	trends[('observations','log_total')] = stats.pearsonr(log_totals,obsR)
	trends[('observations','cnt_mean')] = stats.pearsonr(cnt_means,obsR)
	trends[('observations','max')] = stats.pearsonr(maxV,obsR)
	subplot.add_scatter_trend(obsR,log_totals,R=trends[('observations','log_total')][0]).update({'title': 'Observations vs Total','xlab': 'observation rate','ylab': 'total (logs)'})
	subplot.add_scatter_trend(obsR,cnt_means,R=trends[('observations','cnt_mean')][0]).update({'title': 'Observations vs Cnt Mean','xlab': 'observation rate','ylab': 'Cnt Mean (logs)'})
	subplot.add_scatter_trend(obsR,maxV,R=trends[('observations','max')][0]).update({'title': 'Observations vs Max','xlab': 'observation rate','ylab': 'Max (logs)'})



	trends[('mean','var')] = stats.pearsonr(means,vsx)
	subplot.add_scatter_trend(means,vsx,R=trends[('mean','var')][0]).update({'title': 'Mean vs Variance','xlab': 'Mean (log)','ylab': 'Variance'})

	trends[('mean','cv')] = stats.pearsonr(means,cvs)
	subplot.add_scatter_trend(means,cvs,R=trends[('mean','cv')][0]).update({'title': 'Mean vs CV','xlab': 'Mean (log)','ylab': 'CV'})







	trends[('max','log_total')] = stats.pearsonr(log_totals,maxV)

	subplot.add_scatter_trend(qts,totals).update({'title': 'Upper Quartile vs Total','xlab': 'Upper Quartile','ylab': 'total (logs)'})
#	subplot.add_scatter_trend(qts,maxV).update({'title': 'Upper Quartile vs Max','xlab': 'Upper Quartile','ylab': 'Max (logs)'})
#	subplot.add_scatter_trend(means,qts).update({'title': 'Mean vs Upper Quartile','xlab': 'Mean (logs)','ylab': 'Upper Quartiles (logs)'})
	subplot.add_scatter_trend(means,maxV).update({'title': 'Mean vs Max','xlab': 'Mean','ylab': 'Max (logs)'})
	
	plt.subplots_adjust(left=0.05, bottom=0.04, right=0.95, top=0.90,wspace=0.25,hspace=0.5)
	if X.label == 'samples': 	subplot.save(options.prefix+'_sample_trends.png',{'title': 'Sample Trends'}) 
	elif X.label == 'features':	subplot.save(options.prefix+'_feature_trends.png',{'title': 'Feature Trends'}) 

	return 


def summary_hists(X,Y,options,progress,X_NAME='SAMPLES'):

	#seaborn.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'lightgray'})
	#seaborn.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'w'})
	seaborn.set(rc={'axes.facecolor':'lightcyan', 'figure.facecolor':'whitesmoke'})
	res, p_res, subplot = dd(lambda: {}), dd(lambda: {}), rage_subplots.subplot(3,2,options) 
	cMax = float(sum([y.len for y in Y]))
	for x in X: 
		progress.mark() 
		xMissed = [0 for i in range(len(Y) - len(x.cnts))]
		x_raw = [c for c in x.cnts] + xMissed
		res['#CompIndex'][x] = sum([Y[y].len for y in x.cnts.keys()]) / cMax
		ordered_logs = sorted([log(1.0+c) for c in x.cnts.values()],reverse=True)	

		try: res['#Obs_gtAvg'][x] = len([l for l in ordered_logs if l > np.mean(ordered_logs)])/float(len(ordered_logs))
		except ZeroDivisionError:  continue 


		try: res['#Obs_gtAvg'][x] = len([l for l in ordered_logs if l > np.mean(ordered_logs)])/float(len(ordered_logs))
		except ZeroDivisionError:   res['#Obs_gtAvg'] = 0 


		

		halfE,iX,k =  sum(ordered_logs) * 0.5, 0, -1 		
		p_res['log_total'][x]          =  x.cnt_total
		res['total'][x]          =  x.cnt_total
		res['observations'][x]            = len(x.cnts)
		res['Qrt-75'][x]          = np.percentile(x_raw,75) 
		res['Perc-90'][x]          = round(np.percentile(x_raw,90),4)
		res['Perc-95'][x]          = np.percentile(x_raw,95) 
		res['Perc-99'][x]          = np.percentile(x_raw,99) 
	

		while iX < halfE:
			k+=1;	iX+=ordered_logs[k] 
		res['%Obs_HDepth'][x] = k / float(len(x.cnts)+len(xMissed)) 
		res['CoeffVar'][x] = coVar(ordered_logs) 

	subplot.add_hist(res['total']).update({'xlab':'log(reads)','ylab': 'occurences','title': 'Total Depth'})	
	
	if X.label == 'samples': subplot.add_hist(res['observations']).update({'xlab':'observations','ylab': 'occurences','title': 'Library Diversity (Genes)'})	
	else: subplot.add_hist(res['observations']).update({'xlab':'observations','ylab': 'occurences','title': 'Library Diversity (Samples)'})	
#	subplot.add_hist(res['#Obs_AboveMean']).update({'xlab':'%','ylab': 'occurences','title': 'Percentage of counts above mean'})
	subplot.add_hist(res['Qrt-75']).update({'xlab':'cnts','ylab': 'occurences','title': 'Upper Quartile'})



	subplot.add_hist(res['CoeffVar']).update({'xlab':'CV','ylab': 'occurences','title': 'Coefficient of Variation (Log Space)'})
	subplot.add_hist(res['%Obs_HDepth']).update({'xlab':'%Obs','ylab': 'occurences','title': '% Obs For 50% Read Depth (Log Space)'})
	subplot.add_hist(res['#CompIndex']).update({'xlab':'%Comparisons','ylab': 'occurences','title': 'Comparison Index'})
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


def make_dr_plots2(R,choice='samples'):

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
#	tsne_run = dr.run_tsne() 

#	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,tsne_run,{'title':'TSNE','out': out_name+'tsne.pdf'})  
	kca_gamma = 20
#	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=0.001) 
#	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-001-rbf.pdf','zoom': True})  

#	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=0.0001) 
#	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-0001-rbf.pdf','zoom': True})  

#	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=0.1) 
#	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-1-rbf.pdf','zoom': True})  

	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=0.01) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-01-rbf.pdf','zoom': True})  
	
	kca_run2 = dr.run_kca(r_matrix,kernel='rbf',gamma=0.005) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-005-rbf.pdf','zoom': True})  

	
	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=0.05) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-05-rbf.pdf','zoom': True})  

	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=0.1) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-p1-rbf.pdf','zoom': True})  

	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=1) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-i1-rbf.pdf','zoom': True})  

	kca_run = dr.run_kca(r_matrix,kernel='rbf',gamma=10) 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_run,{'title':'KCA-rbf','out': out_name+'kca-i10-rbf.pdf','zoom': True})  

	kca_lin = dr.run_kca(r_matrix,kernel='linear') 
	dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_lin,{'title':'KCA-linear','out': out_name+'kca-line.pdf','zoom': True})  

	try:
		kca_poly = dr.run_kca(r_matrix,kernel='poly') 
		dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_poly,{'title':'KCA-poly','out': out_name+'kca-poly.pdf','zoom': True})  
	except np.linalg.linalg.LinAlgError:
		kca_poly=None

	try:
		kca_sig = dr.run_kca(r_matrix,kernel='sigmoid') 
		dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_sig,{'title':'KCA-sig','out': out_name+'kca-sig.pdf','zoom': True})  
	except np.linalg.linalg.LinAlgError:
		kca_poly=None
	try:	
		kca_cosine = dr.run_kca(r_matrix,kernel='cosine') 
		dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,kca_cosine,{'title':'KCA-cosine','out': out_name+'kca-cos.pdf','zoom': True})  
	except np.linalg.linalg.LinAlgError:
		kca_poly=None



	return pca_run, kca_run, kca_lin









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





















































