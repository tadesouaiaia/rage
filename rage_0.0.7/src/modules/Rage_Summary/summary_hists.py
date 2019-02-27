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


#from modules.Rage_IO import rage_inputs, rage_counts, rage_outputs, rage_progress
#rom modules.Rage_Plots import rage_subplots

#from modules.Rage_Transforms import rage_KDE

#from modules.Rage_Transforms import rage_DR
#from modules.Rage_Summaries import rage_summarize_features 
#from modules.Rage_Summaries import rage_summarize_dists 



def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
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




def summary_hists(self):

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

	self.res = res
	subplot.add_hist(self.res['#Reads'].values()).update({'xlab':'reads per sample','ylab': 'occurences','title': 'Depth'})	

	print len(self.res['#Reads'])
	sys.exit() 

	subplot.add_hist(res['#Observed_Genes'].values()).update({'xlab':'genes per sample','ylab': 'occurences','title': 'Library Complexity'})	
	subplot.add_hist(res['#Genes_Above_Mean'].values()).update({'xlab':'%','ylab': 'occurences','title': '% genes above mean'})
	subplot.add_hist(res['%Genes_Required_For_HalfDepth'].values()).update({'xlab':'%Obs Genes','ylab': 'occurences','title': '% Genes Required For 50% Read Depth (Log Space)'})
	subplot.add_hist(res['CoeffVar'].values()).update({'xlab':'CV','ylab': 'occurences','title': 'Coefficient of Variation Across Genes (Log Space)'})
	subplot.add_hist(res['#topVals'].values()).update({'xlab':'TopVals','ylab': 'occurences','title': 'Number of maximal values (top5)'})


	plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.90,wspace=0.1,hspace=0.40)
	subplot.save('sample_summary.png',{'title': 'Sample Summary Histograms'}) 
	
	rage_outputs.column_stats(self.args).write(res,self.input.samples,{'suffix': 'sample_stats.out','width': 20})

	self.progress.finish_subtopic() 






































































































































class Samples:
	
        def __init__(self,summary):

		self.progress = summary.progress


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






