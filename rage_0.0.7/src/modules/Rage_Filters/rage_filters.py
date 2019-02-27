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

from random import random
import statsmodels.api as sm
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import fdrcorrection as fdr
import random
from math import fabs
#from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
import pickle
from math import log
import math
import numpy as np 
import pylab 
from matplotlib.patches import Rectangle as Rect
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib._png import read_png
from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import Circle as Circ

from sklearn.decomposition import PCA

from sklearn.manifold import TSNE

from random import shuffle
				
from sklearn.cluster import KMeans	


import seaborn
from sklearn.cluster import KMeans

from sklearn.neighbors import KernelDensity

from sklearn.preprocessing import MinMaxScaler



class dr_pts: 
        def __init__(self,options):
		self.options = options 

	def write(self,pts,col_members,format_key): 

		if 'suffix' in format_key: out_name = self.options.prefix+'_'+format_key['suffix']+'.txt'
		else:			   out_name = self.options.prefix+'.out'

		w = open(out_name,'w') 
		format_str = "%-10s %s\n"
		w.write(format_str % ('---'," ".join([s.name for s in col_members])))


		for i in range(min(10,len(pts[0]))):
			w.write(format_str % ('comp'+str(i+1)," ".join([str(p[i]) for p in pts])))


class dr_matrix: 
        def __init__(self,options):
		self.options = options 

	def write(self,matrix,col_members,format_key): 
		if 'suffix' in format_key: out_name = self.options.prefix+'_'+format_key['suffix']+'.txt'
		else:			   out_name = self.options.prefix+'.out'
		w = open(out_name,'w') 
		format_str = "%-10s %s\n"
		w.write(format_str % ('---'," ".join([s.name for s in col_members])))
		k=0
		for i in range(matrix.shape[1]):
			w.write(format_str % ('comp'+str(i+1)," ".join([str(matrix.flat[j+k]) for j in range(len(col_members))])))
			k+=len(col_members) 
	
					
		
		

			
		


		


class count_file: 
        def __init__(self,options):
		self.options = options 
	def write_row_col_data(self,row_members, col_members, data,format_key): 

		if 'name' in format_key: 
			out_name = format_key['name']
			w=open(out_name,'w') 
		else:
			w = sys.stdout
		format_str = "%-10s %s\n"
		w.write(format_str % ('---'," ".join([s.name for s in col_members])))
		for f,d in zip(row_members,data):
			 
			w.write(format_str % (f.name," ".join([str(x) for x in d])))

class column_coefs:
        def __init__(self,options):
		self.options = options 



	def write(self,coefs,rows,format_key):
		
		if 'suffix' in format_key: out_name = self.options.prefix+'_'+format_key['suffix']+'.txt'
		else:			   out_name = self.options.prefix+'.out'
		w = open(out_name,'w') 
		

		r_key,r_top, format_str= dd(list), ['---'], '%-40s'
		for i,C in enumerate(coefs): 	
			if i > 5: break 
			format_str += '  %s  %s' 
			r_top.extend(['Rank_'+str(i+1),'Val_'+str(i+1)])
			for j,c in enumerate(C):
				vAbs,val,vidx = c 
				r_key[rows[vidx].name].extend([j+1,val]) 

		format_str += '\n'
		w.write(format_str % tuple(r_top))
		for r,v in r_key.items():
			w.write(format_str % tuple([r]+[str(x) for x in v]))












class column_stats:
        def __init__(self,options):

		self.options = options 
		self.prefix = 'foo' 

	def write(self,column_key,row_names,format_key):



		if 'suffix' in format_key: out_name = self.options.prefix+'_'+format_key['suffix']+'.txt'
		else:			   out_name = self.options.prefix+'.out'
		if 'width'  in format_key: format_s, format_e = '%-'+str(format_key['width'])+'s ','%'+str(format_key['width'])+'s'
		else:			   format_s,format_e = '%-20s ','%20s'
		w = open(out_name,'w') 
		column_names = column_key.keys() 
		string_format = format_s+' '.join([format_e for s in column_names])+'\n'
		w.write(string_format % tuple(['---']+column_names))
		for n in row_names: 
			row_data = [str(x) for x in [n.name]+[column_key[k][n] for k in column_names]]	
			w.write(string_format % tuple(row_data))
	
		w.close() 

		


	def calculate_sample_stats(self):
		sample_stats = dd(lambda: dd(bool))
		sample_stats['TOTAL_FEATURES'] = {S: len([c for c in C if c>0]) for S,C in zip(self.samples,self.sample_cnts)}
		sample_stats['TOTAL_READS'] = {S: sum(C) for S,C in zip(self.samples,self.sample_cnts)}
		sample_stats['LOG_TOTAL_READS'] = {S: log(x) for S,x in sample_stats['TOTAL_READS'].items()}
		sample_stats['TOTAL_LOG_READS'] = {S: sum([log(c+1.0) for c in C]) for S,C in zip(self.samples,self.sample_cnts)}
		sample_stats['READS_PER_GENE'] = {S: sum(C)/sample_stats['TOTAL_FEATURES'][S] for S,C in zip(self.samples,self.sample_cnts)}

		self.sample_stats = sample_stats
		return sample_stats



	

	def calculate_feature_stats(self):

		sample_features    = {S: [self.features[j] for j in range(len(C)) if C[j]>0] for S,C in zip(self.samples,self.sample_cnts)}
		sample_read_totals = {S: sum(C) for S,C in zip(self.samples,self.sample_cnts)}

		f_cnts  = {f: sum([c for c in C if c>0]) for f,C in zip(self.features,self.feature_cnts)}
		f_total_cnts = float(sum(f_cnts.values()))
		feature_cnt_rates = {f: c/f_total_cnts for f,c in f_cnts.items()}
		

		s_exp_cnts = dd(lambda: dd(bool)) 
		for s,F in sample_features.items():
			s_rates = [(feature_cnt_rates[f],f) for f in F]
			s_rate_sum = sum([x[0] for x in s_rates])
			s_rel_rates = [(sample_read_totals[s]*(x[0]/s_rate_sum),x[1]) for x in s_rates] 
			s_exp_cnts[s] = {x[1]: x[0] for x in s_rel_rates}

			for f,fr in feature_cnt_rates.items():
				if f not in F: s_exp_cnts[s][f] = sample_read_totals[s]*(fr/(fr+s_rate_sum))


		f_obs  = {f: len([f for c in C if c>0]) for f,C in zip(self.features,self.feature_cnts)}
		f_obs_total = float(sum(f_obs.values()))
		feature_obs_rates = {f: c/f_obs_total for f,c in f_obs.items()}

		s_true_cnts = dd(lambda: {}) 
		for f,C in zip(self.features,self.feature_cnts): 
			for i in range(len(C)):
				s_true_cnts[self.samples[i]][f] = C[i]

		print '---','total_genes','total_reads','feature','fr','f_obsR','f_expC_giveO','true','expected'
		for s,F in sample_features.items():
			for f,r in feature_obs_rates.items():
				#if f not in s_exp_cnts[s]: s_exp_cnts[s][f] = 0.0
				op = 1.0 - ((1.0-r)**len(F)) 
				print s,len(F),sample_read_totals[s],f,r,op, s_exp_cnts[s][f], s_true_cnts[s][f], s_exp_cnts[s][f]*op 


	def read_counts(self):

		ROWS,COLS =  [x.upper() for x in self.options.order.split(',')]
		if COLS in ['SAMPLES','SAMPLES']:
			for line in self.options.counts:
				line = line.split() 
				if line[0] == '---': 
					self.samples = line[1::] 
					self.sample_idxs = range(len(line[1::]))
				else:
					vals = [float(x) for x in line[1::]]
					self.features.append(line[0]) 
					self.feature_cnts.append(vals)
					for i in range(len(vals)):
						if vals[i] > 0: 
							self.sample_vals[self.samples[i]][line[0]] = vals[i] 
							self.feature_vals[line[0]][self.samples[i]] = vals[i]



 
					#self.log_vals.append([math.log(v+1.0) for v in self.vals[-1]])
			
		elif COLS in ['SAMPLES,SAMPLE,CELL,CELLS,ANIMAL,ANIMALS']:
			print 'cool order'
			sys.exit() 

		else:
			print 'bad order'
			sys.exit() 	

		self.sample_cnts = [[self.feature_cnts[i][j] for i in range(len(self.features))] for j in range(len(self.samples))]
		self.sample_sums = [sum(c) for c in self.sample_cnts]
		self.sample_obs = {s: len(self.sample_vals[s]) for s in self.sample_vals}


		mp = max(self.sample_sums) 
	
		

	
#		self.feature_fracs = [[log(1.0+((mp*self.feature_cnts[i][j])/self.sample_sums[j])) for j in range(len(self.sample_sums))] for i in range(len(self.feature_cnts))]

		for s in self.samples: self.sample_key['PRESENT'][s] = True 
		for f in self.features: self.feature_key['PRESENT'][f] = True
		 

	def read_keys(self):
		if self.options.sampleKey: 
			for line in self.options.sampleKey: 
				line = line.split() 
				if line[0] == '---': headers = line
				else:		
					for i in range(1,len(line)):	self.sample_key[headers[i]][line[0]] = line[i] 

		
			
		if self.options.featureKey: 
			for line in self.options.featureKey: 
				line = line.split() 
				if line[0] == '---': headers = line
				else:		
					for i in range(1,len(line)):	self.feature_key[headers[i]][line[0]] = line[i] 

			
		
			




	def prepare_labels(self):
		try: 
			self.sample_color = {k: self.color_key[self.sample_key[k]] for k in self.sample_key.keys()}
		except KeyError:
			self.color_key = {} 
			for i,k in enumerate(list(set(self.key.values()))):
				self.color_key[k] = self.colors[i] 	
			self.sample_color = {k: self.color_key[self.sample_key[k]] for k in self.sample_key.keys()}
		
		self.ncol = len(list(set(self.key.values())))
		self.labels,self.items = [],[] 
		for xI,yI in self.color_key.items():
			if xI == 'NA': continue  
			self.labels.append(xI) 
			self.items.append(Rect((0,0),1,1,fc=yI))







		
