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
from scipy.stats import chisquare


import seaborn
from sklearn.cluster import KMeans

from sklearn.neighbors import KernelDensity

from sklearn.preprocessing import MinMaxScaler












class Classifier_Output:
	def __init__(self,options):
		self.options = options 
		self.f_str = self.options.prefix +'-Classify_L_'+str(self.options.leave)+"_".join(self.options.id)+'_Covary_'+"_".join(self.options.covariates)

		self.res = open(self.f_str+'.result','w') 
		self.coefs = open(self.f_str+'.coefs','w') 
		self.grades = open(self.f_str+'.grades','w') 
#		self.grades = sys.stdout

		self.grades.write('%-30s %15s %10s %10s %10s %10s %10s %10s\n'% ('---','GRADE','AVG','LEN','AVG(PASS)','LEN(PASS)','AVG(FAIL','LEN(FAIL'))
		self.res.write('%-20s %15s %15s %15s %15s %15s\n'% ('---','ID_TRUE','ID_PRED','TRUE_SCORE', 'PRED_SCORE','MATCH'))
		self.coefs.write('%-40s %40s %20s | %20s\n'% ('---','TARGETS','AVG-FRACS/SCORES','RANKS'))



	def add_score(self,p_data,p_tuples):

		s_name,s_known,s_true = p_data
#		['EB1263', True, 'SURE_MEGA~SMALL'] [('SURE_MEGA~SMALL', 0.07293598232705056)]

		self.res.write('%-20s %15s %15s            %15s %15s\n'% (s_name,s_known,s_true,p_tuples[0][0],p_tuples[0][1]))



	def add_coefs(self,f_name, coef_key, target_names): 



		if len(target_names) > 2: 



			t_dict,t_ranks = {target_names[t]: coef_key[t] for t in range(len(target_names))}, [] 
			t_names, t_len = sorted(target_names)  , float(len(t_dict.values()[0]))

			rank_dict = dd(lambda: dd(int)) 
			for i in range(len(t_dict.values()[0])):
				t_srt = sorted([(t_dict[k][i],k) for k in t_dict.keys()])
				if t_srt[0][0] != t_srt[1][0]: 
					rank_dict[t_srt[0][1]][0]  += 1 
				if t_srt[-2][0] != t_srt[-1][0]: 
					rank_dict[t_srt[-1][1]][1] += 1
			for t in t_names:
				t_ranks.extend([str(int(100*round(rank_dict[t][0]/ t_len,2)))+' '+str(int(100*round(rank_dict[t][1] / t_len, 2)))])
			t_means = sorted([(target_names[t],np.mean(coef_key[t])) for t in range(len(target_names))])
			name_sort = " ".join([t[0] for t in t_means])
			avg_sort = " ".join([str(round(t[1],5)) for t in t_means])
			self.coefs.write('%-40s %40s %20s |   %20s\n'% (f_name,name_sort,avg_sort,"   ".join(t_ranks)))

		else:
			x1,x2,xt = 0,0,0 
			for v in coef_key[0]:
				if v < 0: x1+=-1*v
				else:     x2+=v 
			
			v1 = x1 / float(x1+x2+0.00000000001) 
			v2 = x2 / float(x1+x2+0.00000000001) 
			
			 
			self.coefs.write('%-40s %40s %20s %20s %20s %20s |   %20s\n'% (f_name," ".join(target_names),v1,v2,x1,x2,'NA'))
                      
		











	def add_coefs2(self,f_name, coef_key, target_names): 

		t_means = sorted([(target_names[t],np.mean(coef_key[t])) for t in range(len(target_names))])
		name_sort = " ".join([t[0] for t in t_means])
		avg_sort = " ".join([str(t[1]) for t in t_means])
		self.coefs.write('%-40s %40s %40s \n'% (f_name,name_sort,avg_sort))
                      
	def add_gene_grades(self,g,GK): 		
		

		for k,V in GK.items(): 
			self.grades.write('%-30s %15s %10f %10f %10f %10f %10f %10f\n'% (g,k,V[0],V[1],V[2],V[3],V[4],V[5]))

class Marker_Output:
	def __init__(self,options):
		self.options = options 
		self.f_str = self.options.prefix +'-Marker_'+"_".join(self.options.id)+'_Covary_'+"_".join(self.options.covariates)

		self.res = open(self.f_str+'.data','w') 
	#	self.coefs = open(self.f_str+'.coefs','w') 
		self.res.write('%-20s %15s %15s %15s %15s %15s %15s\n'% ('---','MARKER-OPT','obs','enrich', 'ofc','qfc','MARKER'))
	#	self.coefs.write('%-40s %40s %20s | %20s\n'% ('---','TARGETS','AVG-SCORES','RANKS'))



	def add_marker(self,f,K):

		if len(K) > 1:

			
			self.res.write('%-20s %15s %15f %15f %15f %15f %15s\n'% (f,K['opt'],K['obs'],K['enrich'],K['ofc'],K['qfc'],K['MARKER']))
		else:
			self.res.write('%-20s %15s %15s %15s %15s %15s %15s\n'% (f,'NA','NA','NA', 'NA','NA','NA'))

#		s_name,s_known,s_true = p_data
#		['EB1263', True, 'SURE_MEGA~SMALL'] [('SURE_MEGA~SMALL', 0.07293598232705056)]

#		self.res.write('%-20s %15s %15s            %15s %15s\n'% (s_name,s_known,s_true,p_tuples[0][0],p_tuples[0][1]))



class Classifier_Unknown_Output:
	def __init__(self,options):
		self.options = options 
		self.f_str = self.options.prefix +'_UNKNOWN_CLASSIFICATION_LEAVEOUT_'+str(options.leave)+"_".join(self.options.id)+'_Covary_'+"_".join(self.options.covariates)

		self.res = open(self.f_str+'.predictions','w') 
#		self.res = sys.stdout
		self.res.write('%-20s %20s %20s %20s\n'% ('---','ID_TRUE','ID_PRED','WEIGHT'))
		


	def add_pred(self,name,id_true,id_pred,weight):
		
		self.res.write('%-20s %20s %20s %20s\n'% (name,id_true,id_pred,round(weight,5)))
