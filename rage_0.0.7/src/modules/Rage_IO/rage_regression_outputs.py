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






class eval_output: 
	def __init__(self,options):
		self.options = options

	def write(self,M,F,ext=None):
		v_str =  "-".join([p.split('=')[0] for p in self.options.predictors+self.options.covariates])

		S = M.X.segregator


		if len(set([z[0] for z in M.X.zp])) > 1: v_str += '_zpPrecomp'
		else:					 v_str += '_zpConst'



		if ext != None: w = open(self.options.prefix+'_'+M.dist+'_'+v_str+'-'+'eval_dex_brief_sim'+str(ext)+'.txt','w')
		else: w = open(self.options.prefix+'_'+M.dist+'_'+v_str+'-'+'eval_dex_brief.txt','w')

		w.write('%-50s %8s %15s %15s %10s %10s %10s %10s %10s %10s %10s\n' % ('---','VALID','History','pred','beta','p_obs','p_rate','e_rate','FC','cs','pv'))
		for i,f in enumerate(F):

			Y,P,V,H = M.Y[i], M.out['params'][i],M.out['VALID'][i], M.out['history'][i]
			p = sorted([pr for pr in P if pr[-1]])[0]
	
			pObs, pRate,eRate,FC,FS =  S.score(p[3],Y) 
			w.write('%-50s %8s %15s %15s %10f %10d %10.3f %10.3f %10.3f %10.2e %10.2e \n' % (f,V,H,p[2],p[1],pObs,pRate,eRate,FC,FS,p[0]))




class regression_output:

	def __init__(self,options,M):

		self.M = M  
		self.options = options 
		self.variables, self.filename =  [p.split('=')[0] for p in self.options.predictors+self.options.covariates],  ",".join([c.name for c in self.options.counts]) 
		self.gene_summary = open(self.options.prefix+'_featureSummary_'+self.options.dist+'_variables_'+"-".join(self.variables)+'.out' , 'w') 
		self.gene_summary.write('%-50s %8s %8s %8s %8s %8s %8s %8s %8s ' % ('---','LEN','OBS','RATE','AVG', 'AVGr','std','stdmr','stdtr'))
		self.gene_summary.write('%6s %7s %7s %30s %10s\n' % ('BIC','RSQ','RSA','TOP-P','TOP-PV'))
		
		self.gene_preds = open(self.options.prefix+'_featurePreds_'+self.options.dist+'_variables_'+"-".join(self.variables)+'.out' , 'w') 
		self.gene_preds.write('%-40s %7s %7s %7s %7s '% ('---','LEN','RATE','AVG', 'AVGr'))
		self.gene_preds.write('%25s %6s %5s %5s %5s %7s %7s %7s %7s %7s\n' % ("PRED","TYPE","pLen","pRATE","pAVG","pAVGr","BW","PV-RAW","PV-PR","PV"))


	def write_feature_summary(self,f,Yi,i):


		self.f_name = f.name 
		self.yLen,self.yObs,self.yRate,self.yR  = len(Yi),len([y for y in Yi if y>0]), len([y for y in Yi if y>0])/float(len(Yi)),self.M.out['covariate_resids'][i]
		
		cvY,cvR,cvTR = coVar(Yi), coVar(self.yR), coVar(self.M.out['resids'][i])

		cvT = np.std(Yi) 
		cvR =    np.std(self.yR)#/np.mean(Yi)  
		cvTR = np.std(self.M.out['resids'][i])


		self.yMean, self.rMean = np.mean(Yi), np.mean(self.yR) 


		pv,bw,topN,PRED = sorted([x for x in self.M.out['params'][i] if x[2] != 'intercept'])[0]  


		self.gene_summary.write('%-50s %8s %8s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f ' % (f.name,self.yLen,self.yObs,self.yRate,self.yMean,self.rMean,cvY,cvR,cvTR))
		self.gene_summary.write('%6.2f %7.3f %7.3f %30s %10.2e\n' % (self.M.out['bic'][i], self.M.out['rsq'][i], self.M.out['rsa'][i], topN,pv))




	def add_binary_predictor(self,n,raw_pv,pv,bw,perm,yN,rN):

		nLen, nObs,yXN,rXN = len(yN), float(len([yn for yn in yN if yn > 0])) / len(yN) , np.mean(yN), np.mean(rN)

		self.gene_preds.write('%-40s %7s %7.2f %7.2f %7.2f ' % (self.f_name,self.yLen,self.yRate,self.yMean,self.rMean))
		self.gene_preds.write('%25s %6s %5d %5.1f %5.1f %7.2f %7.2f %2.2e %2.2e %2.2e \n' % (n,"BIN",nLen,nObs,yXN,rXN,bw,raw_pv,perm,pv))


	def add_cont_predictor(self,n,raw_pv,pv,bw,perm):


		self.gene_preds.write('%-40s %7s %7.2f %7.2f %7.2f ' % (self.f_name,self.yLen,self.yRate,self.yMean,self.rMean))
		self.gene_preds.write('%25s %6s %5d %5.1f %5.1fs %7.2f %7.2f %2.2e %2.2e %2.2e \n' % (n,"CON",self.yLen,self.yRate,self.yMean,self.rMean,bw,raw_pv,perm,pv))







class predictor_output:

	def __init__(self,options,p,M,F,Y):

		self.M,self.P = M , p 
		self.F,self.Y = F,Y
		self.options = options 
		self.variables, self.filename =  [p.split('=')[0] for p in self.options.predictors+self.options.covariates],  ",".join([c.name for c in self.options.counts]) 


		self.pred_summary = open(self.options.prefix+'_predictorComp_'+self.P+'_covar_'+'-'.join(self.options.covariates)+'.out' , 'w') 
		self.pred_summary.write('%-40s %7s %7s %7s %6s ' % ('---','LEN',"RATE","AVG","SIMS"))
		self.pred_summary.write('%20s %7s %7s %7s %7s %7s %7s | PARAMS \n' % ('PRED','RSQ','ARS','SR','SA','S>R','S>A'))

	def add_sim_keys(self,gt_key,best_key,sims): 


			for i,(f,Yi) in enumerate(zip(self.F,self.Y)):
				self.yLen,self.yObs,self.yRate  = len(Yi),len([y for y in Yi if y>0]), len([y for y in Yi if y>0])/float(len(Yi))
				Mrsa,Mrsq = self.M.out['rsa'][i], self.M.out['rsq'][i]


				self.pred_summary.write('%-40s %7s %7.2f %7.2f %6d ' % (f.name,self.yLen,self.yRate,np.mean(Yi),sims))
				self.pred_summary.write('%20s %7.2f %7.2f %7.2f %7.2f %7d %7d ' % (self.P,Mrsq,Mrsa,best_key[i]['rsq'],best_key[i]['rsa'],gt_key[i]['rsq'],gt_key[i]['rsa']))

				for p in self.M.pv_dict.keys():
					if p != 'intercept': 
						self.pred_summary.write('| %5s %2.1e %3d ' % (p,self.M.pv_dict[p][i],gt_key[i][p]))
				self.pred_summary.write('\n') 



























class dex_result:
	def __init__(self,options,samplestr): 
		self.samplestr = samplestr
		self.options = options 
		self.cutoff = 0.1
		self.variables, self.filename =  [p.split('=')[0] for p in self.options.predictors+self.options.covariates],  ",".join([c.name for c in self.options.counts]) 
		
		 
		self.cnt_out  = open(self.options.prefix+'_dexcnts_variables_'+"-".join(self.variables)+'.out' , 'w') 
		self.res_out  = open(self.options.prefix+'_dexscores_variables_'+"-".join(self.variables)+'.out', 'w') 

		self.res_out.write('%-50s %5s %5s %5s %7s %7s %12s %9s' % ('---', 'len', 'obs', 'rate', 'rS','bic', 'top_param','pv'))
		self.res_out.write('%16s %7s %7s %8s %8s %10s %10s %10s %10s %10s %10s\n' % ('param', 'g_rank', 'p_rank', 'p_len','p_obs','fold', 'raw_fc', 'model_fc', 'model_bw','raw_pv','model_pv'))

	def add_feature(self,feature,Yi,y_new,rsquared = None, bic = None, top_variable=None, top_pv=None): 

		self.fLen,self.fObs,self.fRate, self.Yi, self.Yresid = len(Yi),len([y for y in Yi if y>0]), len([y for y in Yi if y>0])/float(len(Yi)), Yi, y_new
		self.feature = feature  
		self.rsquared, self.bic = rsquared, bic 
		self.top_variable, self.top_pv = top_variable, top_pv

	
	def write_continuous_result(self,variable,(total_rank,p_rank),(bw,raw_pv,model_pv),(raw_cnts,res_cnts),attributes):

	
		raw_means, res_means = [x[1] for x in sorted([(rk,np.mean(rv)) for rk,rv in raw_cnts.items()])], [x[1] for x in sorted([(rk,np.mean(rv)) for rk,rv in res_cnts.items()])]



		if raw_means[0] > raw_means[1]: raw_fold = raw_means[0]/(raw_means[1]+0.0001) 
		else: 				raw_fold = -1 * (raw_means[1]/(raw_means[0]+0.001))
		if res_means[0] > res_means[1]: res_fold = res_means[0]/(res_means[1]+0.0001) 
		else: 				res_fold = -1 * (res_means[1]/(res_means[0]+0.001))

		self.res_out.write('%-50s %5d %5d %5.2f %7.3f %7.2f %12s %9.1e '% (self.feature, self.fLen, self.fObs, self.fRate,self.rsquared,self.bic,self.top_variable.split('~')[0],self.top_pv))
		self.res_out.write('%16s %7d %7d %8d %8.2f %10s %10.2f %10.2f %10.3f %10.2e %10.2e\n' % (variable, total_rank, p_rank,self.fLen, self.fRate,'HI/LO', raw_fold, res_fold, bw, raw_pv,model_pv)) 

		if model_pv < self.cutoff: 
			self.cnt_out.write('%-50s %16s %5.2e %s %s %s %s %s\n' % (self.feature,variable,model_pv,'CONT',self.samplestr,  ",".join([str(a) for a in attributes]),",".join([str(y) for y in self.Yi]),",".join([str(y) for y in self.Yresid])  ))

			


	def write_binary_result(self,parent,results,raw_cnts,res_cnts,(yI,yR,attributes)): 




		raw_sums = {rk: sum(rv) for rk,rv in raw_cnts.items()}
		res_sums = {rk: sum(rv) for rk,rv in res_cnts.items()}	
		MPV = 1
	
		for variable,g_rank,p_rank,bw,raw_pv,model_pv in results:
			raw_obs,raw_total = len([x for x in raw_cnts[variable] if x > 0]), float(len(raw_cnts[variable]))
			raw_mean, raw_others = raw_sums[variable]/raw_total, sum([raw_sums[v] for v in raw_sums.keys() if v != variable]) / float(self.fLen-raw_total)
			res_mean, res_others = res_sums[variable]/raw_total, sum([res_sums[v] for v in res_sums.keys() if v != variable]) / float(self.fLen-raw_total)

			if raw_mean > raw_others: raw_fold = raw_mean /(raw_others+0.0001) 
			else: 				raw_fold = -1 * (raw_others/(raw_mean+0.001))
			if res_mean > res_others: res_fold = res_mean/(res_others+0.0001) 
			else: 				res_fold = -1 * (res_others/(res_mean+0.001))


			self.res_out.write('%-50s %5d %5d %5.2f %7.3f %7.2f  ' % (self.feature, self.fLen, self.fObs, self.fRate,self.rsquared,self.bic)) 
			self.res_out.write('%12s %9.1e ' % (self.top_variable.split('~')[0],self.top_pv))
			self.res_out.write('%16s %7d %7d %8d %8.2f %10s ' % (variable, g_rank, p_rank,int(raw_total), raw_obs/raw_total,'FC'))
			self.res_out.write('%10.2f %10.2f %10.3f %10.2e %10.2e\n' % (raw_fold,res_fold, bw, raw_pv,model_pv)) 

			if model_pv < MPV: MPV = model_pv 

		if MPV < self.cutoff: 
			raw_vals = [",".join([str(x) for x in raw_cnts[k]]) for k in raw_cnts.keys()]
			res_vals = [",".join([str(x) for x in res_cnts[k]]) for k in raw_cnts.keys()]
			
			rawS = ",".join([str(x) for x in yI])
			resS = ",".join([str(x) for x in yR])
			attr = ",".join([str(x) for x in attributes])
			self.cnt_out.write('%-50s %16s %5.2e %s %s %s %s %s\n' % (self.feature, parent, MPV,'BIN',self.samplestr,attr,rawS,resS))
			 
			

class reg_simulate:
	def __init__(self,options):

		self.options = options 

	def write(self,real_pvs, sim_pvs, sims, predictor):
		
		out = open(self.options.prefix+'_'+'_'.join(self.options.predictors)+'_'+'reg_simulate.txt','w') 
		out.write('%-15s %20s %13s %13s %13s %13s\n' % ('---','p-value','obs','sims','sim-obs','empirical-pv'))
		
		#number of real pvalues
		obs = real_pvs

		#pvalues of interest
		pvalues = [0.05, 0.01, 0.001, 0.0001, 0.00001, 0.0000001]
		
		#find number of pvalues 
		pv05 = [sim_pvs[i][0] for i in range(len(sim_pvs)) if sim_pvs[i][0] >= obs[0] and sim_pvs[i][0] != 0]
		pv01 = [sim_pvs[i][1] for i in range(len(sim_pvs)) if sim_pvs[i][1] >= obs[1] and sim_pvs[i][1] != 0]
		pv001 = [sim_pvs[i][2] for i in range(len(sim_pvs)) if sim_pvs[i][2] >= obs[2] and sim_pvs[i][2] != 0]
		pv0001 = [sim_pvs[i][3] for i in range(len(sim_pvs)) if sim_pvs[i][3] >= obs[3] and sim_pvs[i][3] != 0]
		pv00001  = [sim_pvs[i][4] for i in range(len(sim_pvs))if sim_pvs[i][4] >= obs[4] and sim_pvs[i][4] != 0]
		pv0000001 =[sim_pvs[i][5] for i in range(len(sim_pvs)) if sim_pvs[i][5] >= obs[5] and sim_pvs[i][5] != 0] 

		#sim observations for each pvalue
		sim_obs = [len(pv05), len(pv01), len(pv001), len(pv0001), len(pv00001), len(pv0000001)]
		
		# print ' '
		# print 'obs', obs 
		# print ' '
		# print 'pv05',pv05
		# print 'pv01',pv01
		# print 'pv001',pv001
		# print 'pv0001',pv0001
		# print 'pv00001',pv00001
		# print 'pv0000001',pv0000001
		# print ' '
		# print sim_obs, 'list of observations for all sims for each pvalue level'

		#empical p calculation
		empirical_pv = [len(pv05)/float(sims), len(pv01)/float(sims), len(pv001)/float(sims), len(pv0001)/float(sims), len(pv00001)/float(sims), len(pv0000001)/float(sims)]
		# print empirical_pv, 'calculated p values from sims'
		for i in range(len(pvalues)):
			out.write('%-15s %20f %13d %13d %13d %13.5f\n' % (self.options.predictors[0],pvalues[i],obs[i],sims,sim_obs[i],empirical_pv[i]))

		return




class regression_result:
	def __init__(self,options):

		self.options = options 
		 

	def write(self,D,M,Mp=None,suffix='dex'):



		w = open(self.options.prefix+'_'+suffix+'_'+"_".join(self.options.predictors)+'_covar'+"-".join(self.options.covariates)+'.out','w')

		if len(M.reg) != len(D.features): 
			print 'bad'
			sys.exit() 


		Mpreds = [M.X.names[i] for i in M.X.predictor_idx]
		Ppreds = [Mp.X.names[i] for i in Mp.X.predictor_idx]
		

		
		parents = [p for p in M.X.parents if M.X.PREDICTOR[p] and p!='intercept']
		children = [M.X.names[i] for i in M.X.predictor_idx]

		self.seg,self.fracs,self.segMM,self.segLens,self.mKey,self.pKey = {} ,{}, {} ,{},{},{}
		for p in parents:
			self.mKey[p] = {M.X.names[i]: i for i in M.X.predictor_idx if M.X.names[i] in M.X.children[p]}
			self.pKey[p] = {Mp.X.names[i]: i for i in Mp.X.predictor_idx if Mp.X.names[i] in Mp.X.children[p]}
			
			self.seg[p],self.segLens[p] = D.samples.segregate(p)
			self.fracs[p] = {g: float(v)/sum(self.segLens[p].values()) for g,v in self.segLens[p].items()}
#			self.segMM[p] = min(seg_lens.values()),max(seg_lens.values())


		w.write('--- RS CV obs len parent maxG maxMeans maxObs maxChi maxP | params\n')
		for zp,zm,f in zip(Mp.reg,M.reg,D.features):
			
			fObs = len(f.cnts)
			cnts = [f.cnts[i] for i in range(len(D.samples))]
			
			LS = [f.name,round(zm.model.rsquared,3),round(coVar(cnts),3)]
			
			for p in parents:
				
				maxC,maxP = sorted([(np.mean([f.cnts[i] for i in grp]),k) for k,grp in self.seg[p].items()])[-1]
				maxObs= len([x for x in [f.cnts[i] for i in self.seg[p][maxP]] if x >0])/float(self.segLens[p][maxP])
				p_srt = sorted([a for b in [[(f.cnts[i],k) for i in grp] for k,grp in self.seg[p].items()] for a in b],reverse=True)

				if len([x for x in p_srt if x[0]>0]) <5: 
					continue 


				else:
					if p_srt[self.segLens[p][maxP]-1][0] == 0:	p_split = [ps for ps in p_srt if ps[0] > 0]
					else:						p_split = p_srt[0:self.segLens[p][maxP]]	
					maxL,maxN,maxF,maxFn = len([ps[1] for ps in p_split if ps[1] == maxP]),len([ps[1] for ps in p_split if ps[1] != maxP]), self.fracs[p][maxP], 1-self.fracs[p][maxP]
					chiVal =  chisquare([maxL,maxN], f_exp=[len(p_split)*maxF,len(p_split)*maxFn])[1]


				try: line_data = LS+[fObs,len(cnts),p,maxP,round(maxC,3),round(maxObs,3),'%2.2e' % chiVal,'%2.2e' %  zm.model.pvalues[self.mKey[p][maxP]]]
				except KeyError:  line_data = LS+[fObs,len(cnts),p,maxP,round(maxC,2),round(maxObs,3),'%1.1e' % chiVal,'%1.1e' %  np.mean([zm.model.pvalues[zz] for zz in self.mKey[p].values()])]


				for x in self.pKey[p].keys():
					line_data.extend(['|',x,(zm.model.params[self.mKey[p][x]]>0),'%2.2e' % zm.model.pvalues[self.mKey[p][x]],'%2.2e' % zp.model.pvalues[self.pKey[p][x]]])
				w.write(" ".join([str(xx) for xx in line_data])+'\n' )


