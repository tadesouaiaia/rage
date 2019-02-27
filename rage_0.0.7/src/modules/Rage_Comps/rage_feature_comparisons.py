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






	def get_totals(self,RSTAT='MEAN'): 

		
		if self.args.fhtotals == None: 
			if self.args.solver == None: 
				self.progress.start_minor('Selecting DropOut Aware Feature Parameters',len(self.D.features)*len(self.D.features)/2.0,True)
				self.calculate_solver() 
				self.progress.start_minor('Loading DropOut Aware Feature Parameters From File',2*len(self.D.features)*(len(self.D.features)/2.0),True)
			else:	
				self.progress.notate('File Found: Feature Relationships\n')
				self.progress.start_minor('Reading DropOut Aware Feature Parameters From File',len(self.D.features)*(len(self.D.features)/2.0),True)
		     		self.read_solver(self.args.solver) 

   			self.progress.start_minor('Calculating Inferred Totals',len(self.D.features),True)
			self.impute_sample_totals() 		
		else:
			
			self.progress.notate('File Found: FullHouse Totals\n')

			self.progress.start_minor('Reading FullHouse Inferred Totals From File',len(self.D.features)*(len(self.D.features)/2.0),True)
			self.read_imputed_totals(self.args.fhtotals) 

		return self		



	def read_imputed_totals(self,total_file):

		wAmp = open(self.args.prefix+'.ratioAmplification.rates','w') 
		wAmp.write('%-20s %20s %20s\n' % ('---','Inferred-Amp','Imputed-Amp'))
		m_key = {}  
		a_key = {} 
		self.amplification_key = {} 
		self.fraction_key = {} 
		for line in total_file:
			line = line.split() 
			sample = line[0] 
			trueVal, inferAll, inferMissing, imputeAll, imputeMissing = [float(x) for x in line[2::]] 
			inferObs = inferAll - (inferMissing-trueVal)
			imputeObs = imputeAll - (imputeMissing-trueVal)
			
			m_key[sample] = [inferMissing,imputeMissing] 
			a_key[sample] = [inferObs/trueVal, imputeObs/trueVal] 	
			self.fraction_key[sample]      = [trueVal/inferMissing, trueVal/imputeMissing] 		

			wAmp.write('%-20s %20s %20s\n' % (sample,a_key[sample][0],a_key[sample][1]))


		inferAmpMin  =   min([x[0] for x in a_key.values()]) 
		imputeAmpMin =   min([x[0] for x in a_key.values()]) 

		self.amplification_key = {a: [b[0]/inferAmpMin,b[1]/imputeAmpMin] for a,b in a_key.items()}




		self.minInfer  = min([x[0] for x in m_key.values()])
		self.minImpute = min([x[1] for x in m_key.values()])  
	
		

 			




	def read_solver(self,solver_file,S_CUTOFF=0.3,MATCH_CUT=0.05): 
		f_idx = {f.name: f.idx for f in self.D.features}
		r_header = solver_file.readline() 
		self.s_key, self.feature_stats, self.solver_stats = {}, {}, {}  
		self.solver_key, self.ratio_key = dd(lambda: {}), dd(lambda: {}) 

		for line in solver_file:

			self.progress.mark()
			if line[0] == '-': continue 
			line = [x.split() for x in line.split('\t')]	
			f1,f2 = line[0] 
			try: 			idx1,idx2= f_idx[f1],f_idx[f2] 
			except KeyError:	continue 

			obs1,perc1,match1 = line[1] 	
			obs2,perc2,match2 = line[2] 	
			Sr,Sv = [float(x) for x in line[3]]

			if match1 < MATCH_CUT or match2 < MATCH_CUT: continue 
			if Sv > 0.0001: continue 
		
			avg1,med1,sol1 = [float(x) for x in line[4]]
			avg2,med2,sol2 = [float(x) for x in line[5]]
			
			self.solver_key[idx1][idx2] = sol2 
			self.solver_key[idx2][idx1] = sol1
			self.ratio_key[idx1][idx2]  = avg2
			self.ratio_key[idx2][idx1]  = avg1

			

	def calculate_solver(self,RSTAT='MEAN'): 

#		self.s_key = dd(list) 
#		self.progress.start_minor('Selecting DropOut Aware Feature Parameters',len(self.D.features),True)
		wLut = open(self.args.prefix+'.ratioSolutions.out','w') 
		wLut.write('%-40s %40s\t%5s %7s %7s\t%5s %7s %7s\t%5s %8s\t' % ('---','---','obs1','perc1','match2','obs2','perc2','match2','Sc','Spv'))
		wLut.write('%11s %11s %11s\t%11s %11s %11s \n' % ('avg1','med1','lsqsol1','avg2','med2','lsqsol2'))
		SLF = float(len(self.D.samples)) 

		for fi in range(len(self.D.features)): 
			self.progress.mark() 
			f_cnts, f_key, m_key, x_cnts, y_cnts = [], dd(list), dd(list), dd(list), dd(list)  
			fL = float(len(self.D.features[fi].cnts))
			fLR = fL/SLF

			iName,iLen,iRT = self.D.features[fi].name, len(self.D.features[fi].cnts),len(self.D.features[fi].cnts)/SLF

			for si,ci in self.D.features[fi].cnts.items():
				for fj,cj in [st for st in self.D.samples[si].cnts.items() if st[0] > fi]: 
					x_cnts[fj].append(ci) 
					y_cnts[fj].append(cj) 
					f_key[(fi,fj)].append(ci/float(cj))

			for fj in x_cnts.keys(): 
				jName,jLen,jRT = self.D.features[fj].name,len(self.D.features[fj].cnts),len(self.D.features[fj].cnts)/SLF
				xC,yC = x_cnts[fj],y_cnts[fj]
				xCA,yCA = np.array(xC),np.array(yC) 
				matA = len(xC)/float(iLen) 
				matB = len(xC)/float(jLen) 
				SR,SRv = stats.spearmanr(xC,yC) 
				SR2= SR*SR
				
				RM = np.mean([x/y for x,y in zip(xC,yC)]) 
				ZM = np.median([x/y for x,y in zip(xC,yC)]) 
#				print f.name, preds, Sp[0],Sp[1],Rp[0],Rp[1],'|',Se,R
				RY = np.mean([y/x for x,y in zip(xC,yC)]) 
				ZY = np.median([y/x for x,y in zip(xC,yC)]) 
				rawA = np.linalg.lstsq(xCA[:,np.newaxis],yCA,rcond=None)[0][0] 
				rawB = np.linalg.lstsq(yCA[:,np.newaxis],xCA,rcond=None)[0][0] 

				
				wLut.write('%-40s %40s\t%5d %7.4f %7.4f\t%5d %7.4f %7.4f\t' % (iName,jName,iLen,iRT,matA,jLen,jRT,matB))
				wLut.write('%5.3f %8.1e\t%11.4f %11.4f %11.5f\t%11.4f %11.4f %11.5f\n' % (SR2,SRv,RM,ZM,rawA,RY,ZY,rawB))

				#self.s_key[fj].append((fi,iLen,matA,SR2,SRv,RM,rawA))
				#self.s_key[fi].append((fj,jLen,matB,SR2,SRv,RY,rawB))

		wLut.close() 
		wLut = open(self.args.prefix+'.ratioSolutions.out')
		self.read_solver(wLut)









	def impute_sample_totals(self,BREAKLEN=500,MINLEN=100):
		wImpute, wInfer = open(self.args.prefix+'.imputed.cnts','w'), open(self.args.prefix+'.inferred.cnts','w') 
		wTotals, wPerf = open(self.args.prefix+'.fullhouse.totals','w'), open(self.args.prefix+'.feature.perfs','w') 

		wPerf.write('%-50s %6s %6s %6s %6s %6s | %10s %10s\n' % ('---','preds','impR','impRpv','infR','infRpv','impMSE','infMSE'))

		wImpute.write('%s %s\n' % ('---'," ".join([s.name for s in self.D.samples])))
		wInfer.write('%s %s\n' % ('---'," ".join([s.name for s in self.D.samples])))
		s_infer_missing, s_impute_missing, s_inferred, s_imputed = dd(list) , dd(list) , dd(list) , dd(list)

		for i,f in enumerate(self.D.features): 
			self.progress.mark() 
			f_impute, f_infer, preds = [], [], 0 
			comp_real, comp_impute, comp_infer = [], [] ,[] 

			for si,s in enumerate(self.D.samples):
				R_PRED,S_PRED, rP,sP, sVal = [], [] , 0, 0 ,  0 
				if i in s.cnts: sVal += s.cnts[i] 

				for fj in self.solver_key[i].keys(): 
					if fj not in s.cnts: continue 
					fSol =  self.solver_key[i][fj] 
					fAvg =  self.ratio_key[i][fj] 
					R_PRED.append(s.cnts[fj] / fAvg) 
					S_PRED.append(s.cnts[fj] * fSol) 

				if len(S_PRED) > 1: 
					preds+=1
					rP += np.mean(R_PRED)
					sP += np.mean(S_PRED)
				
				f_impute.append(sP) 
				f_infer.append(rP) 
				if sVal>0:      
						s_impute_missing[s].append(sVal)
						s_infer_missing[s].append(sVal)
				else: 		
						s_impute_missing[s].append(sP) 
						s_infer_missing[s].append(rP) 
				s_inferred[s].append(rP) 
				s_imputed[s].append(sP) 
				if preds > 2 and sVal > 0: 
					comp_real.append(sVal) 
					comp_impute.append(sP) 
					comp_infer.append(rP) 
					
			wImpute.write('%s %s\n' % (f.name," ".join([str(int(round(x,1))) for x in f_impute])))
			wInfer.write('%s %s\n' % (f.name," ".join([str(int(round(x,1))) for x in f_infer])))

			if len(comp_real) > 2: 
				Sp,Rp = stats.pearsonr(comp_real,comp_impute), stats.pearsonr(comp_real,comp_infer)
				Se = np.mean([(x-y)*(x-y) for x,y in zip(comp_real,comp_impute)])
				Re = np.mean([(x-y)*(x-y) for x,y in zip(comp_real,comp_infer)])
				wPerf.write('%-50s %6d %6.3f %6.3e %6.3f %6.3e | %10.4f %10.4f\n' % (f.name,preds, Sp[0], Sp[1], Rp[0], Rp[1], Se, Re))

		for s in self.D.samples:
			sInferAll, sInfer, sImputeAll, sImpute     = 	sum(s_inferred[s]), sum(s_infer_missing[s]), sum(s_imputed[s]), sum(s_impute_missing[s]) 
			wTotals.write('%s True/Infer/Infer-Missing/Impute/Impute-Missing %d %d %d %d %d\n' % (s.name,sum(s.cnts.values()),sInferAll,sInfer,sImputeAll,sImpute ))

		wTotals.close() 
		wTotals = open(self.args.prefix+'.fullhouse.totals')
		self.read_imputed_totals(wTotals)
		



