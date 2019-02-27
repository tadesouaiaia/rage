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



	def read_f_ratios(self,RSTAT='MED'): 

		self.progress.start_minor('Reading DropOut Aware Feature Parameters From File',len(self.D.features)*(len(self.D.features)/2.0),True)
		f_idx = {f.name: f.idx for f in self.D.features}
		r_header = self.args.ratios.readline() 

		if self.args.ratios.name.split('.')[-2].upper() == 'RATIOSOLVER': 

			for line in self.args.ratios:
				self.progress.mark()
				line = line.split() 
				iName,jName,_,iLen,iRate,matchA,_,SR2,SRv,_,RM,rawA = line 
				self.s_key[f_idx[jName]].append((f_idx[iName],int(iLen),float(matchA),float(SR2),float(SRv),float(RM),float(rawA)))



		else:
			for line in self.args.ratios:
				self.progress.mark()
				line = line.split() 
				f1,f2,fRate = line[0],line[1],float(line[-1])

				try: 
					self.HOUSEKEEPING[f_idx[f1]] = True
					self.HOUSEKEEPING[f_idx[f2]] = True
					self.r_key[(f_idx[f1],f_idx[f2])] = fRate
				except KeyError: continue
		


	def read_solution(self,solver_file,S_CUTOFF=0.3,MATCH_CUT=0.5): 
		self.progress.start_minor('Reading DropOut Aware Feature Parameters From File',len(self.D.features)*(len(self.D.features)/2.0),True)
		f_idx = {f.name: f.idx for f in self.D.features}
		#r_header = self.args.solver.readline() 

		r_header = solver_file.readline() 

		self.s_key = {} 
		self.solver_key = dd(lambda: {}) 
		self.ratio_key =  dd(lambda: {}) 
		self.feature_stats = {} 

		self.solver_stats = {} 

		#for line in self.args.solver:
		#for line in self.args.solver:
		for line in solver_file:


			self.progress.mark()
			if line[0] == '-': continue 
			#line = [x.split() for x in line.split('|')]	
			line = [x.split() for x in line.split('\t')]	

			f1,f2 = line[0] 

			try: 
				idx1,idx2= f_idx[f1],f_idx[f2] 
			except KeyError:
				continue 



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

			














	def get_f_ratios(self,RSTAT='MEAN'): 


		self.r_key = {}
		self.RSTAT = RSTAT 

		self.s_key = dd(list) 

		self.HOUSEKEEPING = dd(bool) 
		if self.args.ratios != None: 
			self.read_f_ratios()
		elif self.args.solver != None:
			self.read_solution(self.args.solver) 
 
		else:
   			self.progress.start_minor('Selecting DropOut Aware Feature Parameters',len(self.D.features),True)
#			wMean = open(self.args.prefix+'.meanRatios.out','w') 
#			wMean.write('%-40s %40s %5s %5s | %10s %25s %25s\n' % ('---','---','obs1','obs2','matches','medRatio','meanRatio'))
				
#			wSol = open(self.args.prefix+'.ratioSolver.out','w') 
#			wSol.write('%-40s %40s | %5s %10s %10s | %6s %8s | %10s %10s\n' % ('---','---','obs1','rate1','matcheRate','S_2','S_pv','meanRatio','lsqSol'))


			wLut = open(self.args.prefix+'.ratioSolutions.out','w') 
#			wLut.write('%-40s %40s | %5s %10s %10s | %6s %8s | %10s %10s\n' % ('---','---','obs1','rate1','matcheRate','S_2','S_pv','meanRatio','lsqSol'))

			#wLut.write('%-40s %40s | %5s %7s %7s | %5s %7s %7s | %5s %8s | ' % ('---','---','obs1','perc1','match2','obs2','perc2','match2','Sc','Spv'))
			#wLut.write('%11s %11s %11s | %11s %11s %11s \n' % ('avg1','med1','lsqsol1','avg2','med2','lsqsol2'))
			
			wLut.write('%-40s %40s\t%5s %7s %7s\t%5s %7s %7s\t%5s %8s\t' % ('---','---','obs1','perc1','match2','obs2','perc2','match2','Sc','Spv'))
			wLut.write('%11s %11s %11s\t%11s %11s %11s \n' % ('avg1','med1','lsqsol1','avg2','med2','lsqsol2'))
			SLF = float(len(self.D.samples)) 

			for fi in range(len(self.D.features)): 
				self.progress.mark() 
				f_key = dd(list)
				m_key = dd(list)
				f_cnts = []   
				x_cnts = dd(list) 
				y_cnts = dd(list) 

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
					RY = np.mean([y/x for x,y in zip(xC,yC)]) 
					ZY = np.median([y/x for x,y in zip(xC,yC)]) 
					rawA = np.linalg.lstsq(xCA[:,np.newaxis],yCA,rcond=None)[0][0] 
					rawB = np.linalg.lstsq(yCA[:,np.newaxis],xCA,rcond=None)[0][0] 

					
					
#					wSol.write('%-40s %40s | %5d %10.4f %10.4f | %6.4f %8.2e | %10.5f %10.5f \n' % (iName,jName,iLen,iRT,matA,SR2,SRv,RM,rawA))
#					wSol.write('%-40s %40s | %5d %10.4f %10.4f | %6.4f %8.2e | %10.5f %10.5f \n' % (jName,iName,jLen,jRT,matB,SR2,SRv,RY,rawB))


					#wLut.write('%-40s %40s | %5d %7.4f %7.4f | %5d %7.4f %7.4f | %5.3f %8.1e | %11.4f %11.4f %11.5f | %11.4f %11.4f %11.5f\n' % (iName,jName,iLen,iRT,matA,jLen,jRT,matB,SR2,SRv,RM,ZM,rawA,RY,ZY,rawB))
			

					wLut.write('%-40s %40s\t%5d %7.4f %7.4f\t%5d %7.4f %7.4f\t' % (iName,jName,iLen,iRT,matA,jLen,jRT,matB))
					wLut.write('%5.3f %8.1e\t%11.4f %11.4f %11.5f\t%11.4f %11.4f %11.5f\n' % (SR2,SRv,RM,ZM,rawA,RY,ZY,rawB))


					self.s_key[fj].append((fi,iLen,matA,SR2,SRv,RM,rawA))
					self.s_key[fi].append((fj,jLen,matB,SR2,SRv,RY,rawB))



#				for fp,fl in f_key.items():
#					fLen,fMed,fMean = len(fl), np.median(fl), np.mean(fl) 
#					if RSTAT == 'MED':	self.r_key[fp] = fMed
#					else:			self.r_key[fp] = fMean
#					c1 = len(self.D.features[fp[0]].cnts)
#					c2 = len(self.D.features[fp[1]].cnts)
#					self.HOUSEKEEPING[fp[0]] = True 
#					self.HOUSEKEEPING[fp[1]] = True
#					wMean.write('%-40s %40s %5d %5d | %10d %25f %25f\n' % (self.D.features[fp[0]].name,self.D.features[fp[1]].name,c1,c2,fLen,fMed,fMean))


		wLut.close() 
		wLut = open(self.args.prefix+'.ratioSolutions.out')
		
		self.read_solution(wLut)

		

		return self





	def predict_sample_values(self,BREAKLEN=500,MINLEN=100):

#		for fi in self.s_key.keys():
#			self.s_key[fi].sort(key = lambda X: X[3],reverse=True)


		for si,s in enumerate(self.D.samples): 	
			S_ERRS = dd(list) 
			for i,f in enumerate(self.D.features):
				S_TRUE,M_PRED,R_PRED,S_PRED,R_MINI,S_MINI,R_SPEC,S_SPEC = 0,[],[],[],[],[],[],[]
				if i in s.cnts: S_TRUE = s.cnts[i] 
					
				
				
				for fj in self.solver_key[i].keys(): 

					if fj not in s.cnts: continue 

					


					fSol =  self.solver_key[i][fj] 
					fAvg =  self.ratio_key[i][fj] 

					R_PRED.append(s.cnts[fj] / fAvg) 
					S_PRED.append(s.cnts[fj] * fSol) 



				
				if S_TRUE != 0 and len(S_PRED)>1: 
					Rmean,Smean = np.mean(R_PRED),np.mean(S_PRED)
					for pr,pn in zip([Rmean,Smean],['rMean','sMean']):
						S_ERRS[pn].append((S_TRUE-pr)*(S_TRUE-pr))			



			print s.name,'rMean','sMean',np.mean(S_ERRS['rMean'])**0.5,np.mean(S_ERRS['sMean'])**0.5




	def impute_sample_values(self,BREAKLEN=500,MINLEN=100):
		wImpute = open(self.args.prefix+'.imputed.cnts','w') 
		wInfer = open(self.args.prefix+'.inferred.cnts','w') 
		wTotals = open(self.args.prefix+'.inferred.totals','w') 
		wPerf = open(self.args.prefix+'.feature.perfs','w') 

   		self.progress.start_minor('Calculating Ratio Values',len(self.D.features),True)
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
#				print f.name, preds, Sp[0],Sp[1],Rp[0],Rp[1],'|',Se,R

		for s in self.D.samples: 
			wTotals.write('%s True/Infer/Infer-Missing/Impute/Impute-Missing %d %d %d %d %d\n' % (s.name,sum(s.cnts.values()),sum(s_inferred[s]),sum(s_infer_missing[s]),sum(s_imputed[s]),sum(s_impute_missing[s]) ))

		































































	def predict_missing_sample_values(self,STAT='MEAN'): 

		EXP_TOTALS = []
		self.expected_totals = [] 
		s_totals = [sum(s.cnts.values()) for s in self.D.samples] 		
		w = open(self.args.prefix+'_mR_effective.out','w')
                w.write('%-40s %10s %20s %20s\n' % ('---','idx','actual_total','inferred_total'))


		print self.HOUSEKEEPING
			
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



















