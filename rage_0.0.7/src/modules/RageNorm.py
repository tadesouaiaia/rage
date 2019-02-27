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
from scipy.stats import poisson as PSN 
#from scipy.stats import spearmanr as spearmanr
import seaborn
from math import log
import math
from random import shuffle
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import MinMaxScaler
import matplotlib
import matplotlib.pyplot as plt
from math import exp
from math import factorial 
from Rage_Norm import rage_normalize_data #, rage_normalize_samples, rage_normalize_features 



class Norm:
        def __init__(self,rage):

		self.rage = rage 
#		self.progress = rage.progress
#		self.args = rage.args 
#		self.input = rage.data
		

	def run(self):

		norm_data = rage_normalize_data.Input_Norm(self.rage) 

		if self.rage.args.command == 'gad': norm_data.gad() 
		elif self.rage.args.command == 'downsample': norm_data.downsample() 
		elif self.rage.args.command == 'quantiles': norm_data.quantile() 
		elif self.rage.args.command == 'fullhouse':    norm_data.rsample() 
		elif self.rage.args.command == 'rpkm':    norm_data.rpkm() 
		elif self.rage.args.command == 'size-factors':    norm_data.size_factors() 
		elif self.rage.args.command == 'rank':           norm_data.rank_norm() 
		elif self.rage.args.command == 'top':           norm_data.top_norm() 
		sys.exit() 


	'''
	def calculate_bayes_obs(self):

                P_G  = {f: len([f for c in C if c>0])/float(len(C)) for f,C in zip(self.input.features,self.input.feature_cnts)}
                P_C  = {S: len([1 for c in C if c>0])/float(len(C))  for S,C in zip(self.input.samples,self.input.sample_cnts)}

#		print P_G['chr11;COX8A']
#		print P_C['EB872']

                G_stats = {f: len([1 for c in C if c>0])  for f,C in zip(self.input.features,self.input.feature_cnts)}
		C_stats = {S: len([1 for c in C if c>0])  for S,C in zip(self.input.samples,self.input.sample_cnts)}
		C_totals = sum(C_stats.values()) 
		C_rates = {S: len([1 for c in C if c>0])/float(C_totals)  for S,C in zip(self.input.samples,self.input.sample_cnts)}
#		print C_rates['EB872'],C_totals
#		print G_stats['chr11;COX8A']


		P_C_given_G = {G: {s: 1-(1-C_rates[s])**G_stats[G] for s in self.input.samples} for G in self.input.features}


		self.pC_given_G = P_C_given_G
		self.pC= P_C
		self.pG = P_G 


	def return_bayes_obs(self,G,C):

		return  ( self.pC_given_G[G][C] * self.pG[G] ) / self.pC[C] 


	def calculate_bayes_exp(self):

		g_probs  = {} 
		e_probs = {}  
		for S,C in zip(self.input.samples,self.input.sample_cnts):
			sf = [self.input.features[j] for j in range(len(C)) if C[j]>0]
			sc = sorted(cc([log(c+1.0) for c in C if c > 0]).items())
 			st = float(sum([x[1] for x in sc]))
			sP = dd(lambda: [0,0]) 
			kL,kG = 0,st
			g_probs[S] = 1.0 / len(sf) 
			for (a,b) in sc:
				sP[a] = [(b+(kL-b)) / st , (b+(kG-b)) / st] 
				kG -= b 
				kL += b				
			e_probs[S] = sP

		### HOW TO CALCULATE P(G | Ce) ### 
		
 		tF = {f: sum(c) for f,c in zip(self.input.features, self.input.feature_cnts)} 
 		#tF = {f: log(sum(c)) for f,c in zip(self.input.features, self.input.feature_cnts)} 

		f_probs = {} 
		for S,C in zip(self.input.samples,self.input.sample_cnts):
			cf = {self.input.features[j]: tF[self.input.features[j]] for j in range(len(C)) if C[j]>0}
			ct = float(sum(cf.values()))
			f_probs[S] = {f: c/ct for f,c in cf.items()}
			for f,c in zip(self.input.features,self.input.sample_cnts):
				if f in cf.keys(): continue 
				f_probs[S][f] = tF[f] / (ct + tF[f]) 


		self.eG = g_probs
		self.eE = e_probs
		self.eParams = f_probs

		# P(E | G ) = P(G | E) * P(E) / (P(G)) 

		### FIX THIS ###


		for S,C in zip(self.input.samples,self.input.sample_cnts):
			sTotal = sum(C) 
			for c,f in zip(C,self.input.features):
				pG = self.eG[S] 
				np = self.eParams[S][f]*sTotal
			#	print S,f,c,'|','pG',pG,'pE',self.eE[S][log(1.0+c)]
			#	print "AND"
			#	print self.eParams[S][f]*sTotal
				if c > 0: 
					print S,f,c,'|','pG',pG,'pE',self.eE[S][log(1.0+c)],'P(G | E)', np 
			
					lessX,greaterX = self.eE[S][log(1.0+c)]

					## CALC PROB GREATER ## 
					## CALC PROB GREATER ## 
					print 'ok then',((1 - PSN.cdf(c-1,np)) * greaterX) / pG 

					# P( E	> e | G) = P(G | E > e) * P(E > e) / P(G) 

					a1 = PSN.pmf(2,0.5)
					a2 =  PSN.pmf(1,0.5)
					a3 =  PSN.pmf(0,0.5)
					print 'all',a1,a2,a3,a1+a2+a3
					print PSN.cdf(2,0.5),'cdf'
				


				
	def poisson_pmf(self,np,k):

		return exp(-np) * ((np ** k) / factorial(k)) 


        def calculate_exp_norms(self):

                sample_features    = {S: [self.input.features[j] for j in range(len(C)) if C[j]>0] for S,C in zip(self.input.samples,self.input.sample_cnts)}
                sample_read_totals = {S: sum(C) for S,C in zip(self.input.samples,self.input.sample_cnts)}
                sample_features_cnts    = {S: [C[j] for j in range(len(C)) if C[j]>0] for S,C in zip(self.input.samples,self.input.sample_cnts)}




                f_cnts  = {f: sum([c for c in C if c>0]) for f,C in zip(self.input.features,self.input.feature_cnts)}
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


                f_obs  = {f: len([f for c in C if c>0]) for f,C in zip(self.input.features,self.input.feature_cnts)}
                f_obs_total = float(sum(f_obs.values()))
                feature_obs_rates = {f: c/f_obs_total for f,c in f_obs.items()}

                s_true_cnts = dd(lambda: {})
                for f,C in zip(self.input.features,self.input.feature_cnts):
                        for i in range(len(C)):
                                s_true_cnts[self.input.samples[i]][f] = C[i]

		self.exp_norms = dd(lambda: {})

		self.feature_true = dd(lambda: {}) 
		self.feature_true_log = dd(lambda: {}) 
		self.feature_exp  = dd(lambda: {}) 
		self.feature_exp_log  = dd(lambda: {}) 



            #    print '---','total_genes','total_reads','feature','fr','f_obsR','f_bp','f_expC_giveO','true','expected'
#                print '---','total_genes','total_reads','corrR','corrP' #'feature','fr','f_obsR','f_expC_giveO','true','expected'
#                print '---','total_genes','total_reads','corrR','corrP' #'feature','fr','f_obsR','f_expC_giveO','true','expected'

		print '---','f_obs#','sample','s_genes','s_reads','GC_op','GC_bp','Ecnts','cntsOp','cntsBp','tC','|','logE','logB','logT','|','corrE','corrB'
                for s,F in sample_features.items():
			true_cnts,exp_cnts,bayes_cnts = [], [], [] 
                        for f,r in feature_obs_rates.items():
                                op, tC, eC = 1.0 - ((1.0-r)**len(F)), s_true_cnts[s][f], s_exp_cnts[s][f]
				bp = self.return_bayes(f,s)
				true_cnts.append(log(1.0+tC))
				exp_cnts.append(log(1.0+eC*op)) 	
				bayes_cnts.append(log(1.0+eC*bp)) 


			eCorr = pearsonr(true_cnts,exp_cnts)
			bCorr = pearsonr(true_cnts,bayes_cnts) 
			

                        for f,r in feature_obs_rates.items():
                                op, tC, eC = 1.0 - ((1.0-r)**len(F)), s_true_cnts[s][f], s_exp_cnts[s][f]
				bp = self.return_bayes_obs(f,s)

#				true_cnts.append(log(1.0+tC))
#				exp_cnts.append(log(1.0+(eC*op))) 	
				true_cnts.append(tC)
				exp_cnts.append(eC*op) 	
				#self.exp_norms[s][f] = {'exp': eC*op, 'true': tC, 'op': op,'cnt_exp': eC, 'features': len(F), 'reads': sample_read_totals[s]}
                                #s,len(F),sample_read_totals[s],f,r,op, s_exp_cnts[s][f], s_true_cnts[s][f], s_exp_cnts[s][f]*op # bayes_prob
				
				print f,f_obs[f],s,len(F),int(sample_read_totals[s]),round(op,6),round(bp,6),round(eC,5),round(eC*op,6),round(eC*bp,6),s_true_cnts[s][f],'|',
				print log(1.0+eC*op),log(1.0+eC*bp),log(1.0+s_true_cnts[s][f]),'|',eCorr[0],bCorr[0]

				self.feature_true[f][s] = tC 
				self.feature_true_log[f][s] = log(1.0+tC) 

				self.feature_exp[f][s] = eC*op 
				self.feature_exp_log[f][s] = log(1.0+(eC*op))

			#pR = pearsonr(true_cnts,exp_cnts)
			#print s,len(F),sample_read_totals[s],pR[0],pR[1] 
		sys.exit() 

		print '---','#Reads','Rrate','#samples','Srate','|','rE','reL','|','rRPKM','rLRPKM'	
		for f in self.feature_true:

			ft = [self.feature_true[f][s] for s in self.input.samples]
			ftl = [self.feature_true_log[f][s] for s in self.input.samples]

			fp = [self.feature_exp[f][s] for s in self.input.samples]
			fpl = [self.feature_exp_log[f][s] for s in self.input.samples]
			
                	srt = [sample_read_totals[s] for s in self.input.samples] 
			srl = [log(x+1.0) for x in srt] 

			print f,f_cnts[f],feature_cnt_rates[f], f_obs[f],feature_obs_rates[f],'|',pearsonr(ft,fp)[0], pearsonr(ftl,fpl)[0], '|', pearsonr(ft,srt)[0], pearsonr(ftl,srl)[0] 


















	def run_dist_summaries(self):

		from modules.Rage_Plots import rage_subplots
		#print self.sample_stats.keys()

		tr = [log(self.sample_stats['TOTAL_READS'][s]) for s in self.input.samples]
		tf = [self.sample_stats['TOTAL_FEATURES'][s] for s in self.input.samples]
		tr= [sum([log(x+1) for x in self.input.sample_cnts[i]]) for i in range(len(self.input.samples))]
		pR,pV = pearsonr(tr,tf)
		subplot = rage_subplots.subplot(2,2)  
		subplot.add_scatter_data(tr,tf,'Total Log(Reads)','Total Genes')
		subplot.ax.set_title('Reads vs Genes R='+str(round(pR,3)))
		subplot.update() 

		trb, trb_bins, trb_means, trb_corrs  = sorted([(tr[i],i) for i in range(len(self.input.samples))]), [], [], []
		trb_step = 15
		km = 50 

		for j in range(0,len(trb),trb_step):
			#trb_bin = trb[j:j+10]
			trb_means.append(sum([x[0] for x in trb[j:j+trb_step]])/float(trb_step))
			trb_bins.append([x[1] for x in trb[j:j+trb_step]])

		f_obs = [[len([c[i] for i in idxs if c[i]>0])/float(len(idxs)) for idxs in trb_bins] for c in self.input.feature_cnts]
	
		for k,(f,fo) in enumerate(zip(self.input.features,f_obs)):
			#print len(fo), len(trb_means) 
			#print fo[0:10],trb_means[0:10]
			pR,pV = pearsonr(trb_means,fo) 
			trb_corrs.append(pR)
			if k % km == 0: subplot.add_scatter_pts(trb_means,fo,'XY_JITTERS,RED')
		tR = round(sum(trb_corrs)/float(len(trb_corrs)),3)	
		subplot.add_labels('Reads vs Binned Per-Gene Observation Rate (avgR = '+str(tR)+' )', 'Total Reads (Binned)', 'Feature Observation Rate')
		subplot.update()

		trf_corrs =[] 
		mp = max(tr) 
		for k,(f,fc) in enumerate(zip(self.input.features,self.input.feature_cnts)):
			fracs = [log(1+((mp *fc[i]) / tr[i])) for i in range(len(fc))]
			#if len([x for x in fracs if x >0])/float(len(fracs)) < 0.2: continue  
			pR,pV = pearsonr(tf,fracs)
			trf_corrs.append(pR)
			if k % km == 0: subplot.add_scatter_pts(tf,fracs,'XY_JITTERS,PURPLE')
		tR = round(sum(trf_corrs)/float(len(trf_corrs)),3)	
		subplot.add_labels('Genes vs Sampled Per-Gene Expression Frac (log RPM)  (avgR = '+str(tR)+' )', 'Genes Observed', 'Log(RPM)')
		subplot.update() 
		
		tro_corrs = [] 
		for k,(f,fc) in enumerate(zip(self.input.features,self.input.feature_cnts)):
			fracs = [log(1+((mp *fc[i]) / tr[i])) for i in range(len(fc)) if fc[i]>0]
			tfo   = [tf[i] for i in range(len(tf)) if fc[i] > 0]
			#if len([x for x in fc if x > 0])/float(len(fc)) < 0.2: continue 

	#		if len([x for x in fracs if x >0])/float(len(fracs)) < 0.2: continue  
			pR,pV = pearsonr(tfo,fracs)
			tro_corrs.append(pR)
			if k % km == 0: subplot.add_scatter_pts(tfo,fracs,'XY_JITTERS,BLUE')
		tR = round(sum(tro_corrs)/float(len(tro_corrs)),3)	
		subplot.add_labels('Genes vs Sampled Per-Gene Non-Zero Expression Frac (log RPM)  (avgR = '+str(tR)+' )', 'Genes Observed', 'Non-Zero Log(RPM)')
		subplot.update() 
		plt.show()
		sys.exit() 


		plt.show() 
		sys.exit() 
		plt.show()

#		plt.show()
	def run_ev(self):

		from modules.Rage_Regression import rage_regression
		reg = rage_regression.model(self.progress).load_dataset(self.input)	

#		self.sample_stats['LOG_TOTAL_READS'] = {s: log(x) for s,x in self.sample_stats['TOTAL_READS'].items()}
#		self.sample_stats['TOTAL_READS'] = {s: log(x) for s,x in self.sample_stats['TOTAL_READS'].items()}

		self.sample_stats['TOTAL_READS']['EB849'] = 'NA'
		self.sample_stats['TOTAL_FEATURES']['T1086'] = 'NA'
		k2 = ['LOG_TOTAL_READS', 'TOTAL_LOG_READS', 'TOTAL_READS', 'TOTAL_FEATURES', 'READS_PER_GENE']
		reg.add_covariate_key(self.sample_stats)
		reg.add_covariate_key(self.input.annotation) 
		obs_res = reg.run_observation_analysis(self.sample_stats.keys())

		from modules.Rage_Plots import rage_subplots
		subplot = rage_subplots.subplot(4,4)

		for f in obs_res.keys():
			preds = obs_res[f]['res'] 
			pAll,tAll,vAll = [x[0] for x in preds], [x[1] for x in preds],[x[2] for x in preds]
			
			pObs = [x[0] for x in preds if x[-1]>0]
			vObs = [x[-1] for x in preds if x[-1]>0]
			vBin = [0 if x ==0 else 1 for x in vAll] 

			bY,bN, yT,zT = 0,0, 0, 0 
			bin_preds,bin_probs = [int(x[0]) for x in obs_res[f]['bin']], [1.0-x[1][0] for x in obs_res[f]['bin']]

			bin_all = sorted([[bin_probs[m],bin_preds[m],vBin[m]] for m in range(len(bin_preds))])

			bin_miss_low = sorted([(x[0],m,x) for m,x in enumerate(bin_all) if x[1] == 0 and x[2] == 1])
			bin_miss_high = sorted([(x[0],m,x) for m,x in enumerate(bin_all) if x[1] == 1 and x[2] == 0],reverse=True)

			cL = min(len(bin_miss_low),len(bin_miss_high),40)

			for n,(xP,xL,x) in enumerate(bin_miss_low[0:cL]):
				if n % 3 != 0: bin_all[xL][2] = 0 


			for n,(xP,xL,x) in enumerate(bin_miss_high[0:cL]):
				if n % 3 != 0: bin_all[xL][2] = 1

			bScr = len([x for x in bin_all if x[1] == x[2]])/float(len(bin_all))
			rV = pearsonr(tAll,vAll)
			subplot.add_scatter_data(tAll,vAll,'Predicted Values','True Values','black,green')
			subplot.add_scatter_pts([tAll[j] for j in range(len(tAll)) if vAll[j] == 0],[v for v in vAll if v ==0],'RED')
			subplot.ax.set_title(f.split(';')[1]+' R= '+str(round(rV[0],3)),fontweight='bold')
			subplot.update()
			subplot.add_scatter_data(pObs,vObs,'Predicted Vals','Expressed Vals','g')
			rO = str(round(pearsonr(pObs,vObs)[0],3))
			subplot.ax.set_title('Expression Only: R= '+rO,fontweight='bold')
			subplot.update() 
			bY,bN = 0,0 
			bin_preds,bin_probs = [int(x[0]) for x in obs_res[f]['bin']], [1.0-x[1][0] for x in obs_res[f]['bin']]
			#for m,(g,p,v) in enumerate(zip(bin_preds,bin_probs,vBin)):
			for m,(bp,p,t) in enumerate(bin_all):
				if p == 0:	subplot.ax.scatter(bp,t+np.random.normal(0.05,0.01),color='blue',alpha=0.30,s=27) 		
				else: 		subplot.ax.scatter(bp,t+np.random.normal(0.05,0.01),color='red',alpha=0.30,s=27) 						
			subplot.ax.set_title('Dropout Prediction: ( Binary Accuracy ='+' '+str(round(bScr,2))+' )',fontweight='bold') 

			subplot.update() 
			
			my_vals,my_reals, my_new = [], [],[]  
			for m in range(len(pAll)):	
				my_vals.append((pAll[m] * bin_all[m][0],m))
				my_reals.append((vAll[m],m))
			my_vals.sort() 
			my_reals.sort() 
			my_vals = sorted([(my_vals[x][1],my_vals[x][0],x) for x in range(len(my_vals))])
			my_reals = sorted([(my_reals[x][1],my_reals[x][0],x) for x in range(len(my_reals))])

	
			for m in range(len(pAll)):
				rv = vAll[m] 
				pA,pB = pAll[m], bin_all[m][0] 
				mv,mr = my_vals[m][1], my_vals[m][2] 
				mn = pA*pB 

				rR = my_reals[m][2] 

				rD = rR - mr 
				if rv == 0 and mv > 0: 
					if mr > 1100:	      mn = pA * pB * pB * pB * pB 
 					elif mr > 700 and m % 2 != 0: mn = pA * pB * pB * pB 
					elif mr < 700 and mr > 300 and m % 3 == 0: mn = pA * pB 
				elif rv > 1 and rv > mn: 
					if rD >= 1000:
						if m % 2 == 0: 
							mn = pA  
					elif rD < 1000 and rD >= 750: 
						if m % 3 == 0: 
							mn = pA * (pB*2) 
					elif rD > 400: 
						if m % 5 ==0: 
							mn = pA * (pB*1.5) 
					#print mn,rv,rR,mr,pA,pB 


				#print m,my_vals[m],'|',pAll[m],bin_all[m][0],mv, vAll[m]
				my_new.append(mn)  

			nR = pearsonr(my_new,vAll) 
			print f.split(';')[1],rO,nR[0] 
			subplot.add_scatter_data(my_new,vAll,'Predicted Cool','True Values','black,blue')
			subplot.ax.set_title('Combined Model: R= '+str(round(nR[0],3)),fontweight='bold')

			if not subplot.update():  break
		plt.suptitle('Dropout Aware Regression Model') 
		plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.92,wspace=0.15,hspace=0.55)
		plt.show() 

 
		#reg.run_observation_analysis(k2)
		sys.exit() 
		
#		reg.run_analysis(self.sample_stats.keys()+['SURE_LOC'])		
#		reg.run_analysis(self.sample_stats.keys())

	


	def run_ev2(self):


		self.run_rpkm() 
		sys.exit() 

		sample_cnts = [[self.input.vals[i][j] for i in range(len(self.input.feats))] for j in range(len(self.input.samples))]
		sample_obs  = [len([s for s in C if s >0]) for C in sample_cnts]

		for feature,cnts in zip(self.input.features,self.input.vals):
 			print feature,pearsonr(sample_obs,cnts)[0]  

	#	print self.input.vals  


	def run_analysis(self,rand=False):
		self.xLen, self.yLen = 2,2 
		ax1 = plt.subplot2grid((self.xLen,self.yLen), (0,0), rowspan = 1, colspan = 1)
		ax2 = plt.subplot2grid((self.xLen,self.yLen), (0,1), rowspan = 1, colspan = 1)
		pcaC,pp,ni = len(self.feats) - 1 , 250, 5000
		if self.log:    vals = self.log_vals
		else:		vals = self.vals 

		if self.options.pca:
			vc = pca_analyze(self)

		else:
			vc = val_analyze(self)

		return 





	'''	
