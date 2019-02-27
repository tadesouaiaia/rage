#!/usr/bin/env python

import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
import scipy.stats as stats
from scipy.stats import variation as coVar 

from scipy.stats import ttest_ind
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
from scipy.stats import chisquare

from Rage_Classify import rage_classify
from Rage_Regression import rage_regmodels as rt
from Rage_IO import rage_outputs
from Rage_Plots import rage_subplots
#from Rage_Plots import  rage_dimplots as dplot 

#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 



def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))




class Classify:
        def __init__(self,rage):

		self.rage, self.progress, self.options  = rage, rage.progress, rage.args


	def run(self):
		R = self.rage 

		self.progress.start_major('Running Classifier Analysis') 


		if R.args.command == 'groups':

			
			self.options.id = ["*".join(self.options.id)]

			if self.options.raw: 
				self.D = self.rage.data.filter_samples_by_attributes(self.options.id,self.options.covariates).scale_and_transform(LOG=False)
			else: 
				self.D = self.rage.data.filter_samples_by_attributes(self.options.id,self.options.covariates).scale_and_transform(LOG=True)
                        self.V = self.D.set_sample_variables(self.options.id, self.options.covariates)
                        self.X, self.Xp,self.Xc = self.V.select_variables(self.V.variables), self.V.select_variables(self.V.predictors), self.V.select_variables(self.V.covariates)

     #                   if self.options.zero_prob != None and self.options.zero_prob[0:3].upper() == 'OBS': self.X.zp =  np.array([[( 1 - (len(s.cnts)  / float(len(self.D.features))) )] for s in self.D.samples])





			self.F, self.Y = [],[]
                        for f in self.D.features:
                                y,y_obs,y_name = [s.cnts[f.idx] for s in self.D.samples], len(f.cnts)/ float(len(self.D.samples)) , f.name
                                if y_obs > 0.050 and len(set(y)) > 5:
                                        self.Y.append(y)
                                        self.F.append(y_name)
                
                        self.D.samples.create_plot_labels(self.options)


	



			ID_TYPES = list(set([self.D.samples.attribute_class[ID] for ID in self.options.id]))
			if ID_TYPES == ['binary']: 
				classifier = rage_classify.BinaryClassifier(self.options,self.D) 

				
				if len(self.options.covariates) > 0:

					#Mc  = rt.RegModel(self.options,self.Xc,True).run(self.Y).aggregate(True)

					M = rt.RegModel(self.X,'OLS-LOG',self.options,self.progress,True).run(self.Y,self.F) #
					M_resids, C_resids = M.get_resids(COVAR=True)

					Y = C_resids
 				else:
					Y = self.Y 




				self.progress.start_minor('Starting Classification')
				id_locs, id_names =  self.V.select_labels(self.options.id[0])

				


				
				
				if 'UNK' in [x.split('~')[-1][0:3].upper() for x in id_names]:
					classifier.logitU(id_locs,id_names,Y,self.F,self.options.leave) 
				elif self.options.markeronly: 
					

					classifier.findMarkers(id_locs,id_names,Y,self.F,MFC=1.3,MFF=1.1) 
					cY,cF = [y for y,f in zip(Y,self.F) if classifier.marker_stats[f]['MARKER']],[f for y,f in zip(Y,self.F) if classifier.marker_stats[f]['MARKER']]
					classifier.logitR(id_locs,id_names,cY,cF)  
				else:
					classifier.logitR(id_locs,id_names,Y,self.F,self.options.leave)  
				#classifier.logitR(id_locs,id_names,cY,cF)  
				#classifier.logitR(id_locs,id_names,Y,self.F)  
			

			



		elif R.args.command in ['markers']:
			R.progress.start_major('MarkerSearch')
			R.progress.start_minor('Calculating Marker Stats',R.data.features.len)

			det_trends = detail_trends(R.data.features,R.data.samples,R.args,R.progress)
			sys.exit() 

			
			sys.exit()



		elif R.args.command == 'batches':
			print 'cool'
			sys.exit() 


		elif R.args.command == 'logistic':
			print 'logistic'
			sys.exit() 
		else:
			print 'no'
			sys.exit() 
			R.progress.start_minor('Summarizing Dimensional Reduction',R.data.features.len)
			pca, tsne, kca = make_dr_plots(R,'features') 
			rage_outputs.column_coefs(R.args).write(pca['coefs'],R.data.samples,{'suffix': 'PCAcoeffs.samples.out','width': 20})
			rage_outputs.dr_matrix(R.args).write(kca['run'].X_transformed_fit_,R.data.features,{'suffix': 'kcaFit.features.out'})
			for dr in [pca,tsne,kca]:	rage_outputs.dr_pts(R.args).write(dr['pts'],R.data.features,{'suffix': dr['type']+'.featurePts.out'})


				 
				


		sys.exit() 






















			#if self.args.pca: sample_summary.make_pca_and_tsne_plots() 
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


class Groups:
        def __init__(self,Y,gkey):



		self.Y = Y 

		self.pv_max = 0.01
		self.min_dist = 4
		self.colors = {} 
		self.members = {}
		for k,idxs in gkey.items():
			self.members[k] = [Y[i] for i in idxs]
			self.colors[k]  = self.members[k][0].notes['color'] 
		self.names = gkey.keys() 


	def updateFeature(self,x):

		self.vals, self.svals, self.ranks, self.obs, self.means, self.qts, self.smeans, self.size   = {} , {}, {}, {}, {} , {} , {}, {} 

		svals = scale_vals([log(y.cnts[x.idx]+1.0) for y in self.Y],0,100)

		for k,Y in self.members.items():

			self.vals[k]  = [log(1.0+y.cnts[x.idx]) for y in Y]
			self.ranks[k]  = [y.notes['rank'][x.idx] for y in Y]
			self.obs[k] = len([z for z in self.vals[k] if z > 0])/ float(len(self.vals[k]))
			self.means[k] = np.mean(self.vals[k])
			self.qts[k]    = np.percentile(self.vals[k],90) 
				
			self.svals[k] = [svals[y.idx] for y in Y]
			self.smeans[k] = np.mean(self.svals[k])

			self.size[k] = len(self.vals[k]) 


		
		ms = sorted([(self.means[k],k) for k in self.members.keys()])
		self.top = ms[-1][1]
		self.bottom = ms[0][1]



	def marker_filter(self,min_obs = 0.2, pass_obs = 0.5, pass_frac = 1.5):
		self.check_pvs()

			
		if min([p[1] for p in self.tts.values()]) > self.pv_max: return False
		my_obs = self.obs.values() 


		if max(my_obs) > pass_obs: return True
		if max(my_obs) < min_obs: return False
		if max(my_obs) / (0.001+min(my_obs)) > pass_frac: return True 
		return False


	def check_pvs(self):
		self.tts = {} 
		for k,kV in self.vals.items():
			nV = [a for  b in [self.vals[x] for x in self.vals.keys() if x != k] for a in b]
			nvM = np.mean(nV)
			if self.means[k] > nvM: fc = self.means[k] / nvM
			else:			fc = nvM / self.means[k] 
			self.tts[k] = (fc,stats.ttest_ind(kV,nV)[1]) 

	def check_dists(self):
		keys = self.vals.keys()
		means = [self.smeans[k] for k in keys]
		obs   = [self.obs[k] * 100 for k in keys]
		ranks = [np.mean(self.ranks[k]) for k in keys]
		self.dist_key = dd(list) 
		return True 
		for i in range(len(keys)):
			iPt = np.array([means[i],obs[i],ranks[i]])
			for j in range(i+1,len(keys)):
				jPt = np.array([means[j],obs[j],ranks[j]])
				D = round(np.linalg.norm(iPt-jPt),3) 
				self.dist_key[keys[i]].append(D)
				self.dist_key[keys[j]].append(D)	
		for k,D in self.dist_key.items():
			if min(D) > self.min_dist: return True 	
		return False 


	def check_pts(self):

		pts = [x[1] for x in sorted([a for b in  [[(kx,k) for kx in self.vals[k]] for k in self.vals.keys()] for a in b])]

		vals = [x[0] for x in sorted([a for b in  [[(kx,k) for kx in self.vals[k]] for k in self.vals.keys()] for a in b])]
		aS = float(len(pts)) 

		if vals[len(pts) - self.size[self.top]] == 0:
			pres = len([v for v in vals if v>0]) 
			tS,bS,aS = pres, 1-pres, float(len(pts))
		else:
			tS,bS,aS = self.size[self.top], self.size[self.bottom], float(len(pts))


		tPts =  pts[len(pts) - tS::]
		bPts =  pts[0:bS]

		tEH, tEL =  int(( tS / aS ) * tS ), int((tS/aS) * bS) 
		bEH,bEL =  int(( bS / aS ) * tS ), int((bS/aS)* bS )
		tHO, tLO =  cc(tPts)[self.top], cc(bPts)[self.top] 
		bHO, bLO = cc(bPts)[self.bottom], cc(bPts)[self.bottom]
		hi_chi =  chisquare([tHO,tS-tHO],f_exp=[tEH,tS-tEH])[1]
		lo_chi =  chisquare([bLO,bS-bLO],f_exp=[bEL,bS-bEL])[1]

		return tHO-tEH, bLO-bEL,hi_chi,lo_chi

	def check_halfs(self):

		keys = self.vals.keys() 
		vals = [x[0] for x in sorted([a for b in  [[(kx,k) for kx in self.vals[k]] for k in self.vals.keys()] for a in b])]
		pres = len([v for v in vals if v>0]) / float(len(vals))
		aS = float(len(vals))

		if pres < 0.04: return 1.0,1.0 

		if pres > 0.4: pres = 0.4

		pnot = 1-pres

		ec = [int((self.size[k]/aS)*(aS*pres)) for k in keys]
		enc = [int((self.size[k]/aS)*(aS*pnot)) for k in keys]
		 		

		pts = [x[1] for x in sorted([a for b in  [[(kx,k) for kx in self.vals[k]] for k in self.vals.keys()] for a in b],reverse=True)][0:int(aS*pres)]
		ptsN = [x[1] for x in sorted([a for b in  [[(kx,k) for kx in self.vals[k]] for k in self.vals.keys()] for a in b])][0:int(aS*pnot)]
	
		pc = cc(pts) 	
		obs = [pc[k] for k in keys] 

		pn = cc(ptsN) 
		obsN = [pn[k] for k in keys]

			
		P1 = chisquare(obs,f_exp=ec)[1]	
		P2 = chisquare(obsN,f_exp=enc)[1]	


		return P1,P2	













def detail_trends(X,Y,options,progress,min_obs = 0.2):



	seaborn.set(rc={'axes.facecolor':'gainsboro', 'figure.facecolor':'whitesmoke'})

	min_max_diff = 0.75
	min_proj_diff = 0.5
	min_mean_max_diff = 0.5 

	min_max_diff = 0.45
	min_proj_diff = 0.05
	min_mean_max_diff = 0.25 
	min_group_size = 20

	a_key = {a: cc([y.attributes[a] for y in Y]) for a in Y.attributes}


	for y in Y:
		y.notes['obs'] = len(y.cnts) / float(len(X))
		y.notes['total'] = log(1.0+sum(y.cnts.values()))
		y.notes['rank'] = dd(int) 
		for i,(c,f) in enumerate(sorted([(y.cnts[f],f) for f in y.cnts])):	y.notes['rank'][f] = round(100*((i+1) / float(len(y.cnts))),3)
#		for i,(c,f) in enumerate(sorted([(y.cnts[f],f) for f in y.cnts])):	y.notes['rank'][f] = i+1 

		

	interests = ['CELLSPAN'] 

	for interest in interests:



		if Y.attribute_class[interest] != 'binary': continue 

		

		i_names = [a for a,b in a_key[interest].items() if b> min_group_size]


		g_key = {g : [i for i in range(len(Y)) if Y[i].attributes[interest] == g] for g in i_names if g not in ['UHR','NA','IZ_SP','IZ_SVZ','DG_SUB']}

		print g_key.keys() 

		if len(g_key.keys()) < 2: continue 


	#	print g_key
	#	g_labels = g_key.keys()
		Y.colorize(interest,min_group_size) 
	#	g_colors = [Y[g_key[g][0]].notes['color'] for g in g_labels]

		
		groups = Groups(Y,g_key)



		for x_idx,x in enumerate(X):
			groups.updateFeature(x) 

	

			if not groups.marker_filter(min_obs = 0.001, pass_obs = 0.2, pass_frac = 1.2): continue
			s1,s2 = groups.check_halfs()
		#	if s1 > 0.1: continue

		#	print x.name,interest,groups.top,groups.bottom,s1,s2
		#	continue
	#		print s1,s2
	#		continue

#			if not groups.marker_filter(min_obs = 0.2, pass_obs = 0.4, pass_frac = 1.25):
#				continue
			if not groups.check_dists(): 
				continue 
			hiE,lowE,hi_chi, low_chi = groups.check_pts()
			#if hi_chi > 0.15: continue

			print x.name, interest, groups.top, groups.bottom, 'enrich', hiE,lowE,'chis', hi_chi, low_chi, 'split-chis',s1,s2,"|",groups.obs[groups.top], groups.obs[groups.bottom]
			continue


		jCands = dd(list) 
		for x_idx,x in enumerate(X):

#			print x.name

			jKey,allObs,allDiffs,allLocs = {},[0],[0],[]
			svals = scale_vals([log(x.cnts[g]+1.0) for g in range(len(Y))])
			sranks = scale_vals([s.notes['rank'][x.idx] for s in Y])

			for grp,G in g_key.items():
				LC = [svals[g] for g in G]
				LR = [sranks[g] for g in G] 
				gMean =  np.mean(LC) 
				gQ    =  np.percentile(LC,90)
				gObs = len([lc for lc in LC if lc>0]) / float(len(G))
				gRank = np.mean([Y[g].notes['rank'][x.idx] for g in G])
				gRank = np.mean(LR)

				jKey[grp] = [np.array([gObs,gRank,gMean,gQ]),G,Y[G[0]].notes['color']]
				allLocs.append((np.linalg.norm(np.array([0,0,0,0]) - jKey[grp][0]),grp))
				allObs.append(gObs) 


			if max(allObs) < 0.2: continue 
			if (max(allObs) < 0.75) and (max(allObs) / (min(allObs[1::])+0.01)) < 1.5:  continue 

			print allObs,interest,x.name,g_key.keys()
			print min(allObs[1::]), max(allObs)


#			if max(allObs) < min_obs: continue 
			allLocs.sort(reverse=True) 

			maxG = allLocs[0][1]
			dists = dd(list) 
			grps = jKey.keys()
			big,small,med = 0,0,0

			for i in range(len(g_labels)-1):
				for j in range(i+1,len(g_labels)):
					gI,gJ = g_labels[i],g_labels[j] 
					d = np.linalg.norm(jKey[gI][0] - jKey[gJ][0])
					allDiffs.append(d) 
					dists[gI].append(d) 
					dists[gJ].append(d) 

			mD = dists[maxG]

			

		#	if max(allDiffs) < min_max_diff: continue
		#	if np.mean(mD) < min_mean_max_diff: continue 
		#	if np.median(mD) < min_mean_max_diff: continue 
		#	if max(mD) < min_max_diff: continue 
		

	
			jDiff = np.mean(mD) 
			jCands[maxG].append((jDiff,x,jKey.items()))


		jFound = sorted([(len(jG),jk) for (jk,jG) in jCands.items()],reverse=True)


		if len(jFound) == 0: continue 
		



		

		xLen,yLen = len(jFound),min(jFound[0][0],6)

		subplot = rage_subplots.subplot(xLen,yLen,options) 

		for jL,jName in jFound:


			jData = sorted(jCands[jName],reverse=True)[0:yLen]

			for jIdx,(jDiff,fX,jKeys) in enumerate(jData):
			
	
				for jkey in jKeys:
					cName = jkey[0]
					(cObs,cRank,cMean,cQ),cLocs,cColor = jkey[1]			

					LC = [svals[g] for g in cLocs]
					LR = [sranks[g] for g in cLocs] 
					


					for xz,yz in zip(LR,LC):
						subplot.ax.scatter(xz,yz,color=cColor,s=20) 
						



	
					subplot.ax.scatter(cRank,cQ,color=cColor,s=150) 
					subplot.ax.scatter(cObs,cMean,color=cColor,s=150,marker='s') 
					#plt.show() 

				if jIdx + 1 < len(jData):

					subplot.update({'title':fX.name})
			
			subplot.skip_row() 
		

		subplot.add_legend(g_labels,g_colors)
		#plt.subplots_adjust(left=0.03, bottom=0.08, right=0.85, top=0.90,wspace=0.25,hspace=0.3)
		subplot.save(options.prefix+'_markers_'+interest+'.pdf',{'title': interest})

 		sys.exit() 


		#plt.show() 
	sys.exit() 



			


