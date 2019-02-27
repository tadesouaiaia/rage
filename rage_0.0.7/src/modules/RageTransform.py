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
import random
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
import random
from operator import mul

from scipy.signal import savgol_filter as svgf 
from math import exp
from math import factorial 
import copy 

from Rage_IO import rage_outputs, rage_inputs
#from Rage_Plots import rage_subplots
from Rage_Plots import  rage_dimplots as dplot 

#from Rage_Transforms import rage_KDE
from Rage_Transforms import rage_DR
from Rage_Transforms import rage_dim
#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 



def scale_vals(vals,f1=0,f2=1):
        scaler = MinMaxScaler(feature_range=(f1,f2))
        return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0]


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))




class Transform:
        def __init__(self,rage):

		self.rage = rage 
		self.options = rage.args
		self.attributes = rage.args.color+[x for x in [rage.args.marker,rage.args.size] if x != None]
		 
	



	def run(self):
		R = self.rage
		R.progress.start_major('Rage Transforms') 
		
		if R.args.command == 'condensation': 	self.condense()  
		elif R.args.command == 'fabrication':	self.fabricate() 

		
		else: 
			if len(self.attributes) > 0: R.data.filter_samples_by_attributes(covariates=self.attributes).scale_and_transform() 

			R.data.samples.create_plot_labels(self.options) 

			if R.args.command == 'dim': 	self.run_dim(R) 
			 





			elif R.args.condensedcnts:	
				condense_data = self.rage.condensed_data
	#			if len(R.args.color) > 0:	R.condensed_data.filter_samples_by_attributes(self.options.color).scale_and_transform() 

				R.condensed_data.samples.create_plot_labels(self.options,SIZE=50) 

				if R.args.command == 'pca': self.run_condense_pca(R) # .condense_data) 
		

			else: 

				if R.args.command == 'pca' and self.options.iterative: self.run_iterative_pca(R) 
				elif R.args.command == 'pca':                          self.run_pca(R) 

				elif R.args.command == 'tsne': 				self.run_tsne(R)  





				else:
					print 'wtf', R.args.command




	def write_coeffs(self,coefs,F,out): 

		F_key = dd(list) 
		for i,C in enumerate(coefs):
			for j,(vs,vl,vi) in enumerate(C): F_key[F[vi].name].append((j,vl))


		w = open(out,'w') 
		w.write("%-50s %10s %20s %10s %20s %10s %20s\n" % ('---','R1','V1','R2','V2','R3','V3'))
		for k,C in F_key.items(): 
			w.write("%-50s %10d %20f %10d %20f %10d %20f\n" % (k,C[0][0],C[0][1],C[1][0],C[1][1],C[2][0],C[2][1]))

	def run_condense_pca(self,R): 
	#	R = self.rage
		D = R.data


		out_name = R.args.prefix+'_condensepca_'
		

		shared_features = [f.name for f in R.condensed_data.features if f.name in [fm.name for fm in D.features]]

		cF = [f for f in R.condensed_data.features if f.name in shared_features] 
		rF = [f for f in D.features if f.name in shared_features] 


		Yc = [[s.cnts[f.idx] for s in R.condensed_data.samples] for f in cF]
		Y = [[s.cnts[f.idx] for s in D.samples] for f in rF]

		dr = rage_DR.DR(R.args,R.progress)
		out_name = R.args.prefix+'_transformsamples_'
		R.progress.start_minor('PCA') 
		r_members = R.data.samples


		condense_run = dr.set_y_matrix(Yc, LOG_TRANSFORM=True).pca(req='FULL') 
		self.write_coeffs(condense_run['coefs'],cF,out_name+'coefs.out') 


		pca_run = dr.set_y_matrix(Y, LOG_TRANSFORM=True, SET_TRANSFORM=True).pca(req='FULL') 

		dimplot = dplot.dimplot(1,1,R.args,R.progress)
		
		shared_run = {'pts': condense_run['pts']+pca_run['pts'], 'axes': condense_run['axes']}
		shared_samples = [s for s in R.condensed_data.samples] + [s for s in R.data.samples] 


		dimplot.add_data(shared_samples,[shared_run],{'title': 'SHARED_PCA', 'out': out_name+'pca.pdf'}) 



#		dimplot.add_data(R.condensed_data.samples,[condense_run],{'title':'PCA','out': out_name+'pca.pdf'})
#		dimplot.add_data(R.data.samples,[pca_run],{'title':'PCA','out': out_name+'pca.pdf'})
		dimplot.finish(out_name+'pca') 
		R.progress.end() 



		sys.exit() 



	def precomp_pca(self,R,coeffs,PLEN=1500,MAX_COEFS=8):


		coeff_key = dd(lambda: {}) 
		scale_key = dd(lambda: {}) 
		projection_key = dd(lambda: {}) 
		for line in coeffs: 
			line = line.split() 
			if line[0] == '---': continue

			for i in range(2,len(line),2):

				coeff_key[(i/2)-1][line[0]] = float(line[i])
				if i >= 40: break 
			


		for f in self.D.features: 
			if f.name not in coeff_key[0]: 
				for i in coeff_key.keys(): coeff_key[i][f.name] = 0.0 
		

		for s in self.S: projection_key[s] = sorted([[c,self.D.features[i].name,[coeff_key[n][self.D.features[i].name] for n in range(len(coeff_key))]] for i,c in s.cnts.items()],reverse=True) 
				

		prj_len =  max(min([len(X) for X in projection_key.values()]),PLEN) 
		self.plt_name+='_projected_'+str(prj_len)
		pca_key = dd(list) 
		tsne_key = {} 

		for s in self.S:

			LK_DATA = [[pk[2][n]*log(pk[0],2) for pk in projection_key[s]] for n in range(len(coeff_key.keys()))]
			RK_DATA = [[pk[2][n]*log(pk[0],2) for pk in projection_key[s]] for n in range(len(coeff_key.keys()))]


			RK_DOT = [sum(rk) for rk in RK_DATA]
			LK_DOT = [sum(lk) for lk in LK_DATA]
			RK_PROJ = [sum(rk[0:prj_len]) for rk in RK_DATA]
			LK_PROJ = [sum(lk[0:prj_len]) for lk in LK_DATA]

			pca_key['LOGDOT'].append(LK_DOT)
			pca_key['RAWDOT'].append(RK_DOT)

			pca_key['LOGPRJ'].append(LK_PROJ)
			pca_key['RAWPRJ'].append(RK_PROJ)



		rawdot = {'pts': pca_key['RAWDOT'], 'axes': ['PC'+str(x+1)+'-RAWDOT' for x in range(len(coeff_key.keys()))]}
		logdot = {'pts': pca_key['LOGDOT'], 'axes': ['PC'+str(x+1)+'-LOGDOT' for x in range(len(coeff_key.keys()))]}
		rawprj = {'pts': pca_key['RAWPRJ'], 'axes': ['PC'+str(x+1)+'-RAWPRJ' for x in range(len(coeff_key.keys()))]}
		logprj = {'pts': pca_key['LOGPRJ'], 'axes': ['PC'+str(x+1)+'-LOGPRJ' for x in range(len(coeff_key.keys()))]}



		for kp,kpts in pca_key.items(): 
			w_name = self.plt_name+'_'+kp+'_pca_proj.pts'
			w= open(w_name,'w') 
			w.write("%-50s %10s %10s %10ss %10s %10s \n" % ('---','PC1','PC2','PC3','PC4','PC5'))
			for si,p in enumerate(kpts):
				s = self.S[si] 
				w.write("%-50s %10f %10f %10f %10f %10f\n" % (s.name,p[0],p[1],p[2],p[3],p[4]))
			w.close() 


		for dc in [(0,1),(2,3)]: 
			p_name = self.plt_name+'_'+'-'.join([str(ss) for ss in dc])+'_'
			dimplot = dplot.dimplot(R.args,R.progress).add_data(self.S,[rawdot,logdot,rawprj,logprj],dim_comps=[dc,dc,dc,dc],NAMES=False).finish(p_name+'pca')
			dimplot = dplot.dimplot(R.args,R.progress).add_data(self.S,[rawdot,logdot,rawprj,logprj],dim_comps=[dc,dc,dc,dc],NAMEOUTLIERS=True).finish(p_name+'exnamed_pca')

		dr = rage_DR.DR(R.args,R.progress)
		tsne_key = {k: dr.tsne(pca_pts=vals,axes_prefix=k) for k,vals in pca_key.items()}
		t_runs = [tsne_key['RAWDOT'],tsne_key['LOGDOT'],tsne_key['RAWPRJ'],tsne_key['LOGPRJ']]

		for kp,kpts in tsne_key.items(): 
			w_name = self.plt_name+'_'+kp+'_tsne_proj.pts'	
			w= open(w_name,'w') 
			w.write("%-50s %10s %10s\n" % ('---','TSNE1','TSNE2'))
			for si,p in enumerate(kpts['pts']):
				s = self.S[si] 
				w.write("%-50s %10f %10f\n" % (s.name,p[0],p[1]))
			w.close() 



		dimplot = dplot.dimplot(R.args,R.progress).add_data(self.S,t_runs,dim_comps=[(0,1),(0,1),(0,1),(0,1)],NAMES=False).finish(self.plt_name+'_tsne')
		dimplot = dplot.dimplot(R.args,R.progress).add_data(self.S,t_runs,dim_comps=[(0,1),(0,1),(0,1),(0,1)],NAMEOUTLIERS=True).finish(self.plt_name+'_exnamed_tsne')
		R.progress.end() 
		sys.exit() 




	def run_dim(self,R): 

		

		dim = rage_dim.DIM(R)
		if 'test' in R.args.scaling: 
			dim.run_scaling_demo(TEST=True) 
		else:
			dim.run_scaling_demo(INFO=R.args.scaling[0]) 
			

		



	def run_pca(self,R):
		self.LT, self.SC, self.STD, ttype, cstr = False, False, False, '_', '_'

		if len(self.options.color)>0: 	  cstr+= '-'.join(self.options.color) 
		if self.options.marker != None:   cstr+= '-'+self.options.marker		

		if 'log' in self.options.notes:     
			self.LT= True 
			ttype += 'LOG_'
		if 'scale' in self.options.notes:   
			self.SC = True
			ttype += 'SCALE_'
		elif 'std'   in self.options.notes: 
			self.STD = True 
			ttype += 'STD_'

		if len(ttype) < 3: ttype = '_RAW'

		self.out_name = R.args.prefix+'_'+cstr+'_transformsamples'+ttype
		self.plt_name = R.args.prefix+'_'+cstr+'_transformsamples'+ttype 
		self.D = R.data 
		self.S = self.D.samples 
		self.Y = [[s.cnts[f.idx] for s in self.D.samples] for f in self.D.features]

		if self.options.coeffs:
			self.precomp_pca(R,self.options.coeffs) 
			return

		else: 


			dr = rage_DR.DR(R.args,R.progress)
			R.progress.start_minor('PCA') 
			dr.set_y_matrix(self.Y, LOG_TRANSFORM=self.LT,SCALE=self.SC,STD_SCALE=self.STD)


			pca_run = dr.pca(req='FULL') 

			F_key = dd(list) 
			for i,C in enumerate(pca_run['coefs']):
				for j,(vs,vl,vi) in enumerate(C): 
					F_key[self.D.features[vi].name].append((j,vl))

			w= open(self.out_name+'pca_coefs.out','w') 
			w.write("%-50s %5s %10s %5s %10s %5s %10s\n" % ('---','R1','V1','R2','V2','R3','V3'))
			for k,C in F_key.items(): 
				w.write("%-50s" % (k))
				for i in range(len(C)): 	w.write(" %5d %10f" % (C[i][0],C[i][1]))
				w.write('\n') 
			w.close() 

			w= open(self.out_name+'pca_pts.out','w') 
			w.write("%-50s %10s %10s %10ss %10s %10s \n" % ('---','PC1','PC2','PC3','PC4','PC5'))
			for p,s in zip(pca_run['pts'],R.data.samples):	w.write("%-50s %10f %10f %10f %10f %10f\n" % (s.name,p[0],p[1],p[2],p[3],p[4]))
			w.close() 

				
			dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples,[pca_run],dim_comps=[(0,1),(2,3),(4,5),(6,7)],NAMES=False).finish(self.plt_name+'_sample_pca')
			dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples,[pca_run],dim_comps=[(8,9),(10,11),(12,13),(14,15)],NAMES=False).finish(self.plt_name+'_sample_hipca')
				
			R.progress.end() 
			
			
			tsne_run = dr.tsne() 	
			w= open(self.out_name+'tsne_pts.out','w') 
			w.write("%-50s %10s %10s \n" % ('---','TS1','TS2'))
			for p,s in zip(tsne_run['pts'],R.data.samples):	w.write("%-50s %10f %10f \n" % (s.name,p[0],p[1]))
			w.close() 

			dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples,[tsne_run],dim_comps=[(0,1)],NAMES=False).finish(self.plt_name+'_sample_tsne')

			R.progress.end() 






































	def run_tsne(self,R): 

		dr = rage_DR.DR(R.args,R.progress)
		out_name = R.args.prefix+'_transformsamples_TSNE'
		R.progress.start_minor('TSNE') 
		r_members = R.data.samples
		r_matrix = R.data.matrix('log')

		tsne_run = dr.run_tsne(r_matrix)
		dimplot = dplot.dimplot(2,2,R.args,R.progress).add_data(r_members,[tsne_run],{'title':'TSNE','out': out_name+'tsne.pdf'})

		dimplot.finish(out_name+'tsne.pdf') 
		R.progress.end() 




















  
	def condense(self):
		R = self.rage
		condenser = R.args.predictors[0]


		self.D 	     = R.data.filter_samples_by_attributes(R.args.predictors,keep=True)
		self.seg,self.seg_lens      =         self.D.samples.segregate(condenser)
		grps = {} 




		for k,grp in self.seg.items():
			random.shuffle(grp)
			for i,g in enumerate([grp[i:i+R.args.condense] for i in range(0,len(grp),R.args.condense)]):	grps[str(i+1)+'~'+k] = g



		group_keys = grps.keys() 
		group_cnts = []
		for f in R.data.features:

			f_cnts =  [f.cnts[i] for i in range(len(self.D.samples))]
			#f_cnts =  [log(f.cnts[i]+1.0,2) for i in range(len(self.D.samples))]

			group_cnts.append([sum([f_cnts[i] for i in grps[k]]) for k in group_keys])
			
		rage_outputs.count_file(R.args).write_row_col_data(R.data.features,group_keys,group_cnts,{'suffix': "_".join(R.args.predictors)+'_Econdense_'+str(R.args.condense)+'.cnts'})
		R.progress.end() 

	def fabricate(self):
		R = self.rage
		fabricator = R.args.predictors[0] 

		self.D 	     = R.data.filter_samples_by_attributes(R.args.predictors,keep=True)
		self.seg,self.seg_lens      =         self.D.samples.segregate(fabricator)
		grps={} 
		trials = 10
		for k,grp in self.seg.items():
			for n in range(trials): 
				random.shuffle(grp)
				for i,g in enumerate([grp[i:i+int(R.args.size)] for i in range(0,len(grp),int(R.args.size))]):	grps[str((n+1)*(i+1))+'~'+k] = g


		group_keys = grps.keys() 
		group_cnts = []

		for f in R.data.features:

			f_cnts =  [f.cnts[i] for i in range(len(self.D.samples))]
			#f_cnts =  [log(f.cnts[i]+1.0,2) for i in range(len(self.D.samples))]

			group_cnts.append([sum([f_cnts[i] for i in grps[k]]) for k in group_keys])
			
		rage_outputs.count_file(R.args).write_row_col_data(R.data.features,group_keys,group_cnts,{'suffix': "_".join(R.args.predictors)+'_fabricate'+str(R.args.size)})
		R.progress.end() 





	def write_pca_data(self,pca_run, cols, rows, out_name, SIZE = 10, PROJECTION_SIZE=2 ): 

		sample_key = dd(list) 

		coeff_key = dd(list) 
		for i,C in enumerate(pca_run['coefs']):
			for j,(vs,vl,vi) in enumerate(C): coeff_key[rows[vi].name].append((j,vl))

		w= open(out_name+'pca_coefs.out','w') 
		w.write("%-50s %5s %10s %5s %10s %5s %10s\n" % ('---','R1','V1','R2','V2','R3','V3'))
		for k,C in coeff_key.items(): 
			w.write("%-50s" % (k))
			for i in range(len(C)): w.write(" %5d %10f" % (C[i][0],C[i][1]))
			w.write('\n') 
		w.close() 




		P_LEN = max(100,min([len(s.cnts) for s in cols]))
		s_project = [[0 for i in range(SIZE)] for s in cols]
		s_rank = [[0 for i in range(SIZE)] for s in cols]
		for i in range(SIZE): 
			coeff_order = sorted([(coeff_key[r.name][i],r.name,r.idx) for r in rows])
			p_key = {r.idx: coeff_key[r.name][i] for r in rows} 
			for j,s in enumerate(cols):
				sc_srt = sorted([(s.cnts[c],c) for c in s.cnts],reverse=True)[0:P_LEN]  
				s_rank[j][i] = sum([p_key[cIdx][1] * P_LEN-k for k,(cVal, cIdx) in enumerate(sc_srt)])
				s_project[j][i] = sum([p_key[cIdx][1] * cVal for cVal, cIdx in sc_srt])

		s_projection = {s.name: [s_rank[j][i] for i in range(PROJECTION_SIZE)] for j,s in enumerate(cols)}
		#s_projection = {s.name: [s_project[j][i] for i in range(PROJECTION_SIZE)] for j,s in enumerate(cols)}
		#is_projection = {s.name: [s_project[j][i] for i in range(PROJECTION_SIZE)] for j,s in enumerate(cols)}
#		print "THIS RANK" 
#		s_projection = {s.name: [s_project[j][i] for i in range(PROJECTION_SIZE)] for j,s in enumerate(cols)}

		w= open(out_name+'pca_project.out','w') 
		w.write("%-50s %10s %10s %10ss %10s %10s \n" % ('---','PR1','PR2','PR3','PR4','PR5'))
		for p,s in zip(s_project,cols):	
			w.write("%-50s %50s\n" % (s.name," ".join([str(ps) for ps in p])))
		w.close() 

		w= open(out_name+'pca_rankproject.out','w') 
		w.write("%-50s %10s %10s %10ss %10s %10s \n" % ('---','PR1','PR2','PR3','PR4','PR5'))
		for p,s in zip(s_rank,cols):	
			w.write("%-50s %50s\n" % (s.name," ".join([str(ps) for ps in p])))
		w.close() 


		w= open(out_name+'pca_pts.out','w') 
		w.write("%-50s %10s %10s %10ss %10s %10s \n" % ('---','PC1','PC2','PC3','PC4','PC5'))
		for p,s in zip(pca_run['pts'],cols):	
			#w.write("%-50s %10f %10f %10f %10f %10f\n" % (s.name,p[0],p[1],p[2],p[3],p[4]))
			w.write("%-50s %50s\n" % (s.name," ".join([str(ps) for ps in p[0:50]])))
		w.close() 


		w= open(out_name+'vars_exp.out','w') 
		w.write("%-30s %10s \n" % ('---', 'var'))
		for v in pca_run['axes']: 	w.write("%-30s %10s \n" % (v.split()[0],v.split()[1]))
		w.close() 
		
		return coeff_key, s_projection


	def run_iterative_pca(self,R,GROUPS=2):
		self.R = R
		D, F, S = R.data , R.data.features, R.data.samples 
		Y = [[s.cnts[f.idx] for s in D.samples] for f in D.features]


		
		self.D = R.data 

		### PARAMS ### 


		GROUPS = 2 
		SELECTION_FRACTION = 0.24	


		dr = rage_DR.DR(R.args,R.progress)
		out_name = R.args.prefix+'_transformsamples_'
		R.progress.start_minor('PCA') 

		LT = ('log' in self.options.notes)
		SC = ('scale' in self.options.notes)
		STD = ('std' in self.options.notes)

	

		if GROUPS == 200: 


			iterplot = dplot.iterplot(self.R.data.samples,self.R.data.features,R.args,R.progress,6)
			#sc = MinMaxScaler()
                	#Y = [[x[0] for x in sc.fit_transform(np.array(y).reshape(-1,1))] for y in Y]
#			dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples,[pca_run],dim_comps=[(0,1),(2,3),(4,5),(6,7)],NAMES=False).finish(R.args.prefix+'_transformsamples_sample_pca')

			PCA_1, Y_1, S_1 = self.iterative_run(Y,[s.name for s in S],iteration=1,SELECTION_FRACTION=0.02,MAX_SIZE=2,DIMS=1) #LOWEST=True) 
			for s in [s for s in S_1 if len(s.split('@')) > 1]:	print 1,s  


			iterplot.add_data(PCA_1, [s.name for s in self.R.data.samples]) 
			


			PCA_2, Y_2, S_2 = self.iterative_run(Y_1,S_1,iteration=2,SELECTION_FRACTION=0.05,MAX_SIZE=2,DIMS=1,LOWER=True) 
			for s in [s for s  in S_2 if len(s.split('@')) > 1]:	print 2,s  
 
			iterplot.add_data(PCA_2, S_1) 

			PCA_3, Y_3, S_3 = self.iterative_run(Y_2,S_2,iteration=2,SELECTION_FRACTION=0.025,MAX_SIZE=3,DIMS=1,LOWEST=True) 
			for s in [s for s  in S_3 if len(s.split('@')) > 1]:	print 3,s  
 
			iterplot.add_data(PCA_3, S_2) 
			
			PCA_4, Y_4, S_4 = self.iterative_run(Y_3,S_3,iteration=2,SELECTION_FRACTION=0.50,MAX_SIZE=2,DIMS=1,TAIL=False) 
			for s in [s for s  in S_4 if len(s.split('@')) > 1]:	print 4,s  
			iterplot.add_data(PCA_4, S_3) 

#			PCA_5, Y_5, S_5 = self.iterative_run(Y_4,S_4,iteration=2,SELECTION_FRACTION=0.12,MAX_SIZE=3,DIMS=1,TAIL=True,LOWER=True) 
#			for s in [s for s  in S_5 if len(s.split('@')) > 1]:	print 5,s  
#			iterplot.add_data(PCA_5, S_4) 

			
#			PCA_6, Y_6, S_6 = self.iterative_run(Y_5,S_5,iteration=2,SELECTION_FRACTION=0.40,MAX_SIZE=4,DIMS=2,TAIL=False) 
#			for s in [s for s  in S_6 if len(s.split('@')) > 1]:	print 6,s  
#			iterplot.add_data(PCA_6, S_5) 
			

			plt.savefig('PCA_ITERATIVE.pdf') 


				
			
						#print cn, pn, np.linalg.norm(c_pts-p_pts)


			
		#dimplot = dplot.dimplot(1,1,R.args,R.progress).add_data(R.data.samples,[pca_run],dim_comps=[(0,1),(2,3)],NAMES=True).finish(out_name+'named.pca')
		#dimplot = dplot.dimplot(1,1,R.args,R.progress).add_data(R.data.samples,[pca_run],dim_comps=[(0,1),(1,2),(3,4),(5,6)],NAMES=True).finish(out_name+'named.pca')
		#dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples,[pca_tran],dim_comps=[(0,1),(2,3),(4,5),(6,7)],NAMES=False).finish(out_name+'Cpca')

	
		dr.set_y_matrix(Y, LOG_TRANSFORM=LT,SCALE=SC,STD_SCALE=STD) 
		pca_run = dr.pca(req='FULL') 


		self.coeff_key, sample_data  = self.write_pca_data(pca_run, S, F, R.args.prefix+'_sampletransform_')		
		dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples, [pca_run],dim_comps=[(0,1),(2,3),(4,5),(6,7)],NAMES=False).finish(R.args.prefix+'_transformsamples_sample_pca')
		s_pts = [[p[1] for p in self.coeff_key[f.name]] for f in R.data.features]
		pca_run['pts'] = s_pts 
		dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.features,[pca_run],dim_comps=[(0,1),(2,3),(4,5),(6,7)], NAMES=True).finish(R.args.prefix+'_transformsamples_feature_pca')




		dr.set_y_matrix(Y, TRANSPOSE = True, SCALE=True) 
		pca_tran = dr.pca(req='FULL') 

		self.s_key, feature_data = self.write_pca_data(pca_tran,  F,S, R.args.prefix+'_featuretransform_')
		dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.features,[pca_tran],dim_comps=[(0,1),(2,3),(4,5),(6,7)], NAMES=True).finish(R.args.prefix+'_transformfeatures_feature_pca')

		s_pts = [[p[1] for p in self.s_key[s.name]] for s in R.data.samples]
		pca_tran['pts'] = s_pts 
		dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples,[pca_tran],dim_comps=[(0,1),(2,3),(4,5),(6,7)], NAMES=False).finish(R.args.prefix+'_transformfeatures_samples_pca')


		
#		print sample_data 


			
		#dimplot = dplot.dimplot(1,1,R.args,R.progress).add_data(R.data.samples,[pca_run],dim_comps=[(0,1),(2,3)],NAMES=True).finish(out_name+'named.pca')
		#dimplot = dplot.dimplot(1,1,R.args,R.progress).add_data(R.data.samples,[pca_run],dim_comps=[(0,1),(1,2),(3,4),(5,6)],NAMES=True).finish(out_name+'named.pca')
		#dimplot = dplot.dimplot(R.args,R.progress).add_data(R.data.samples,[pca_tran],dim_comps=[(0,1),(2,3),(4,5),(6,7)],NAMES=False).finish(out_name+'Cpca')





	def iterative_run(self,Y,S_NAMES,iteration=1,GROUPS=2,SELECTION_FRACTION=0.24,MAX_SIZE=4,DIMS=3,TAIL=False,LOWER=False,LOWEST=False): 


		LT = ('log' in self.options.notes)
		SC = ('scale' in self.options.notes)
		STD = ('std' in self.options.notes)
		DIMS=2
	
		if GROUPS == 2: 


			dr = rage_DR.DR(self.R.args,self.R.progress)
			#dr.set_y_matrix(Y, LOG_TRANSFORM=LT,SCALE=SC,STD_SCALE=STD) 
			dr.set_y_matrix(Y, CENTER=False,SCALE=True) 
			pca_run = dr.pca(req='FULL') 
			pca_pts = pca_run['pts'] 
			coeff_key = dd(list) 
			for i,C in enumerate(pca_run['coefs']):
				for j,(vs,vl,vi) in enumerate(C): coeff_key[self.R.data.features[vi].name].append((j,vl))



			S_CNTS = [sorted([(Y[i][j],self.R.data.features[i].name) for i in range(len(Y)) if Y[i][j] > 0],reverse=True) for j in range(len(S_NAMES))]


			S_OBS =  {S_NAMES[j]: len(S_CNTS[j]) for j in range(len(S_NAMES))}
			S_MED =  np.median(S_OBS.values())
			S_25 =  np.percentile(S_OBS.values(),25)

			P_LEN  =  max(min([len(sc) for sc in S_CNTS]),100)
			S_PROJECTIONS = dd(lambda: dd(list)) 
			for i,C in enumerate(S_CNTS): 
				for v,f in C[0:P_LEN]: 
					for j in range(len(coeff_key[f])):
						S_PROJECTIONS[S_NAMES[i]][j].append(v*coeff_key[f][j][1])
			
			S_PROJECTIONS =  {s_name: [np.mean(S_PROJECTIONS[s_name][j]) for j in range(DIMS)] for s_name in S_PROJECTIONS.keys()}
			
			dr = rage_DR.DR(self.R.args,self.R.progress)
			#dr.set_y_matrix(Y, TRANSPOSE = True, SCALE=True) 
			dr.set_y_matrix(Y, TRANSPOSE=True,CENTER=False,SCALE=True) 
			pca_tran = dr.pca(req='FULL') 
			s_key = dd(list) 
			for i,C in enumerate(pca_tran['coefs']): 
				for j,(vs,vl,vi) in enumerate(C): s_key[S_NAMES[vi]].append((j,vl))

			print DIMS

#			print '--- name2 idx1 idx2 obs1 obs2 | aCoef bCoef'
			print '--- idx1 obs1 | aCoef aPCA aProj',
			print '| idx2 obs2 | bCoef bPCA bProj'
			for i in range(len(S_NAMES)):
				aName = S_NAMES[i] 
				aCoefs,aPCA, aProj,aObs = [x[1] for x in s_key[aName]][0:DIMS], pca_pts[i][0:DIMS], S_PROJECTIONS[aName], S_OBS[aName]


				for j in range(i+1,len(S_NAMES)): 

					print aName,i,aObs,'|',	
					print ",".join([str(round(x,3)) for x in aCoefs]),
					print ",".join([str(round(x,3)) for x in aPCA]),
					print ",".join([str(round(x,3)) for x in aProj]),'|',


					bName = S_NAMES[j] 
					bCoefs,bPCA, bProj,bObs = [x[1] for x in s_key[bName]][0:DIMS], pca_pts[j][0:DIMS], S_PROJECTIONS[bName], S_OBS[bName]

					print bName,j,bObs,'|',	
					print ",".join([str(round(x,3)) for x in bCoefs]),
					print ",".join([str(round(x,3)) for x in bPCA]),
					print ",".join([str(round(x,3)) for x in bProj]),

					coefD = np.linalg.norm(np.array(aCoefs) - np.array(bCoefs))
					coefP = np.linalg.norm(np.array(aPCA) - np.array(bPCA))
					coefPr = np.linalg.norm(np.array(aProj) - np.array(bProj))
				

					print '|', coefD, coefP,coefPr








			sys.exit() 



			if TAIL: s_order =  sorted([(s_key[s][0][1],s) for s in s_key.keys()])
			else: 	 s_order =  sorted([(s_key[s][0][1],s) for s in s_key.keys()],reverse=True)


			PAIRS, S_STOP, FOUND,SCAN, ADDED = [], len(s_order) * SELECTION_FRACTION, dd(bool) , 10, 0 
			for i in range(len(s_order)): 
				i_val, i_name = s_order[i]
				if LOWER and S_OBS[i_name] > S_MED: continue 
				elif LOWEST and S_OBS[i_name] > S_25: continue 
			 
				i_size, i_loc,i_dists, j  = len(i_name.split('@')), np.array(S_PROJECTIONS[i_name] ),[], i + 1
				if FOUND[i] or i_size >= MAX_SIZE: continue 
				while j < len(s_order): 
					if not FOUND[j] and (i_size+len(s_order[j][1].split('@'))) <= MAX_SIZE: 
						j_val, j_name = s_order[j]
						i_dists.append((np.linalg.norm(i_loc - np.array(S_PROJECTIONS[j_name])),j_name,j))
					if len(i_dists) >= SCAN: break  
					j+=1
				if len(i_dists) > 0:
					dist,j_name,j_idx = sorted(i_dists)[0] 
					PAIRS.append((i_name,j_name))
					FOUND[j_idx], FOUND[i_name], FOUND[j_name] = True, True, True 
					
					ADDED +=1
				if ADDED > S_STOP: break 

	


			NEW_Y, NEW_IDX, IDX_KEY = [], {}, {s: i for i,s in enumerate(S_NAMES)} 
			for a,b in PAIRS: NEW_IDX[a+"@"+b] = [IDX_KEY[a],IDX_KEY[b]]
			for s,i in [(s,i) for (s,i) in IDX_KEY.items() if s not in FOUND]: NEW_IDX[s] = [i] 

			NEW_SAMPLE_NAMES = NEW_IDX.keys() 
			#for y in Y: NEW_Y.append([max([y[j] for j in NEW_IDX[ns]]) for ns in NEW_SAMPLE_NAMES]) 
			#for y in Y: NEW_Y.append([sum([y[j] for j in NEW_IDX[ns]]) for ns in NEW_SAMPLE_NAMES]) 
			for y in Y: NEW_Y.append([np.mean([y[j] for j in NEW_IDX[ns]]) for ns in NEW_SAMPLE_NAMES]) 
	

			return pca_run, NEW_Y, NEW_SAMPLE_NAMES


			





	def select_iterative_pairs(self,pca_run,ROUND): 
		self.V_COMPS = 1	
		self.STD_THRESH = 2.0
		self.M_THRESH = 2.0
		self.D_THRES = 1.1

	
		for k,C in self.coeff_key.items():
			print k 
			for i in range(len(C)): print C[i] 
				

		sys.exit() 	

		cand_key, p_key = dd(list), dd(list) 
		for i,P in enumerate(pca_run['pts']):
			for j in range(len(P)):
				if j == self.V_COMPS: break 
				p_key[j].append((P[j],i))


		CAND_SIZE = len(pca_run['pts']) / 4
		if CAND_SIZE < 10: CAND_SIZE = 10 
		for j  in range(len(p_key)): 
			pp = sorted(p_key[j]) 
			pn = [p[0] for p in pp]
			pJumps = [pn[i] - pn[i-1] for i in range(1,len(pn))]
			jMean, jMed = np.mean(pJumps), np.median(pJumps) 
			pMean, pMed, pStd = np.mean(pn), np.median(pn) , np.std(pn) * self.STD_THRESH
			pSmall = pp[0:CAND_SIZE] 
			pLarge = pp[-1::-1][0:CAND_SIZE]
			
			for pi,p in enumerate(pSmall):				 
				if pMean - p[0] < pStd: break
			pSmall = pSmall[0:pi]					
			for pi,p in enumerate(pLarge):		 
				if  p[0] - pMean < pStd: break
			pLarge = pLarge[0:pi]


			print j, len(pSmall), len(pLarge) 
			pSmall = [] 
			smallPass = self.search_select(pSmall,jMean,jMed)
			bigPass   = self.search_select(pLarge,jMean,jMed)
		

			for (xA,xB) in smallPass + bigPass:
				cand_key[xA].append(xB) 

			if j > self.V_COMPS:
				break 


		for c,partners in cand_key.items(): 
			cn = self.samples[c].name 
			for p in partners: 
				pn = self.samples[p].name
				c_pts = pca_run['pts'][c][0:self.V_COMPS+1]  
				p_pts = pca_run['pts'][p][0:self.V_COMPS+1]
				print cn, pn, np.linalg.norm(c_pts-p_pts)
		return 


	def search_select(self,p_list,jMean,jMed): 

		my_pairs = [] 
		i=1
		pJumps = [999]+[fabs(p_list[k][0] - p_list[k-1][0]) for k in range(1,len(p_list))]+[999]
		gMean, gMed = np.mean(pJumps[1:-1]), np.median(pJumps[1:-1]) 
		gSep = min([np.mean(pJumps[1:-1]), np.median(pJumps[1:-1])])
		while i +1 < len(p_list): 
			a,b,c = p_list[i-1],p_list[i],p_list[i+1] 
			pA,pB,pC = a[0],b[0],c[0] 
			nA,nB,nC = self.samples[a[1]].name, self.samples[b[1]].name, self.samples[c[1]].name

			jA,jB,jC,jD = pJumps[i-1],pJumps[i],pJumps[i+1],pJumps[i+2] 
			j1,j2 =    fabs(pB-pA),fabs(pB-pC)

			if jB * self.M_THRESH < jA and jB * self.M_THRESH < jC and jB < gSep / self.D_THRES:
				my_pairs.append(sorted((a[1],b[1]))) 
				i+=2
			elif jC * self.M_THRESH < jD and jC * self.M_THRESH < jB and jC < gSep / self.D_THRES:
				my_pairs.append(sorted((b[1],c[1]))) 
				i+=2
			else:
				i+=1
		return my_pairs













	
