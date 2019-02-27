#!/usr/bin/env python


import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
from scipy.stats import variation as coVar 
import numpy as np
import random
import seaborn
from math import log
import math
from random import shuffle

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Rage_IO import rage_regression_outputs
from Rage_Plots import rage_regression_plots
#from Rage_Transforms import rage_KDE
from Rage_Transforms import rage_DR
from Rage_Regression import rage_regmodels as rt 
#from Rage_Filters    import rage_filters 
import statsmodels.stats.multitest as mlt 
#import statsmodels.sandbox.stats.multicomp as mpt 
import copy 


#from Rage_Summary import summary_hists # rage_summarize_features, rage_summarize_samples, rage_summarize_dists 
import warnings
warnings.filterwarnings("ignore")


def regression_error(msg):
        sys.stderr.write('RageRegressionError: '+msg+'\n')
        sys.exit()









class Regression:
        def __init__(self,rage):

		self.rage, self.progress, self.options  = rage, rage.progress, rage.args 
		if len([b for b in rage.args.predictors if b in rage.args.covariates])>0: regression_error('Predictors and Covariates must be distinct (54)') 

	def run(self): 

		if len(self.options.predictors) == 0: regression_error('No predictors')


		self.key_nick = self.options.sampleKey.name.split('/')[-1].split('.key')[0]
		self.dists = [] 
		for dist in self.options.dist: 
			if dist.upper() not in ['OLS-LOG','OLS','ZGP','ZNB','ZPO','NB','GP','PO']:
				regression_error('Unsupported Distribution ( '+dist+' )\n'+ 
						  'Choose From: \n'+
						  '\t    OLS         = Ordinary Least Squares\n'+
						  '\t    OLS-LOG     = Ordinary Least Squares (Log Transformed)\n'+
						  '\t    PO          = Poisson Regression\n'+
						  '\t    NB          = Negative Binomial Regression\n'+
						  '\t    GP          = Generalized Poisson Regressioon\n'+
						  '\t    ZPO         = Zero Inflated Poisson Regression\n'+
						  '\t    ZNB         = Zero Inflated Negative Binomial Regression\n'+
						  '\t    ZGP         = Zero Inflated Generalized Poisson Regressioon\n')
			else:
				self.dists.append(dist.upper()) 





		if self.options.command == 'eval-predictors': 

			self.progress.start_major('Running Regression Predictor Analysis') 
			self.eval_predictors() 
		else:
			self.D = self.rage.data.filter_samples_by_attributes(self.options.predictors,self.options.covariates).scale_and_transform(LOG=False)	
			self.V = self.D.set_sample_variables(self.options.predictors, self.options.covariates) 


			self.X, self.Xp,self.Xc = self.V.select_variables(self.V.variables), self.V.select_variables(self.V.predictors), self.V.select_variables(self.V.covariates) 	
			if self.options.zero_prob != None and self.options.zero_prob[0:3].upper() == 'OBS': self.X.zp =  np.array([[( 1 - (len(s.cnts)  / float(len(self.D.features))) )] for s in self.D.samples])
			self.feature_names, self.Y = [],[] 
			for f in self.D.features: 
				y,y_obs,y_name = [s.cnts[f.idx] for s in self.D.samples], len(f.cnts)/ float(len(self.D.samples)) , f.name
				if y_obs > 0.025 and len(set(y)) > 5: 
					self.Y.append(y) 
					self.feature_names.append(y_name) 

			self.D.samples.create_plot_labels(self.options)


			if self.options.command == 'dex':
				self.progress.start_major('Running Dex Evaluation') 
				self.evaluate_dex() 

			elif self.options.command == 'eval-model':
				self.progress.start_major('Running Model Evaluation') 
				self.evaluate_model() 	


			elif self.options.command == 'eval-covariates':
				self.evaluate_model2() 	
			








	def evaluate_dex(self,FDR_CUTOFF=0.05):
		
		pred_format_str,pred_format_num  =  '%-50s %30s %20s %15s %10s %10s %10s %10s %10s' ,'%-50s %30s %20s %15s %10s %10.2f %10.2f %10.2e %10.2e' 
		vals_format_str,vals_format_num  =  '%-40s %15s %10s %20s %25s %30s %40s %15s %15s %15s\n' , '%-40s %15s %10.1e %20s %25s %30s %40s %15.3f %15.1e %15.1e \n'

		for dist in self.dists: 
			pv_summary = dd(lambda: dd(list))  

			if len(self.options.covariates) == 0: dex_prefix = self.options.prefix+'-'+self.key_nick+'-'+'-'.join(self.options.predictors)+'_C'+"-None."+dist
			else: 				      dex_prefix = self.options.prefix+'-'+self.key_nick+'-'+'-'.join(self.options.predictors)+'_C-'+"-".join(self.options.covariates)+'.'+dist
					
			wPreds, wVals = open(dex_prefix+'.dex.preds','w'), open(dex_prefix+'.dex.stats','w') #open(dex_prefix+'.dex.summary','w')
			wPreds, wVals, wSummary = open(dex_prefix+'.dex.preds','w'), open(dex_prefix+'.dex.stats','w'), open(dex_prefix+'.dex.summary','w')
			
			if self.options.permutepvs: wPreds.write(str(pred_format_str+' %15s %10s\n') % ('---','variable','model-desig','v_type','v_size','bw','tval','pval','fdr','permutepv','trials'))
			else:  			    wPreds.write(str(pred_format_str+'\n') % ('---','variables','model-desig','v_type','v_size','bw','tval','pval','fdr'))
			wVals.write(vals_format_str  %  ('---','pred','min-pv','min-label','groups','obs_rates','averages','Fold-Change','Chi-Lo','Chi-Hi'))	


			M = rt.RegModel(self.X,dist,self.options,self.progress,True).run(self.Y,self.feature_names) 
			fdrs = [mlt.fdrcorrection([M.out['params'][i][n][0] for i in range(len(self.feature_names))])[1] for n in range(len(M.out['params'][0]))]
			M_params = [[[fdrs[j][i]]+list(M.out['params'][i][j]) for j in range(len(fdrs))] for i in range(len(self.feature_names))]
			
			if self.options.permutepvs: 	M.run_permutation_tests(self.V)
			if self.options.saveresids:	self.write_resids(M.get_resids(),dex_prefix)

			for j,f in enumerate(self.feature_names):	
				for fdr,pv,tv,bw,n,pb in sorted(M_params[j],key=lambda X:X[-2]):
					if self.options.permutepvs:
						if pb: wPreds.write(str(pred_format_num+' %15.3e %10s\n') % (f,n,self.X.i_type[n],self.X.v_type[n],self.X.group_size[n],bw,tv,pv,fdr,M.permutations[n][j][0],M.permutations[n][j][1]))
						else:  wPreds.write(str(pred_format_num+' %15s %10s\n') % (f,n,self.X.i_type[n],self.X.v_type[n],self.X.group_size[n],bw,tv,pv,fdr,'NA','NA'))
					else: 
						wPreds.write(str(pred_format_num+'\n') % (f,n,self.X.i_type[n],self.X.v_type[n],self.X.group_size[n],bw,tv,pv,fdr))
				for p in [px for px in self.X.parents if px != 'intercept']:


					fdr,pv,tv,bw,n,pb  = sorted([xp for xp in M_params[j] if xp[-2].split('~')[0] == p])[0] 
					if dist[-3::].upper() == 'LOG':         ps = M.X.segregator.summarize(p,[math.log(yc+1.0,2) for yc in self.Y[j]])
					else: 					ps = M.X.segregator.summarize(p,self.Y[j])		
					wVals.write(vals_format_num%(f,p,pv,n,ps.out['names'],ps.out['obs'],ps.out['avg'],ps.fc,ps.chiLo,ps.chiHi))

				for pvf,pvr,tv,pdir,pname,pbool in  M_params[j][1::]:
					if pvf < FDR_CUTOFF:	pv_summary[pname.split('~')[0]][pname.split('~')[-1]].append((pvf,(pdir>0),j))
 
			self.write_model_summary(M.out,pv_summary,dist,dex_prefix) 
			return








	def evaluate_model(self):

                self.progress.start_minor('Running Model Regressions',len(self.D.features),False)
		for dist in self.dists:
			M 	= rt.RegModel(self.X,dist,self.options,self.progress,True).run(self.Y,self.feature_names).aggregate(True) 

			#M_out = rage_regression_outputs.eval_output(self.options).write(M,self.feature_names) #[f.name for f in self.D.features])		
			M_resids, C_resids = M.get_resids() 
			Mc 	= rt.RegModel(self.Xc,dist,self.options).run(self.Y,self.feature_names).aggregate(True)
			#Mp 	= rt.RegModel(self.Xp,dist,self.options).run(self.Y,self.feature_names).aggregate(True) 

			sims = dd(list) 
			self.progress.start_minor('Running Model PCA',len(self.D.features),False)

			dim = rage_DR.DR(self.options,self.progress)#.set_fit_matrix(self.D.matrix('log'))
			pca_init =    dim.set_y_matrix(self.Y,  LOG_TRANSFORM = True,SCALE=True).pca(req='brief') 
			pca_c_resid = dim.set_y_matrix(C_resids,LOG_TRANSFORM = dist[-3::] != 'LOG',SCALE=True).pca(req='brief') 
			pca_resid =   dim.set_y_matrix(M_resids,LOG_TRANSFORM = dist[-3::] != 'LOG',SCALE=True).pca(req='brief')




#			pca_detail = rage_regression_plots.detail_pca(self.D.samples,self.options,3,2)
#			pca_detail.add_pts(pca_init,{'comps': [0,1], 'LEGEND': True,'title':'INIT'}) 
#			pca_detail.add_pts(pca_init,{'comps': [2,3], 'title':'INIT'}) 
#			pca_detail.add_pts(pca_c_resid,{'comps': [0,1], 'title':'COVAR-RESIDS'}) 
#			pca_detail.add_pts(pca_c_resid,{'comps': [2,3], 'title':'COVAR-RESIDS'}) 
#			pca_detail.add_pts(pca_resid,{'comps': [0,1], 'title':'MODEL'}) 
#			pca_detail.add_pts(pca_resid,{'comps': [2,3], 'title':'MODEL'}) 
#			pca_detail.save(dist,self.options.predictors,self.options.covariates)

			for n in range(self.options.simulations): 
				self.progress.start_major('Running Simulation '+str(n+1),False)
				Xs = self.V.select_variables(self.V.variables, permute=self.V.predictors)
				Xs.zp =  self.X.zp 
				Ms = rt.RegModel(Xs,dist,self.options,self.progress).run(self.Y,self.feature_names).aggregate()  #AGGREGATE=True) 
				S_out = rage_regression_outputs.eval_output(self.options).write(Ms,self.feature_names,n+1) #[f.name for f in self.D.features],n+1)
				sims['pv'].append(Ms.pv_cnt)
				sims['rs'].append(Ms.rs_cnt) 
				sims['v_exp'].append(np.mean(Ms.out['v_exp'])) 	
			




			self.progress.start_minor('Plotting Results  ',100,False)
			mplot = rage_regression_plots.model_plot(self.D.samples,self.X,self.options,3,2,{'p_key': M.pv_key,'r_key': M.rs_key}) 

			#mplot.add_model_table(M, total_var = [M.total_var,Mc.total_var,np.mean(sims['var']) ]).update() #[pca_init['total_var'], pca_c_resid['total_var'], pca_resid['total_var'], sims['var']]).update()
			#mplot.add_model_table(M, total_var = [100,100,100]).update() #,Mc.total_var,np.mean(sims['var']) ]).update() #[pca_init['total_var'], pca_c_resid['total_var'], pca_resid['total_var'], sims['var']]).update()
			mplot.add_model_table(M, total_var =  [np.mean(Mc.out['v_exp']), np.mean(M.out['v_exp']), np.mean(sims['v_exp'])]).update() 
			#mplot.add_model_table(M, skrees = skrees).update() #[pca_init['total_var'], pca_c_resid['total_var'], pca_resid['total_var'], sims['var']]).update()
			#mplot.add_model_table(M,self.X, skrees = [pca_init['total_var'], pca_c_resid['total_var'], pca_resid['total_var'], sims['var']]).update()
			#mplot.add_predictor_table(M,Mp,Mc,self.X,self.options,{'sim_pvs': sims['pv']}).update() 
			mplot.add_predictor_table(M,self.X,self.options,{'sim_pvs': sims['pv']}).update() 

			mplot.add_rs_bars(M.rs_cnt,self.options,sims['rs']).update({'title': '$'+"\ ".join(self.V.predictors)+'$ '+'$\  R^2\ Values$'})
			mplot.add_pv_bars(M.pv_cnt,self.options,sims['pv']).update({'title': '$'+"\ ".join(self.V.predictors)+'$ '+'$\  P\ \ Values$'})
			mplot.add_pca_pts(pca_init,{'colspan':2}).update({'title': 'PCA Initial Values','yadd': 2,'colspan':2})
			mplot.add_pca_pts(pca_c_resid,{'colspan':2}).update({'title': 'PCA Covariate Residuals','yadd': 2,'colspan':2})
			mplot.add_pca_pts(pca_resid,{'colspan':2}).update({'title': 'PCA Model Residuals','yadd': 2,'colspan':2})
			mplot.save(dist,self.options.predictors,self.options.covariates)
		
			rage_regression_outputs.reg_simulate(self.options).write(M.pv_cnt, sims['pv'], self.options.simulations, self.V.predictors)

			
			self.progress.end() 



























































































	def eval_covariates(self):

                self.progress.start_minor('Testing Base Predictor Model',len(self.D.features))
		sims = dd(list) 			

		M 	= rt.RegModel(self.options,self.X,True).run(self.Y).aggregate(True) 

		
		Mc 	= rt.RegModel(self.options,self.Xc,True).run(self.Y).aggregate(True)
		Mp 	= rt.RegModel(self.options,self.Xp,True).run(self.Y).aggregate(True) 




#		rage_outputs.regression_result(self.options).write(self.D,M,Mp,suffix='dexcovar') 
		pv_discs,minPv = [0.05,0.025], min(Mp.pv_mins)
		if minPv < 0.01:
			if np.percentile(Mp.pv_mins,0.5) > 0.001: 	pv_discs = [0.05,0.01]
			elif np.percentile(Mp.pv_mins,0.5) > 0.0001: 	pv_discs = [0.05,0.001]
			else: 					  	pv_discs = [0.05,0.0001] 

		sims,covar_data = [rt.RegModel(self.options).score(self.Y,self.V.select_variables(self.V.predictors, permute = self.V.predictors)) for s in range(self.options.simulations)],{}
		for i in self.X.covariate_idx:

                	self.progress.start_minor('Testing Covariate: '+ self.X.names[i],len(self.D.features))
			cLen,cType,cName = len(self.D.samples) , 'continuous',self.X.names[i]
			cParent = self.X.parent[self.X.names[i]]
			cSeg,cLens = self.D.samples.segregate(cParent) 

			if self.D.samples.attribute_class[cParent] == 'binary': cLen,cType = cLens[self.X.names[i]],'categorical'
			Xc = self.V.select_variables([self.X.names[i]])
			Ci = rt.RegModel(self.options).score(self.Y,Xc) #.summarize()

			Xi= self.V.select_variables(self.V.predictors+[self.X.names[i]])
			Cm = rt.RegModel(self.options).score(self.Y,Xi) #.summarize()
			mDisc = Mp.pv_discovery(Cm.pv_mins,pv_discs)	

			cSims = [rt.RegModel(self.options).score(self.Y,self.V.select_variables(self.V.predictors+[self.X.names[i]], permute = self.V.predictors)) for s in range(self.options.simulations)]
			sDisc = {pv_discs[j]: round(np.mean(sD[pv_discs[j]]),3) for sD in [mS.pv_discovery(cS.pv_mins,pv_discs) for mS,cS in zip(sims,cSims)] for j in range(len(pv_discs))}
			covar_data[cName] = [cType,cLen,Cm.vif[cName],round(Ci.stats['rs'].mean,3),round(Cm.stats['rs'].mean,3)],mDisc,sDisc

		self.progress.end() 
                self.progress.start_minor('Plotting Results  ',100,False)

		eplot = rage_regression_plots.eplot(self.options,pv_discs)
		eplot.add_base_model(Mp,self.Xp.names,sims)

		eplot.add_covariate_data(covar_data)
		eplot.save_efig(self.options.model,self.V.predictors,self.V.covariates)














	def write_resids(self,residuals,dex_prefix): 
		M_resids, C_resids = residuals 
		wC, wM = open(dex_prefix+'.cResids.cnts','w'), open(dex_prefix+'.mResids.cnts','w')
		wC.write('%-30s %s\n' % ('---',"\t".join([s.name for s in self.D.samples])))
		wM.write('%-30s %s\n' % ('---',"\t".join([s.name for s in self.D.samples])))
		for j,f in enumerate(self.feature_names):
			resid_str = "\t".join([str(round(x,3)) for x in C_resids[j]])
 			wC.write('%-30s %50s\n' % (f,resid_str))
			m_str = "\t".join([str(round(x,3)) for x in M_resids[j]])
 			wM.write('%-30s %50s\n' % (f,resid_str))
		wC.close(); wM.close() 



	def write_model_summary(self,model_output,model_summary,dist,dex_prefix):
			
			
		wSummary = open(dex_prefix+'.dex.summary','w')
		FL = range(len(self.feature_names)) 
		if len(self.options.covariates) > 0:  wSummary.write('#MODEL-SUMMARY  DIST,PREDICTORS,COVARIATES        %10s %20s %20s\n' % (dist,','.join(self.options.predictors),",".join(self.options.covariates))) 
		else: 				      wSummary.write('#MODEL-SUMMARY  DIST,PREDICTORS,COVARIATES        %10s %20s %20s\n' % (dist,','.join(self.options.predictors),'None'))
		wSummary.write('#MODEL-RESULT   Rsq,Rsa,Var-Exp,BIC,ZeroInflation') 
		for k_key,k_name in zip(['rsq','rsa','v_exp','bic','zero_infl'],['R-Squared','R-squared-adj','Variance-Explained','BIC','Zero-Inflation']):
			try: 	 wSummary.write(' %10.5f ' %  (np.mean([model_output[k_key][j] for j in FL])))
			except:  wSummary.write(' %10s ' %  ('NA'))
		wSummary.write('\n') 
		for K in model_summary.keys():
			k_dex,k_set = dd(lambda: dd(list)), [] 
			for k in model_summary[K].keys():
				k_sort =  sorted(model_summary[K][k]) 
				k_set.extend([kx[-1] for kx in k_sort])
				k_dex[k]['UP'],k_dex[k]['DOWN']   = [(kd[0],kd[-1]) for kd in k_sort if kd[1]], [(kd[0],kd[-1]) for kd in k_sort if not kd[1]]
			wSummary.write('%-20s %20s %10d \n' % ('#VARIABLE:'+K,'DEX-GENES',len(list(set(k_set)))))	
			for k in model_summary[K].keys(): 
				for st in ['UP','DOWN']: 
					kBrief = k_dex[k][st][0:10]
					uStr = ",".join([self.feature_names[ki[-1]].split(';')[-1] for ki in kBrief])
					try: 	 wSummary.write('%-20s %20s %10d %120s  %12.2e...%5.2e\n' % ('#VARIABLE:'+K+'='+k,st,len(k_dex[k][st]),uStr,kBrief[0][0],kBrief[-1][0]))
					except:  wSummary.write('%-20s %20s %10d %120s %12s\n' % ('#VARIABLE:'+K+'='+k,st,len(k_dex[k][st]),'None','NA...NA'))				
		wSummary.write('\n') 

			












































































 















	def eval_predictors(self):



		dim = rage_DR.DR(self.options,self.progress)#.set_fit_matrix(self.D.matrix('log'))
		predictor_plot = rage_regression_plots.predictor_plot(self.options,len(self.rage.args.predictors))

#		reg_out = rage_outputs.regression_output(self.options,M_full) 
#		pred_out = rage_outputs.predictor_output(self.options)

		for p in self.rage.args.predictors: 

			self.D = copy.deepcopy(self.rage.data) 
			self.D.rage.progress.reset()	

			#self.D.filter_samples_by_attributes([p],[]).normalize() 			
			#self.V = self.D.set_sample_variables([p]) 

			self.D.filter_samples_by_attributes([p],self.options.covariates).normalize() 			
			self.V = self.D.set_sample_variables([p],self.options.covariates) 

			self.Y = [[s.cnts[f.idx] for s in self.D.samples] for f in self.D.features]
			self.X = self.V.select_variables(self.V.variables)

			self.options.color = [p] 
			self.D.samples.create_plot_labels(self.options) 

                	self.progress.start_minor('Running Predictor Regression: '+p,len(self.D.features),False)
			self.progress.mark() 
			M = rt.RegModel(self.options,self.X,True).run(self.Y).aggregate(True) 
			pca_init    = dim.pca(self.D.matrix(),req='brief') #['total_var'] 
			pca_resid =	dim.pca(np.matrix(M.out['resids']).getT(),req='brief') #@['total_var']


			sims = dd(list) 
			SUMMARIZE=True
			
			preds = [sp for sp in M.pv_dict.keys() if sp != 'intercept'] 
			gt_key = dd(lambda: dd(int)) 
			best_key = dd(lambda: dd(float)) 
			for n in range(self.options.simulations): 
                		self.progress.start_minor('Running Simulation '+str(n+1),False)
				Ms = rt.RegModel(self.options,self.V.select_variables(self.V.variables, permute = [p])).run(self.Y).aggregate(True) 
				sims['var'].append(dim.pca(np.matrix(Ms.out['resids']).getT(),req='brief')['total_var'])
				sims['pv'].append(Ms.pv_cnt)
				sims['rs'].append(Ms.rs_cnt) 

				if SUMMARIZE: 
					for i,(f,Yi) in enumerate(zip(self.D.features,self.Y)): 
						rsq,rsa = Ms.out['rsq'][i], Ms.out['rsa'][i] 
						if rsq > M.out['rsq'][i]: gt_key[i]['rsq'] += 1
						if rsa > M.out['rsa'][i]: gt_key[i]['rsa'] += 1
						if rsa > best_key[i]['rsa']: best_key[i]['rsa'] = rsa
						if rsq > best_key[i]['rsq']: best_key[i]['rsq'] = rsq
						for sp in M.pv_dict.keys(): 
							spV = Ms.pv_dict[sp][i] 
							if spV < M.pv_dict[sp][i]: gt_key[i][sp] += 1
							if spV < best_key[i][sp]: best_key[i][sp] = spV
						

			pred_out = rage_outputs.predictor_output(self.options,p,M,self.D.features,self.Y)
			pred_out.add_sim_keys(gt_key,best_key,self.options.simulations) 

			predictor_plot.add_predictor_row(p,self.D.samples,pca_init,pca_resid,M,sims) 	
			self.progress.end() 
		predictor_plot.save(self.options.prefix+"-predictorplot-"+"_".join(self.rage.args.predictors)+'-cov-'+'_'.join(self.rage.args.covariates))
		sys.exit() 





















	def evaluate_dex2(self,permutations=100): 

                self.progress.start_minor('Running Dex Tests',len(self.D.features),False)

#		if self.options.dist == 'OLS':  Y = [[math.log(s.cnts[f.idx]+1.0,2) for s in self.D.samples] for f in self.D.features]
#		else:				Y = [[s.cnts[f.idx] for s in self.D.samples] for f in self.D.features]
		
		if self.options.dist == 'OLS': 
			self.D = self.D.normalize() 
		Y = [[s.cnts[f.idx] for s in self.D.samples] for f in self.D.features]
		M_full	= rt.RegModel(self.options,self.X,True).run(Y).aggregate(True) 
		M_raw 	= rt.RegModel(self.options,self.Xp,True).run(Y).aggregate(True) 
		self.segregation_key = {p: self.D.samples.segregate(p) for p in self.V.predictors} 
#		output = rage_outputs.dex_result(self.options,",".join([s.name for s in self.D.samples]))  

		reg_out = rage_outputs.regression_output(self.options,M_full) 
		
		M_full.run_permutation_tests(Y,self.V) 

		for i,(f,Yi) in enumerate(zip(self.D.features,Y)): 



			yLen,yObs,yRate,y_new = len(Yi),len([y for y in Yi if y>0]), len([y for y in Yi if y>0])/float(len(Yi)), M_full.out['covariate_resids'][i] 
			M_ALL, M_PRED = sorted(M_full.out['params']), [x for x in sorted(M_full.out['params']) if x[-1]]
			reg_out.write_feature_summary(f,Yi,i) 
#			continue 
			for k,(pv,bw,n,PRED) in enumerate(sorted(M_full.out['params'][i])):
				if not PRED: continue 

				if len(n.split('~')) > 1:
					nI = self.segregation_key[n.split('~')[0]][0][n]
					yN,rN = [Yi[ni] for ni in nI], [y_new[ni] for ni in nI]
					nLen, nObs,yXN,rXN = len(yN), float(len([yn for yn in yN if yn > 0])) / len(yN) , np.mean(yN), np.mean(rN)  
					#reg_out.add_binary_predictor(n,M_raw.pv_dict[n][i],pv,bw,M_full.permutations[n][i],yN,rN)
					#reg_out.add_binary_predictor(n,M_raw.pv_dict[n][i],pv,bw,M_full.permutations[n][i],yN,rN)
					reg_out.add_binary_predictor(n,M_raw.pv_dict[n][i],pv,bw,M_full.permutations[n][i],yN,rN)
				
				else:
					reg_out.add_cont_predictor(n,M_raw.pv_dict[n][i],pv,bw,M_full.permutations[n][i])
										

























	def evaluate_model2(self):

                self.progress.start_minor('Running Model Regressions',len(self.D.features),False)
		for dist in self.options.dist:


			print 'yo' 
			M 	= rt.RegModel(self.X,dist,self.options,self.progress,True).run(self.Y,self.feature_names).aggregate(True) 

			print 'yo' 

			M_resids, C_resids = M.get_resids(COVAR=True) 


			print 'yo' 
			M_out = rage_regression_outputs.eval_output(self.options).write(M,self.feature_names) #[f.name for f in self.D.features])		
			Mc 	= rt.RegModel(self.Xc,dist,self.options).run(self.Y,self.feature_names).aggregate(True)
			#Mp 	= rt.RegModel(self.Xp,dist,self.options).run(self.Y,self.feature_names).aggregate(True) 

			print 'yo' 
			sims = dd(list) 
			self.progress.start_minor('Running Model PCA',len(self.D.features),False)

			dim = rage_DR.DR(self.options,self.progress)#.set_fit_matrix(self.D.matrix('log'))
			pca_init =    dim.set_y_matrix(self.Y,  LOG_TRANSFORM = dist[-3::] != 'LOG').pca(req='brief') 
#			pca_c_resid = dim.set_y_matrix(C_resids,LOG_TRANSFORM = dist[-3::] != 'LOG').pca(req='brief') 
			pca_resid =   dim.set_y_matrix(M_resids,LOG_TRANSFORM = dist[-3::] != 'LOG').pca(req='brief')



