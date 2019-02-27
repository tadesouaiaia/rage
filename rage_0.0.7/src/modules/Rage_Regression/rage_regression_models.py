# 2019.02.06 16:19:20 EST
#Embedded file name: /home/tade/RAGE_PACKAGE/code/rage_0.0.7/rage_0.0.7/src/modules/Rage_Regression/rage_regmodels.py



from collections import Counter as cc
from collections import defaultdict as dd
from scipy.misc import factorial 
from scipy.stats import nbinom
from scipy.stats import poisson
from sklearn.preprocessing import MinMaxScaler
from statsmodels.base.model import GenericLikelihoodModel
from statsmodels.stats import power as smp
from statsmodels.stats.multitest import fdrcorrection as fdr
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif
import math
import numpy as np
import os
import random
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.discrete.count_model as scm 
import statsmodels.discrete.discrete_model as sdm
import statsmodels.genmod.families as sfams
import statsmodels.genmod.families.links as slinks
import statsmodels.stats.multitest as mlt
import sys
import warnings
warnings.filterwarnings("ignore")



















def scale_vals(vals, f1 = 0, f2 = 1):
    scaler = MinMaxScaler(feature_range=(f1, f2))
    return scaler.fit_transform(np.array(vals, dtype=float).reshape(-1, 1)).reshape(1, -1)[0]


def vif_test(X):
    vd, vd_out = dd(list), {}
    for i, n in enumerate(X.names):
        try:
            vd_out[n] = round(vif(X.array, i), 4)
        except:
            vd_out[n] = 0.0

    return vd_out



class RegModel():
    def __init__(self, X, dist, options, progress = False, FULL = True):
        self.X, self.dist, self.options, self.progress, self.vif = X,dist,options,progress,vif_test(X) 
        self.alp1, self.alp2 = (0.05, 0.001)
        self.out = dd(list)
        self.rs_key = [0.01,0.03,0.05,0.1,0.25]
        self.pv_key = [0.05,0.01,0.001,0.0001,1e-05,1e-07]
        self.permutations = dd(list)
        self.resid = []
        self.regress = RegTest(X, self.dist, alphas=[self.alp1, self.alp2])
	self.residuals = RegResids(X, self.dist) 





    





    def run(self, Y, Ynames = None):
        t_var, d_var, self.y_names = 0, 0, Ynames

        if self.dist[-3::].upper() == 'LOG':	self.Y = [ [ math.log(yc + 1.0, 2) for yc in y ] for y in Y ]
        else:	self.Y = Y


        if self.progress:	self.progress.start_minor('Running ' + self.dist + ' Regression', len(self.Y))
        for yi, y in enumerate(self.Y):
            if self.progress:	self.progress.mark()
            t = self.regress.test(y)
            self.out['params'].append(t.output)
            if len(self.X.predictor_idx) > 0:	self.out['pvs'].append([ p[0] for p in sorted(t.output) if p[-1] ][0])
            for r, k in zip([t.valid,t.history,t.zero_infl,t.v_explained,t.rsq,t.rsa,t.bic,t.pwr[self.alp1],t.pwr[self.alp2]], ['VALID','history','zero_infl','v_exp','rsq','rsa','bic','pwr1','pwr2']): self.out[k].append(r)

            self.out['run'].append(t.res)

        fdrs = [ mlt.fdrcorrection([ self.out['params'][i][n][0] for i in range(len(self.Y)) ])[1] for n in range(len(self.X.names)) ]
        self.results = [ [ [fdrs[j][i]] + list(self.out['params'][i][j]) for j in range(len(fdrs)) ] for i in range(len(self.Y)) ]

        self.param_key = ['fdr','pval','tval','bw','name','pred_bool']
        return self



    def aggregate(self, PARAMS = False):
        if len(self.X.predictor_idx) > 0: pv_list = [ sorted([ p[0] for p in P if p[-1] ])[0] for P in self.out['params'] ]
        else:				  pv_list = []
        self.pv_cnt = [ len([ p for p in pv_list if p < self.pv_key[j] ]) for j in range(len(self.pv_key)) ]
        self.rs_cnt = [ len([ p for p in self.out['rsq'] if p > self.rs_key[j] ]) for j in range(len(self.rs_key)) ]
        if PARAMS:	self.pv_dict = {self.X.names[i]:[ p[i][0] for p in self.out['params'] ] for i in range(len(self.X.names))}
        return self



    def run_permutation_tests(self, V):
        maxT = self.options.maxtrials
        if maxT < 100: maxT = 100 
        if maxT > 1000: 
            minT = 100
            midT = 200
        elif maxT > 499:
            minT = 50
            midT = 250
        else:
            minT = 10
            midT = 50
        minCut, midCut = 0.05 * minT, 0.01 * (midT + minT)


        for i, y in enumerate(self.Y):
            p_cands, self.p_key = [ [a[3], a[0]] for a in [ r for r in self.out['params'][i] if r[-1] ] ], {}
            for c, pv in p_cands:
                if pv > 0.1:
                    self.permutations[c].append((round(pv, 2), 1))
                else:
                    self.p_key[c] = [pv, 0, 0]

            if len(self.p_key) > 0:
                self.regress.permute(y, [ V.select_variables(V.variables, permute=V.predictors) for j in range(minT) ], self.p_key)
                for c, (pv, L, G) in self.p_key.items():
                    if L > minCut:
                        self.permutations[c].append((self.p_key.pop(c)[1] / float(L + G), L + G))

                if len(self.p_key) > 0:
                    self.regress.permute(y, [ V.select_variables(V.variables, permute=V.predictors) for j in range(midT) ], self.p_key)
                    for c, (pv, L, G) in self.p_key.items():
                        if L > midCut:
                            self.permutations[c].append((self.p_key.pop(c)[1] / float(L + G), L + G))

                    if len(self.p_key) > 0:
                        self.regress.permute(y, [ V.select_variables(V.variables, permute=V.predictors) for j in range(maxT) ], self.p_key)
                        for c, (pv, L, G) in self.p_key.items():
                            self.permutations[c].append((self.p_key.pop(c)[1] / float(L + G), L + G))

















    def get_resids(self, COVAR = True, SCALE = True):
        #self.resids, self.p_resids, self.c_resids = [], [], []
        #rIdx, cIdx, pIdx = range(len(self.X.names)), [ i for i in range(len(self.X.names)) if self.X.COVARIATE[self.X.names[i]] ], [ i for i in range(len(self.X.names)) if self.X.PREDICTOR[self.X.names[i]] ]


	self.residuals.retrieve(self.Y) 
	sys.exit() 

        if self.dist[0:3].upper() == 'OLS':
            for i, y in enumerate(self.Y):

                print y
                sys.exit()
                mR = [ y[s] - sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in rIdx ]) for s in range(len(y)) ]
                cR = [ y[s] - sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in cIdx ]) for s in range(len(y)) ]
                if SCALE:
                    sc = MinMaxScaler(feature_range=(0, max(y)))
                    mR = [ x[0] for x in sc.fit_transform(np.array(mR).reshape(-1, 1)) ]
                    cR = [ x[0] for x in sc.fit_transform(np.array(cR).reshape(-1, 1)) ]
                self.resids.append(mR)
                self.c_resids.append(cR)

            return (self.resids, self.c_resids)
        SCALE = True
        for i, y in enumerate(self.Y):
            if not self.out['VALID'][i]:
                self.resids.append(y)
                self.p_resids.append(y)
                self.c_resids.append(y)
                continue
            r_pred = [ np.exp(sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in rIdx ])) for s in range(len(y)) ]
            p_pred = [ np.exp(sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in pIdx ])) for s in range(len(y)) ]
            c_pred = [ np.exp(sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in cIdx ])) for s in range(len(y)) ]
            if SCALE:
                try:
                    sc = MinMaxScaler(feature_range=(0, max(y)))
                    r_pred = [ x[0] for x in sc.fit_transform(np.array(r_pred).reshape(-1, 1)) ]
                    p_pred = [ x[0] for x in sc.fit_transform(np.array(p_pred).reshape(-1, 1)) ]
                    c_pred = [ x[0] for x in sc.fit_transform(np.array(c_pred).reshape(-1, 1)) ]
                except ValueError:
                    print 'FAIL'
                    print self.y_names[i]
                    print max(y), min(y)
                    sys.exit()

            r_resid = [ ys - r_pred[s] for s, ys in enumerate(y) ]
            p_resid = [ ys - r_pred[s] for s, ys in enumerate(y) ]
            c_resid = [ ys - r_pred[s] for s, ys in enumerate(y) ]
            if min(r_resid) < 0:
                r_resid = [ rr + -1 * min(r_resid) for rr in r_resid ]
            if min(p_resid) < 0:
                p_resid = [ rr + -1 * min(p_resid) for rr in p_resid ]
            if min(c_resid) < 0:
                c_resid = [ rr + -1 * min(p_resid) for rr in c_resid ]
            self.resids.append(r_resid)
            self.p_resids.append(p_resid)
            self.c_resids.append(c_resid)

        return (self.resids, self.c_resids)







class RegResids:
	def __init__(self,X,dist='OLS'): 
        	

		#self.resids, self.p_resids, self.c_resids = [], [], []
        	#rIdx, cIdx, pIdx = range(len(self.X.names)), [ i for i in range(len(self.X.names)) if self.X.COVARIATE[self.X.names[i]] ], [ i for i in range(len(self.X.names)) if self.X.PREDICTOR[self.X.names[i]] ]

		self.cIdx =  [i for i,n in enumerate(X.names) if n != 'intercept' and X.COVARIATE[n]]
		self.mIdx =  [i for i,n in enumerate(X.names) if n != 'intercept'] 

		if 'OLS' in dist: 
			self.retrieve = self.get_OLS 




	def get_OLS(self, Y, COVAR = True, SCALE = True):


		print 'hi' 
		sys.exit() 

		if self.dist[0:3].upper() == 'OLS':
		    for i, y in enumerate(self.Y):

			print y
			sys.exit()
			mR = [ y[s] - sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in rIdx ]) for s in range(len(y)) ]
			cR = [ y[s] - sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in cIdx ]) for s in range(len(y)) ]
			if SCALE:
			    sc = MinMaxScaler(feature_range=(0, max(y)))
			    mR = [ x[0] for x in sc.fit_transform(np.array(mR).reshape(-1, 1)) ]
			    cR = [ x[0] for x in sc.fit_transform(np.array(cR).reshape(-1, 1)) ]
			self.resids.append(mR)
			self.c_resids.append(cR)

		    return (self.resids, self.c_resids)


















class RegTest:
	def __init__(self,X,dist='OLS',alphas=[0.05,0.01],log=True):

		self.X,self.xLen,self.dist,self.permute,self.zero_infl,self.alphas,self.dfd,self.dfn, self.LOG =  X, len(X.names), dist.upper(), self.permute_REG, 0.0, alphas, len(X.names) -1,  len(X.array) - len(X.names),False
		self.test_history = []  



		if self.dist[0] == 'Z': 
			Z_KEY = {'ZGP':(scm.ZeroInflatedGeneralizedPoisson, scm.GeneralizedPoisson),'ZNB':(scm.ZeroInflatedNegativeBinomialP,scm.NegativeBinomialP),'ZPO':(scm.ZeroInflatedPoisson,scm.Poisson)}
			self.execute, self.permute = self.execute_zif, self.permute_ZIF
			self.zi_model, self.nz_model = Z_KEY[self.dist]	
		elif self.dist in ['NB','GP','PO']:
			F_KEY = {'GP': scm.GeneralizedPoisson,'NB':scm.NegativeBinomialP,'PO':scm.Poisson}
			self.execute, self.nz_model, self.permute = self.execute_cm , F_KEY[self.dist], self.permute_CM

		elif self.dist[0:3] == 'OLS': 
			self.execute, self.LOG, self.permute = self.execute_ols, (dist[-3::].upper() == 'LOG'), self.permute_OLS 
			
		else:
			print "UNSUPPORTED"
			sys.exit() 


	def test(self,y): 

		#print self.LOG,'huh' 

		#if self.LOG: y = [math.log(yi+1.0,2) for yi in y]



		self.output,self.zero_infl, self.rsq, self.rsa, self.bic, self.aic  = [],0.0,'NA','NA','NA','NA'
		self.valid, self.y, self.yA, self.yLen, self.history  = True, y, np.array(y), len(y), '' 
		self.execute() 

		if self.valid: 
			self.v_explained = 1-(np.var(self.res.resid) / np.var(self.yA)) 
			try : self.pwr = {a: smp.FTestPower().solve_power(effect_size=np.sqrt(self.v_explained/(1-self.v_explained)),df_num=self.dfn,df_denom=self.dfd,alpha=a) for a in self.alphas}
			except: self.pwr = {a: 0.5 for a in self.alphas}
			if any([np.isnan(pw) for pw in self.pwr.values()]): self.pwr = {a: 0.5 for a in self.alphas}
		else:
			self.v_explained = 0 
			self.pwr = {a: 0.0 for a in self.alphas} 
			self.output = [(0.5,b,t,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues,self.res.tvalues, self.res.params, self.X.names))]
			self.bic, self.aic, self.rsq, self.rsa = 0, 0, 0, 0 
			 
			self.tvalues = self.res.tvalues

			#self.bic, self.aic, self.rsq, self.rsa = self.res.bic, self.res.aic, self.res.prsquared, 1- (((1-self.res.prsquared)*(self.yLen-1)) / self.dfn) #(self.yLen-self.X.len-1))
			#self.output = [(p,b,x,i in self.X.predictor_idx) for i,(p,b,x) in enumerate(zip(self.res.pvalues, self.res.params, self.X.names))]

		return self


	def execute_ols(self): 

		self.history += 'ols-fit'
		self.res = sm.OLS(self.yA, self.X.array).fit(disp=0) 

		self.fval,self.fpv, self.tvalues = self.res.fvalue, self.res.f_pvalue, self.res.tvalues
		self.rsq, self.rsa, self.bic, self.aic  = round(self.res.rsquared,5), round(self.res.rsquared_adj,3), round(self.res.bic,3), round(self.res.aic)
		 
		self.output = [(p,t,b,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues, self.res.tvalues, self.res.params, self.X.names))]
		self.history += '-sucess'

	def execute_zif(self): 
		self.reg = self.zi_model(self.yA, self.X.array,exog_infl=self.X.zp) 
		self.history += 'zif-fit'

		try: 
			self.res = self.reg.fit(disp=0) 	
			if any([np.isnan(pv) for pv in self.res.pvalues]): 
				self.history += '-fail,zif-regularized'
				self.res = self.reg.fit_regularized(disp=0) 
				if any([np.isnan(pv) for pv in self.res.pvalues]): raise np.linalg.linalg.LinAlgError('NAN Vals') 
		except np.linalg.linalg.LinAlgError:	self.history += '-fail'
		except AssertionError:			self.history += '-fail'
		if self.history[-4::] == 'fail': 
			self.execute_cm() 
			return 
	
		
		self.history += '-sucess'
		self.v_explained = 1-(np.var(self.res.resid) / np.var(self.yA)) 
		self.zero_infl = 1.0 - (1.0 / (1 + np.exp(self.res.params[0])))	
		self.tvalues = self.res.tvalues
		self.bic, self.aic, self.rsq, self.rsa = self.res.bic, self.res.aic, self.res.prsquared, 1- (((1-self.res.prsquared)*(self.yLen-1)) / self.dfn) #(self.yLen-self.X.len-1))
		self.output = [(p,t,b,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues[1::],self.res.tvalues,self.res.params[1::], self.X.names))]

			



	def execute_cm(self):

		self.history += ',cm-fit'
		FAIL = False
		try: 
			self.res = self.nz_model(self.yA, self.X.array).fit(disp=0)
			if any([np.isnan(pv) for pv in self.res.pvalues]): raise np.linalg.linalg.LinAlgError('NAN Vals') 
		except np.linalg.linalg.LinAlgError:	self.history += '-fail'
		except AssertionError:			self.history += '-fail'
			

		if self.history[-4::] == 'fail': 
			self.history += ',cm-regfit'
			try: 
				self.res = self.nz_model(self.yA, self.X.array).fit_regularized(disp=0)
				if any([np.isnan(pv) for pv in self.res.pvalues]): raise np.linalg.linalg.LinAlgError('NAN Vals') 
			except np.linalg.linalg.LinAlgError:	self.history += '-fail'
			except AssertionError:			self.history += '-fail'

		
		if self.history[-4::] == 'fail':
			self.valid = False 
			return 
		
		if self.history[-4::] != 'fail': 
 
			self.history += '-sucess'
			self.tvalues = self.res.tvalues

			self.bic, self.aic, self.rsq, self.rsa = self.res.bic, self.res.aic, self.res.prsquared, 1- (((1-self.res.prsquared)*(self.yLen-1)) / self.dfn) #(self.yLen-self.X.len-1))
			self.output = [(p,t,b,x,i in self.X.predictor_idx) for i,(p,t,b,x) in enumerate(zip(self.res.pvalues, self.res.tvalues,self.res.params, self.X.names))]
 



	def permute_OLS(self,y,Xp,key):



 
		for P in Xp:
			
			#print P.array
 
			#res = sm.OLS(np.array(y),Xp[0].array).fit(disp=0) 
			res = sm.OLS(np.array(y),P.array).fit(disp=0) 
			for i,(p,n) in enumerate(zip(res.pvalues,  self.X.names)):
				if n in key:

 
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1

				
	def permute_CM(self,y,Xp,key): 
		for P in Xp:
			try:  
				res = self.nz_model(np.array(y), P.array).fit(disp=0)
			except:	
				try: res = self.nz_model(np.array(y), P.array).fit_regularized(disp=0)

				except: return 
					

			
			for i,(p,n) in enumerate(zip(res.pvalues, self.X.names)):
				if n in key: 
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1


	def permute_ZIF(self,y,Xp,key): 
		for P in Xp: 
			try: res = self.zi_model(self.yA, P.array,exog_infl=self.X.zp).fit(disp=0) 
			except: 
				try: res = self.nz_model(np.array(y), P.array).fit(disp=0)
				except: 
					try:    res = self.nz_model(np.array(y), P.array).fit_regularized(disp=0)
					except: res.pvalues = [0.5 for x in self.X.names]



			for i,(p,n) in enumerate(zip(res.pvalues, self.X.names)):
				if n in key: 
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1





	def permute_REG(self,y,Xp,key): 		
		for P in Xp: 

			


			model = self.reg(y, Xp[0].array).fit(disp=0)
			for i,(p,b,n) in enumerate(zip(model.pvalues,model.params,self.X.names)):
				if n in key:
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1

	def permute_GLM(self,y,Xp,key): 
		for P in Xp: 
			model = sm.GLM(y, P.array, family= self.family).fit(disp=0)
			for i,(p,b,n) in enumerate(zip(model.pvalues,model.params,self.X.names)):
				if n in key:
					if p < key[n][0]: key[n][1]+=1
					else:		  key[n][2]+=1





