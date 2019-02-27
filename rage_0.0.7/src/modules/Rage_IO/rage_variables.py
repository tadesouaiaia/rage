import sys
from collections import defaultdict as dd
from collections import Counter as cc
import numpy as np
from math import log
import mmap 
from sklearn.preprocessing import MinMaxScaler
from random import shuffle
import copy 
from scipy.stats import chisquare






def command_line_error(msg):
	sys.stderr.write('\n')
        sys.stderr.write('RageCommandLineError: '+msg+'\n')
        sys.exit()







			
		


class RegVariables:
	def __init__(self,sample_attributes,variable_options,variable_key,predictors,covariates):



		sample_num =  len(sample_attributes.values()[0]) 
		self.group_sizes = {'intercept': sample_num} 
		self.BIN = dd(bool) 



 
		for v in variable_options: 
			if len(variable_options[v][1]) == 0 or type(variable_options[v][1][0]) == str: 
				self.BIN[v] = True 
				for n,x in cc(sample_attributes[v]).items(): 
					self.group_sizes[n] = x 
					
			else:
				self.group_sizes[v] = sample_num 
	



		self.predictors, self.covariates, self.variables = predictors,covariates,predictors+covariates 
		self.PREDICTOR, self.COVARIATE = dd(bool), dd(bool) 
		self.PREDICTOR['intercept'],self.COVARIATE['intercept'] = True, True
		self.names = ['intercept'] 
		self.options, self.types, self.inferred = {'intercept': ['intercept']}, {},  dd(list) 
		
		self.vals = {'intercept' :  [1.0 for s in sample_attributes.values()[0]] }
		self.key = {'intercept': {1.0: [1.0]}}
		self.sample_vals = [[1.0] for s in sample_attributes.values()[0]] 

		for opt in variable_options:


			self.names.append(opt) 
			self.options[opt], self.inferred[opt]  = variable_options[opt][0] , variable_options[opt][1]  
			self.vals[opt] = sample_attributes[opt]
			for i,v in enumerate(self.vals[opt]):  self.sample_vals[i].append(v) 
			self.key[opt] = variable_key[opt] 
			if opt in predictors: self.PREDICTOR[opt] = True
			else:		      self.COVARIATE[opt] = True 
	

	




	def select_variables(self,members,permute=[]): 



		variables = list(set(members+permute))
		shuffle_variables, assigned_variables = [], dd(list) 
		s_parents, s_opts, s_lists = [],[],[[] for s in range(len(self.sample_vals))]	




		for j,v in enumerate(self.names):

					
	
			if v == 'intercept' or v in variables: 
				if v in permute: 	shuffle_variables.append((j,v)) 
				else:


					s_parents.append(v) 
					s_opts.append(self.options[v]) 
					for i in range(len(self.sample_vals)): 
						s_lists[i].append(self.key[v][self.sample_vals[i][j]])
						assigned_variables[v].append(self.sample_vals[i][j]) 



		for v in [v for v in variables if v not in self.names]:
			if v.split('~')[0] in self.names:
				j,V =self.names.index(v.split('~')[0]), v.split('~')[0]
				k=self.options[V].index(v)
				s_parents.append(v)
				s_opts.append([v]) 
				for i in range(len(self.sample_vals)): s_lists[i].append([self.key[V][self.sample_vals[i][j]][k]])


		for (j,v) in shuffle_variables:
			sv =  [self.sample_vals[i][j] for i in range(len(self.sample_vals))] 



			shuffle(sv) 
			s_parents.append(v) 
			s_opts.append(self.options[v]) 
			for i in range(len(self.sample_vals)): 
				s_lists[i].append(self.key[v][sv[i]])
				assigned_variables[v].append(sv[i]) 



		return RegVariable(self.PREDICTOR,self.COVARIATE,assigned_variables,self.BIN,self.group_sizes,permute).create(s_parents,s_opts,s_lists) 



	def select_labels(self,ID,MINSIZE=5,MAXGROUPS=20):

		id_vals =  self.vals[ID]
		id_opts = list(set(self.vals[ID]))
		id_cc   =  sorted(cc(self.vals[ID]).items(), key = lambda X: X[1],reverse=True)

		p_opts = [c[0] for i,c in enumerate(id_cc) if (i==0 or (i<MAXGROUPS and c[1]>MINSIZE))]
	
		f_opts = [opt for opt in id_opts if opt not in p_opts] 



		if len(f_opts) > 0: 
			p_opts.append('UNAVAIL') 
			
	
		v_locs =  [p_opts.index(v) if v in p_opts else p_opts.index('UNAVAIL') for v in id_vals]

		

		return v_locs , p_opts 


	def __str__(self):
        	return str(self.names)
			
			

class RegVariable:
	def __init__(self,PREDICTOR,COVARIATE,vals,BIN,group_sizes,permute): 
		self.PREDICTOR,self.COVARIATE, self.vals, self.BIN, self.group_size, self.permuted, self.parent, self.children, self.notes =  PREDICTOR, COVARIATE, vals, BIN, group_sizes, permute,{} , dd(list), dd(list)  
		self.segregator = VariableSegregator(self.vals,self.BIN) 

		



	def create(self,parents,options,lists):

		self.parents, self.options, self.lists =  parents, options , lists
		self.data = [[a for b in x for a in b] for x in lists]
	 	self.array = np.array([np.array(x) for x in self.data])
		self.null  = np.array([[1.0] for x in self.data])
	


		for i in range(len(self.options)):
			for opt in self.options[i]:

 


				self.parent[opt] = self.parents[i] 
				self.children[self.parents[i]].append(opt) 
			

				

	



		self.names = [a for b in self.options for a in b] 
		self.len = len([n for n in self.names if n != 'intercept'])
		

		self.predictor_idx, self.covariate_idx = [] , [] 

		for i,n in enumerate(self.names): 

			self.PREDICTOR[n] = self.PREDICTOR[self.parent[n]] 
			self.COVARIATE[n] = self.COVARIATE[self.parent[n]] 
			if n != 'intercept' and self.PREDICTOR[n]:   self.predictor_idx.append(i) 
			elif n != 'intercept' and self.COVARIATE[n]: self.covariate_idx.append(i) 

		self.zp = np.array([[1.0] for n in range(len(self.array))])


		
		self.v_type = {} 
		self.i_type = {} 


		for i,n in enumerate(self.names): 

			if self.PREDICTOR[n]: self.i_type[n] = 'PREDICTOR'
			else: 		      self.i_type[n] = 'COVARIATE'

			if self.BIN[self.parent[n]]:   self.v_type[n] = 'BIN' 
			else: 			       self.v_type[n] = 'CONT'



		return self




	def isolate_covariate(self,covariate):

		s = copy.deepcopy(self)
		


	def __str__(self):
        	return str(self.names)







class VariableSegregator:
	def __init__(self,vals,BIN): 
		self.vals = vals 
		self.key, self.type, self.opts, self.lens,self.parent, self.fracs = {}, {}, {}, {} , {} , {} 
		for k,V in self.vals.items(): 



			if BIN[k]: 
				self.type[k] = 'binary'


				self.opts[k] = list(set(V)) 	



				self.key[k] =  {opt: [i for i in range(len(V)) if V[i] == opt] for opt in self.opts[k]}
				self.lens[k] = {opt: len(self.key[k][opt]) for opt in self.opts[k]}


				 
 
				for opt in self.opts[k]: 
					self.parent[opt] = k 


				self.fracs[k] = {opt: self.lens[k][opt] / float(len(V)) for opt in self.lens[k]} 
			else:
				self.parent[k] = k 
				self.type[k] = 'continuous'
				self.opts[k] = [min(V),np.mean(V),np.median(V),max(V)]
				ASSIGNED = dd(bool) 
				val_srt = sorted([(v,i) for i,v in enumerate(V)]) 
				ASSIGNED[val_srt[0][1]] = 'LO' 
				ASSIGNED[val_srt[-1][1]] = 'HI' 
				LO,HI = [val_srt[0][0]], [val_srt[-1][0]] 
				val_data = val_srt[0:-1] 
				i = 1
				while True: 
					vL,iL = val_data[i]  
					vH,iH = val_data[-1*i] 		
					if iL == iH: 
						if vH == HI[-1]: 
							HI.append(vH) 
							ASSIGNED[iH] = 'HI' 
						elif vL == LO[-1]: 
							LO.append(vL) 
							ASSIGNED[iL] = 'LO' 
						elif vL <= np.mean(V): 
							LO.append(vL) 
							ASSIGNED[iL] = 'LO' 
						else: 
							HI.append(vH) 
							ASSIGNED[iH] = 'HI' 
						break 
					elif ASSIGNED[iL] and ASSIGNED[iH]: break  
					else: 
	
						if vL < vH: 
							LO.append(vL) 
							HI.append(vH) 
							ASSIGNED[iL],ASSIGNED[iH] = 'LO','HI' 
						elif vH == HI[-1]: 
							HI.append(vH) 
							HI.append(vL) 
							ASSIGNED[iL], ASSIGNED[iH] = 'HI','HI' 
						elif vL == LO[-1]: 
							LO.append(vH) 
							LO.append(vL) 
							ASSIGNED[iL], ASSIGNED[iH] = 'LO','LO' 
					i+=1 
			

				self.key[k] = {'LO': [i for i  in ASSIGNED.keys() if ASSIGNED[i] == 'LO'],'HI': [i for i  in ASSIGNED.keys() if ASSIGNED[i] == 'HI']}



				#self.key[k]  = {'LO': [i for i in range(len(V)) if V[i] <= self.opts[k][2]], 'HI': [i for i in range(len(V)) if V[i] > self.opts[k][2]]}
				self.lens[k] = {opt: len(self.key[k][opt]) for opt in ['LO','HI']} 
				self.fracs[k] = {opt: self.lens[k][opt] / float(len(V)) for opt in self.lens[k]} 





	def score(self,p,Y):




		parent = self.parent[p] 

		vals, key, lens, opts, fracs = self.vals[parent], self.key[parent], self.lens[parent], self.opts[parent], self.fracs[parent] 
		yKey  = {k: [Y[i] for i in key[k]] for k in key} 
		if self.type[parent] == 'binary': 
			yP    =  yKey[p] 
			yElse = [a for b in [yKey[n] for n in yKey if n != p] for a in b] 
			yMeans = {k: np.mean(v) for k,v in yKey.items()} 
			yMeans['_ELSE_'] = np.mean(yElse) 
			yObs = {k: len([vi for vi in v if vi>0])/float(len(v)) for k,v in yKey.items()}
			yObs['_ELSE_'] = len([vi for vi in yElse if vi > 0]) / float(len(yElse))
			avgList = sorted(yMeans.items(),key = lambda X: X[1]) 	
			yMix =  sorted([(Y[i],vals[i]) for i in range(len(Y))],reverse=True)
			pCnts, eCnts = dd(int) , dd(int) 
			for i,(v,n) in enumerate(yMix): 
				if v == 0 or i == len(yElse): break 
				elif i <= lens[p]: pCnts[n] += 1
				eCnts[n] += 1 
			pList, eList   = [pCnts[k] for k in opts] , [eCnts[k] for k in opts] 
			pExp,  eExp    = [fracs[k]*sum(pList) for k in opts], [fracs[k]*sum(eList) for k in opts]
			pChi,  eChi    = chisquare(pList,f_exp=pExp)[1], chisquare(eList,f_exp=eExp)[1]
			if yMeans[p] < yMeans['_ELSE_']: 	FC = -1* (yMeans['_ELSE_']/(yMeans[p]+0.1))
			else:					FC = yMeans[p] / (yMeans['_ELSE_']+0.1)
			kObs = len([vi for vi in yKey[p]]) 
			return kObs,yObs[p], yObs['_ELSE_'],round(FC,3),np.mean([pChi, eChi]) 


		else:
			kObs = len(Y) 
			yKey = {YK: [Y[i] for i in self.key[p][YK]] for YK in self.key[p].keys()}
			yMeans = {k: np.mean(v) for k,v in yKey.items()} 
			yObs = 	 {k: len([vi for vi in v if vi>0])/float(len(v)) for k,v in yKey.items()}
			vLow, vTop =  [x[0] for x in sorted(yMeans.items(), key=lambda X: X[1])]
			topV, lowV, tLen = yMeans[vTop], yMeans[vLow], len(yKey[vTop])
			topO, lowO = yObs[vTop], yObs[vLow]
			FC = round(topV / (lowV+0.01),3) 

			hLen = len([x for x in sorted([a for b in [[(v,YK) for v in yKey[YK]]for YK in yKey.keys()] for a in b],reverse=True)[0:tLen] if x[1] == vTop])
			chsq = chisquare([hLen,tLen-hLen], f_exp = (tLen/2.0,tLen/2.0))[1] 

			return kObs, topO, lowO, FC, chsq

 

	def summarize(self,parent,Y): 


		self.summaries = {}
		vals, key, lens, opts, fracs = self.vals[parent], self.key[parent], self.lens[parent], self.opts[parent], self.fracs[parent] 


		yKey  = {k: [Y[i] for i in key[k]] for k in key} 
		yMeans = {k: np.mean(v) for k,v in yKey.items()} 
		yObs = {k: len([vi for vi in v if vi>0])/float(len(v)) for k,v in yKey.items()}

		yObsMean = sorted([(np.mean(v),len([vi for vi in v if vi>0])/float(len(v)),k) for k,v in yKey.items()])

		yMix = sorted([a for b in [[(v,k) for v in yKey[k]] for k in yKey.keys()] for a in b])
		yZero = [ym for ym in yMix if ym[0] == 0] 
		yNZ   = [ym for ym in yMix if ym[0] > 0]  
		yRanks = yZero + [(i+1,yNZ[i][1]) for i in range(len(yNZ))]
		yChis = {p : [yRanks[i][0] for i in range(len(yRanks)) if yRanks[i][1] == p] for p in yKey.keys()}
		chiSort = sorted([(np.mean(yChis[p]),p) for p in yChis.keys()])
		low_len, hi_len =  max(lens[chiSort[0][1]],len(yZero)), lens[chiSort[-1][1]]
		lowExp,  hiExp    = [fracs[k]*low_len for k in yChis.keys()], [fracs[k]*hi_len for k in yChis.keys()]
		lowCC = cc([yRanks[i][1] for i in range(0,low_len)])
		hiCC  =  cc([yRanks[i][1] for i in range(len(yRanks)-hi_len,len(yRanks))])
		lowObs,hiObs = [lowCC[k] for k in yChis.keys()], [hiCC[k] for k in yChis.keys()] 
		lowChi,  hiChi    = chisquare(lowObs,f_exp=lowExp)[1], chisquare(hiObs,f_exp=hiExp)[1]






		if self.type[parent] == 'binary':  

			VS = VariableSummarize('binary').add_data(yObsMean,lowChi,hiChi)
			#VS = VariableSummary('binary') #.add_data(yObsMean,chiSort,(lowChi,hiChi))

			

			#VS.add_data(sorted(yMeans.items(),key=lambda X:X[1]),sorted(yObs.items(),key=lambda X: X[1]),chiSort,(lowChi,hiChi))
		else:
			VS = VariableSummarize('contiinuous').add_data(yObsMean,lowChi,hiChi)
			#VS = VariableSummary('continuous')
			#VS.add_data(sorted(yMeans.items(),key=lambda X:X[1]),sorted(yObs.items(),key=lambda X: X[1]),chiSort,(lowChi,hiChi))


		return VS  












class VariableSummarize:
	def __init__(self,variable_type):

		self.type = variable_type 


	def add_data(self,data,lowChi,hiChi):

		self.chiLo, self.chiHi = lowChi,hiChi

		avgs,obs,self.names = [x[0] for x in data],[x[1] for x in data],[x[2] for x in data]

		self.fc = data[-1][0] / float(data[0][0]) 
		

		self.out = {'names': ",".join([x[2].split('~')[-1] for x in data])+',', 'obs': ','+','.join([str(int(100*round(a[1],2))) for a in data])+',', 'avg':  ','+','.join([str(round(a[0],2)) for a in data])+','}

		return self 




class VariableSummary:
	def __init__(self,variable_type):

		self.type = variable_type 


	def add_data(self,avgs,obs,chis,(lowChi,hiChi)):

		avg_str = ','+','.join([a[0].split('~')[-1] for a in avgs])+','
		avg_vals = ','+','.join([str(round(a[1],3)) for a in avgs])+','
		obs_str = ','+','.join([a[0].split('~')[-1] for a in obs])+','
		obs_vals = ','+','.join([str(round(a[1],3))[1::] if a[1]<1 else str(1) for a in obs])+','
		obs_cov = ','+','.join([str(1) for a in obs])+','
		chi_str = ','+','.join([a[1].split('~')[-1] for a in chis])+','
		chi_vals = ','+','.join([str(round(a[0])) for a in chis])+','


		FC =str(round(avgs[-1][1] / np.mean([a[1] for a in avgs][0:-1]),2))+' '+str(round(-1*(np.mean([a[1] for a in avgs][1::]) / (0.01+avgs[0][1])),2))

		self.vals = {'avg': avg_vals, 'obs': obs_vals, 'chi': chi_vals,'rbs': obs_cov} 
		self.strs = {'avg': avg_str, 'obs': obs_str, 'chi': chi_str} 
		self.FC = FC
		self.chis = lowChi, hiChi 













