import sys
from collections import defaultdict as dd
from collections import Counter as cc
from collections import MutableSequence
import numpy as np
from math import log
import mmap 
from sklearn.preprocessing import MinMaxScaler
from random import shuffle
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib.lines import Line2D as Line
from matplotlib.patches import Circle, Wedge, Polygon
import copy 
import rage_members 
import rage_variables







def command_line_error(msg):
	sys.stderr.write('\n')
        sys.stderr.write('RageCommandLineError: '+msg+'\n')
        sys.exit()


class cnt_key(dd):
    	def __missing__(self, key):
        	return 0.0

	def update(self,mylist):
		for (a,b) in mylist: self[a] = b 
		return self



MIN_PRED_SIZE=10
MIN_GRP_SIZE=20
MAX_GRP_MEMBERS=30
MAX_GRP_MEMBERS=20














class read: 
	def __init__(self,rage,data_idx=0,data_type='raw'):

		self.rage = rage
		self.options = rage.args  
		rage.progress.start_major('Reading Input') 	

		if data_type == 'CONDENSED': 
			self.cnt_file = rage.args.condensedcnts[data_idx].name 
			self.read_primary_counts(rage.args.condensedcnts[data_idx]) 
			self.add_name_attributes() 
			self.add_summary_attributes() 
			#if not rage.args.sampleKey:	self.add_name_attributes() 
		else:

			self.cnt_file = rage.args.counts[data_idx].name
			self.read_primary_counts(rage.args.counts[data_idx]) 
			if 'sampleKey' in vars(rage.args) and rage.args.sampleKey != None: 
				self.add_sample_attributes(rage.args.sampleKey)


			self.add_name_attributes() 
			self.add_summary_attributes() 
		
			if 'featureKey' in vars(rage.args) and rage.args.featureKey != None: self.add_feature_attributes(rage.args.featureKey)
			self.add_summary_attributes(FEATURES=True) 



		self.TRANSFORMED = None 
		self.SCALING     = None 



	
	def read_primary_counts(self,cnt_file,msg='Reading Primary Counts'): 
		self.rage.progress.start_minor(msg,5000,True) 	


		f_line =  cnt_file.readline().split()[1::]


		self.features, kv, kHi,kMid,kLo = rage_members.Members(), 0 , 0 , 0 , 0 
		if self.rage.args.test > 0:
			self.test_idx = [i for i in range(1,len(f_line)) if i % self.rage.args.test == 0]
			self.samples = rage_members.Members([f_line[i] for i in self.test_idx]) 
			MIN_SAMPLES =  max(3,int(0.999 + (len(self.samples) * self.options.min_obs_rate)))
			#MIN_SAMPLES = 0 
			

			for j,line in enumerate(cnt_file): 
				line = line.split()
 				fr,test_cnts,PASS = line[0],[float(line[i]) for i in self.test_idx], True
#				fr,cnts = line[0],[float(line[i]) for i in self.test_idx in line[1::]]
#				test_cnts = [cnts[i] for i in self.test_idx] 
				obs, obs_rate, nums = len([x for x in test_cnts if x > 0]), len([x for x in test_cnts if x > 0]) /float(len(test_cnts)), len(set(test_cnts))
				if obs < MIN_SAMPLES: PASS = False
				elif kv > 10 and (obs_rate < 0.25 or nums < 5): PASS = False 
			

				if PASS: 
					self.features.append(rage_members.Member(fr,kv).add_line_cnts(test_cnts))
					for i,c in enumerate(test_cnts): 
						if c > 0: self.samples[i].cnts[kv] = c
					kv+=1


				if kv > self.samples.len: break 


		else: 

			MIN_SAMPLES =  max(2,int(0.999 + (len(f_line) * self.options.min_obs_rate)))

			self.samples = rage_members.Members(f_line) 

			for j,line in enumerate(cnt_file):

 
				self.rage.progress.mark() 
				line = line.split() 
				

				fr,cnts,PASS = line[0],[float(x) for x in line[1::]],True
			#	if fr in ["ENSG00000251562;chr11;MALAT1","ENSG00000281344;chr12;HELLPAR"]: continue 

				obs = len([x for x in cnts if x > 0]) 
				nums = len(set(cnts)) 


				if obs < MIN_SAMPLES: PASS = False 
				elif nums < 4: 	      PASS = False 
				# PASS = False 
				# if obs < 2: continue 
				# if len(self.features) < 5: PASS=True 
				# elif (nums>3) and (obs > 10) and (obs/float(len(cnts)) > self.rage.args.min_obs_rate): PASS=True 
				


				if PASS: 
					self.features.append(rage_members.Member(fr,kv).add_line_cnts(cnts))
					

					for i,c in enumerate(cnts): 
						if c > 0: self.samples[i].cnts[kv] = c
					kv+=1

		
		self.samples.collate('samples') 
		self.features.collate('features')
		self.rage.progress.end() 



	def add_name_attributes(self):



		try: 
			cands =  list(set([s.name.split('~')[1] for s in self.samples]))


			if len(cands) == 1: 
				cand = cands[0] 

				if cand not in self.samples.attributes: 
					self.samples.add_attribute(cand,'binary') 
					for i,s in enumerate(self.samples): 
						s.attributes[cand] = cand+'~'+s.name.split('~')[-1] 




					return 


			elif len(cands) < 10:

 
				vals =  [s.name.split('~')[-1] for s in self.samples]
				self.samples.add_attribute('C1','binary') 
				for i,s in enumerate(self.samples): 
					s.attributes['C1'] = s.name.split('~')[-1] 
				return 	

			else:
				return 

		except IndexError: 
			return 				




	def add_sample_attributes(self,sampleKey):

		self.rage.progress.start_minor('Reading Sample Annotation',1000,True) 	
		samples, opt_key   =  [], [(h,[]) for h in sampleKey.readline().split()[1::]]
		for line in sampleKey: 
			self.rage.progress.mark() 
			line = line.split()
			samples.append(line[0]) 
			for i,v in enumerate(line[1::]): opt_key[i][1].append(v) 
		for opt,vals in opt_key:
			try:  
				for s_name,s_val in zip(samples,['NA' if p == 'NA' else float(p) for p in vals]):
					if s_name in self.samples.lookup: self.samples[self.samples.lookup[s_name]].attributes[opt] = s_val 	
				self.samples.add_attribute(opt,'continuous') 
			except ValueError: 
				for s_name,s_val in zip(samples,vals):	
					if s_name in self.samples.lookup: self.samples[self.samples.lookup[s_name]].attributes[opt] = opt+'~'+s_val 	
				self.samples.add_attribute(opt,'binary') 
		sampleKey.seek(0)


	def add_summary_attributes(self,FEATURES=False): 

		if FEATURES: S = self.features
		else:        S = self.samples


		if 'TOTAL' not in S.attributes: 
			S.add_attribute('TOTAL','continuous') 	
			for s in S: 	s.attributes['TOTAL'] = sum(s.cnts.values()) 

		if 'OBS' not in S.attributes: 

			S.add_attribute('OBS','continuous') 	
			for s in S: 	s.attributes['OBS'] = len(s.cnts) 
		


		
			

	def filter_samples_by_attributes(self,predictors=[],covariates=[],keep=False):

		predictors = list(set([a for b in [p.split('*') for p in predictors] for a in b]))
		covariates = list(set([a for b in [p.split('*') for p in covariates] for a in b]))
		variables  = list(set(predictors + covariates))

		missing_attributes = [v for v in variables if v.split('=')[0] not in self.samples.attributes]



		if len(missing_attributes) > 0: 
			for m in missing_attributes: 
				if m[0:3].upper() == 'OBS': 
					self.samples.add_attribute(m,'continuous') 	
					for s in self.samples: 	s.attributes[m] = len(s.cnts) 
				elif m[0:3].upper() == 'TOT': 
					self.samples.add_attribute(m,'continuous') 	
					for s in self.samples: 	s.attributes[m] = sum(s.cnts.values()) 
				else:	
					command_line_error('Supplied variable '+m+' not included in key files')



		if self.rage.args.allowNA:   	self.samples.filter(predictors,self.rage.args.prune) 
		else: 				self.samples.filter(variables,self.rage.args.prune) 


		

		for i in range(len(self.features)):	
			new_key = cnt_key(None) 
			for a,b in self.features[i].cnts.items(): 
				if a in self.samples.swap: new_key[self.samples.swap[a]] = b 
			self.features[i].cnts = new_key
		return self 






	def matrix(self,data=None,log_transform = False,center=True,scaling=None,TRANSPOSE=False):


		if data == None: sample_cnts = [[0 if f.idx not in s.cnts else s.cnts[f.idx] for f in self.features] for s in self.samples]
		else:		 sample_cnts = [[data[i][j] for i in range(len(data))] for j in range(len(data[0]))]




		if log_transform: 	dMatrix = np.matrix([[log(1.0+x,2) for x in sample] for sample in sample_cnts])
		else: 		        dMatrix = np.matrix(sample_cnts)
#		if transform == 'log': 	dMatrix =  np.matrix([[0 if f.idx not in s.cnts else log(1.0+s.cnts[f.idx],2) for f in self.features] for s in self.samples])
#		else:			dMatrix =  np.matrix([[0 if f.idx not in s.cnts else s.cnts[f.idx] for f in self.features] for s in self.samples])


		
		if TRANSPOSE: 
			dMatrix = dMatrix.getT() 

		
		if center:	return dMatrix - dMatrix.mean(axis=0)
		else: 		return dMatrix	




	def scale_and_transform(self,LOG=False,SCALE=False): 
	
		

		min_obs = min(100,int(len(self.samples)*0.05))
		min_obs = 1
		f_swap = {}

		for i in range(len(self.features)): 
			f_swap[self.features[i].idx] = i 
			self.features[i].idx = i
			if LOG: 	  self.features[i].cnts = cnt_key(None).update([(j,log(1.0+c,2)) for j,c in self.features[i].cnts.items()]) 
			else: 		  self.features[i].cnts = cnt_key(None).update([(j,c) for j,c in self.features[i].cnts.items()]) 
 


		for i,s in enumerate(self.samples):	self.samples[i].cnts = cnt_key(None).update([(f_swap[a],self.features[f_swap[a]].cnts[i]) for (a,b) in self.samples[i].cnts.items() if a in f_swap])
		

		if LOG:  self.TRANSFORMED = 'log' 

		return self


#	def normalize(self,log_transform = True):
#
#
#		if self.NORMALIZATION != None: 
#			print 'huh'
#			sys.exit() 
#
#		min_obs = min(100,int(len(self.samples)*0.05))
#		min_obs = 0 
#
#
#		f_swap = {}
#		self.features = self.features.prune_cnts(min_obs) 
#
#		for i in range(len(self.features)):
#			f_swap[self.features[i].idx] = i 
#			self.features[i].idx = i
#			if log_transform: self.features[i].cnts = cnt_key(None).update([(j,log(1.0+c,2)) for j,c in self.features[i].cnts.items()]) 
#			else: 		  self.features[i].cnts = cnt_key(None).update([(j,c) for j,c in self.features[i].cnts.items()]) 
#
#		for i,s in enumerate(self.samples):
#			self.samples[i].cnts = cnt_key(None).update([(f_swap[a],self.features[f_swap[a]].cnts[i]) for (a,b) in self.samples[i].cnts.items() if a in f_swap])
#
#
#		self.NORMALIZATION = True 
#
#		return self	



	def set_sample_variables(self,predictors, covariates = []):


		self.variable_options, self.variable_key, self.value_key, reg_predictors, reg_covariates = {}, {}, {}, [], []  


		MINSIZE   = self.options.mingroupsize
		MAXGROUPS = self.options.maxgroups




		for V in predictors+covariates: 



			self.rage.progress.start_minor('Setting Sample Values For: '+V+'..',False)	
			s_vars  = [[s.attributes[v.split('=')[0]] for v in V.split('*')] for s in self.samples]
			s_types = [self.samples.attribute_class[v.split('=')[0]] for v in V.split('*')] 



			modV  = "*".join([v.split('=')[0] for v in V.split('*')])
			if V in predictors: reg_predictors.append(modV) 
			else:		    reg_covariates.append(modV) 
			V = modV 




			if list(set(s_types)) == ['binary']: 


				

				if len(s_types) == 1: 
					self.value_key[V]  =  ["_*_".join(sv) for sv in s_vars]

				else:

					self.value_key[V] =  ["*".join([sx.split('~')[0] for sx in sv])+'~'+"*".join([sx.split('~')[-1] for sx in sv]) for sv in s_vars]





				opts = list(set(self.value_key[V]))
				cnts = sorted(cc(self.value_key[V]).items(),key=lambda X: X[1], reverse=True)

				p_opts = [c[0] for i,c in enumerate(cnts) if (i==0 or (i<MAXGROUPS and c[1]>MINSIZE))]
				f_opts =  [opt for opt in opts if opt not in p_opts] + [p_opts.pop() for x in range(1) if len(p_opts) == len(opts)] 

				if len(p_opts) == 0: 
					p_opts = [f_opts[0]]
					f_opts = f_opts[1::] 

				self.variable_options[V] = [p_opts, f_opts] 
				self.variable_key[V] = {val: [1 if p_opts[i] == val else 0 for i in range(len(p_opts))] for val in list(set(self.value_key[V]))}


				 
			elif list(set(s_types)) == ['continuous']:

				valid_mean = np.mean([np.prod(sv) for sv in s_vars if 'NA' not in sv]) 
				self.value_key[V] = [np.prod(sv) if 'NA' not in sv else valid_mean for sv in s_vars]
				self.variable_options[V] = [[V],[self.value_key[V],np.mean(self.value_key[V])] ]
				self.variable_key[V] = {sv: [sv] for sv in self.value_key[V]}
				
			elif s_types == ['continuous','binary'] or s_types == ['binary','continuous']: 

				b_mod = {x: -0.5+i if i <2 else 0 for i,x in enumerate([c[0] for c in sorted(cc([sv[s_types.index('binary')] for sv in s_vars]).items(), key = lambda X: X[1],reverse=True)])}
				self.value_key[V] = [np.prod([b_mod[x] if s_types.index('binary') == i else x for i,x in enumerate(sv)]) for sv in s_vars]
				self.variable_options[V] = [[V],[self.value_key[V],np.mean(self.value_key[V])] ]
				self.variable_key[V] = {sv: [sv] for sv in self.value_key[V]}




			if len(s_types) > 1: 

				if 'continuous' in s_types: self.samples.add_attribute(V,'continuous')
				else: 			    self.samples.add_attribute(V,'binary')
				for i,s in enumerate(self.samples): s.attributes[V] = self.value_key[V][i] 
	


		


		self.rage.progress.end()

		return  rage_variables.RegVariables(self.value_key,self.variable_options,self.variable_key,reg_predictors,reg_covariates)
		









		
