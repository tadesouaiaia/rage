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


def get_color_list(list_size,MISSING,rank=None):

	basics  = ['blue','green', 'red',     'black',   'grey', 'cornflowerblue'] 
	brights = ['lime','cyan',  'magenta', 'purple',  'orange','yellow']
	pastels = ['gold','peru',  'salmon',  'olive',    'khaki','indigo','crimson']



	if rank == None: 
		if len(list_size) <= 10: color_list = basics[0:3]+brights[0:5]+pastels[0:5]
		else: 			 color_list = basics+brights+pastels+basics+brights+pastels
	else: 
		if rank[0] < rank[1]: 	color_list = basics+pastels+basics+basics
		else:			color_list = brights+pastels[::-1]+brights+brights

	if MISSING: return [c for c in color_list[0:list_size] if c != 'black'],'black'
	else: 	    return color_list[0:list_size],None





def get_marker_list(list_size,MISSING,rank=None): 


	marks = ["o","v","^","<",">","8","s","p","H","D","d"]
	
	if MISSING: return [m for m in marks if m != "s"][0:list_size],"s"
	else:       return marks[0:list_size], None 















def predictor_warning(msg):
        sys.stderr.write('RagePredictorWarning: '+msg+'\n')



def predictor_error(msg):
	sys.stderr.write('\n')
        sys.stderr.write('RagePredictorError: '+msg+'\n')
        sys.exit()


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


class attribute_key(dd):
    	def __missing__(self, key):
        #value = self.default_factory(key)
        	self[key] = 'NA'
        	return 'NA'


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
			if not rage.args.sampleKey:	self.add_name_attributes() 
		else:

			self.cnt_file = rage.args.counts[data_idx].name
			self.read_primary_counts(rage.args.counts[data_idx]) 
			if 'sampleKey' in vars(rage.args) and rage.args.sampleKey != None: self.add_sample_attributes(rage.args.sampleKey)
			self.add_name_attributes() 
		
		self.predictor_minsize = 15
		self.predictors, self.covariates =[],[]


		self.cutoffs = dd(int) 



#		if 'predictors' in self.rage.args: self.predictors = self.rage.args.predictors
#		if 'covariates' in self.rage.args: self.covariates = self.rage.args.covariates
#		self.variables = self.predictors + self.covariates 
	
	def read_primary_counts(self,cnt_file,msg='Reading Primary Counts'): 
		self.rage.progress.start_minor(msg,5000,True) 	


		f_line =  cnt_file.readline().split()[1::]


		#f_line =  cnt_file.readline().split()[1::]
		self.features, kv, kHi,kMid,kLo = rage_members.Members(), 0 , 0 , 0 , 0 
		if self.rage.args.test > 0:
			self.test_idx = [i for i in range(len(f_line)) if i % self.rage.args.test == 0]
			self.samples = rage_members.Members([f_line[i] for i in self.test_idx]) 

			
			for j,line in enumerate(cnt_file): 
				line = line.split()
 
				fr,cnts = line[0],[float(x) for x in line[1::]]
				test_cnts = [cnts[i] for i in self.test_idx] 
				obs_rate = len([x for x in test_cnts if x >0]) / float(self.samples.len) 
				kPass = False
				if obs_rate > 0.5: kPass = True
				elif obs_rate > 0.25 	and  len(self.features) < self.samples.len * 2: kPass=True
				elif obs_rate > 0.10    and  len(self.features) < self.samples.len:     kPass = True 
				elif obs_rate > 0  	and  len(self.features)  < 100: 		kPass = True

				if kPass: 

					self.features.append(rage_members.Member(fr,kv).add_line_cnts(test_cnts))
					for i,c in enumerate(test_cnts): 
						if c > 0: self.samples[i].cnts[kv] = c
					kv+=1
				if kv > self.samples.len*3: break 


		else: 
			self.samples = rage_members.Members(f_line) 
			for j,line in enumerate(cnt_file): 
				self.rage.progress.mark() 
				line = line.split() 
				fr,cnts,PASS = line[0],[float(x) for x in line[1::]],False
				obs = len([x for x in cnts if x > 0]) 
				nums = len(set(cnts)) 
				if obs < 2: continue 
				if len(self.features) < 5: PASS=True 
				elif (nums>3) and (obs > 10) and (obs/float(len(cnts)) > self.rage.args.min_obs_rate): PASS=True 
				if PASS: 
					self.features.append(Member(fr,kv).add_line_cnts(cnts))
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
						s.attributes[cand] = s.name.split('~')[-1] 
					return 
		
			elif len(cands) < 10: 
				vals =  [s.name.split('~')[-1] for s in self.samples]
				self.samples.add_attribute('C1','binary') 
				for i,s in enumerate(self.samples): 
					s.attributes['C1'] = s.name.split('~')[-1] 
				return 	

			else:
				return 

		except IndexError: return 				

















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

	def matrix(self,transform = None,center=True,scaling=None):
		if transform == 'log': 	dMatrix =  np.matrix([[0 if f.idx not in s.cnts else log(1.0+s.cnts[f.idx],2) for f in self.features] for s in self.samples])
		else:			dMatrix =  np.matrix([[0 if f.idx not in s.cnts else s.cnts[f.idx] for f in self.features] for s in self.samples])

		if center:	return dMatrix - dMatrix.mean(axis=0)
		else: 		return dMatrix	

			

	def filter_samples_by_attributes(self,predictors=[],covariates=[],keep=False):


		#self.predictors = list(set([a for b in [p.split('*') for p in predictors] for a in b]))
		#self.covariates = list(set([a for b in [p.split('*') for p in covariates] for a in b]))
		#self.variables  = list(set(self.predictors + self.covariates))

		
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


	def normalize(self,log_transform = True):
		min_obs = min(100,int(len(self.samples)*0.05))
		min_obs = 0 

		f_swap = {}
		self.features = self.features.prune_cnts(min_obs) 

		for i in range(len(self.features)):
			f_swap[self.features[i].idx] = i 
			self.features[i].idx = i
			if log_transform: self.features[i].cnts = cnt_key(None).update([(j,log(1.0+c,2)) for j,c in self.features[i].cnts.items()]) 
			else: 		  self.features[i].cnts = cnt_key(None).update([(j,c) for j,c in self.features[i].cnts.items()]) 

		for i,s in enumerate(self.samples):
			self.samples[i].cnts = cnt_key(None).update([(f_swap[a],self.features[f_swap[a]].cnts[i]) for (a,b) in self.samples[i].cnts.items() if a in f_swap])

		return self	








	def set_sample_values(self,min_size=MIN_GRP_SIZE,min_pred_size=MIN_PRED_SIZE,maxgroup=MAX_GRP_MEMBERS):

		
		

		S, s_preds, s_covars = self.samples, dd(list), dd(lambda: {}) 
		self.variable_options = dd(list) 
	
		for a in self.variables: 
			self.rage.progress.start_minor('Setting Sample Values For: '+a+'..',False)
			if a.split('=')[0] not in self.samples.attribute_class: predictor_error('Uknown variable '+a)
			cdic = cc([s.attributes[a] for s in self.samples]) 

			

			vals = [s.attributes[a] for s in self.samples if str(s.attributes[a]).split('~')[-1] != 'NA'] 
			cnts = sorted(cc(vals).items(),key=lambda G: G[1])
			opts = list(set([s.attributes[a] for s in self.samples]))

			if self.samples.attribute_class[a.split('=')[0]] == 'binary': 


				if a in self.predictors: 
					pred,p_opts = a,[x[0] for x in cnts if x[1] > 2 or x[1] == cnts[-1][1]]
				
					if len(p_opts) > 2: 		   p_opts = [x[0] for x in cnts if x[1] > min_pred_size/2.0 or x[1] == cnts[-1][1]]					
					if len(p_opts) > min(10,maxgroup): p_opts = [x[0] for x in cnts if x[1] > min_pred_size or x[1] == cnts[-1][1]]					
				


					self.variable_options[a] = [p_opts,[opt for opt in opts if not opt in p_opts]]
					for i,s in enumerate(self.samples): 
						s_preds[i].append([s.attributes[pred],s.attributes[pred] in p_opts])
				else:
					#pred,p_opts = a,[x[0] for x in cnts if x[1] > min_size or x[1] == cnts[-1][1]]
					pred,p_opts = a,[x[0] for x in cnts if x[1] > 2 or x[1] == cnts[-1][1]]


					if len(p_opts) > min(10,maxgroup): p_opts = [x[0] for x in cnts if x[1] > min_pred_size or x[1] == cnts[-1][1]]					
#					if len(p_opts) > maxgroup:  	   p_opts = [x[0] for x in cnts if x[1] > min_size or x[1] == cnts[-1][1]]					
					if len(p_opts) > maxgroup: p_opts =  p_opts[len(p_opts)-maxgroup::]
					self.variable_options[a] = [p_opts,[opt for opt in opts if not opt in p_opts]]
					for i,s in enumerate(self.samples): 
						s_covars[pred][i] = (s.attributes[pred],s.attributes[pred] in p_opts)	

				if len(self.variable_options[a][1]) > 0:
					missed = self.variable_options[a][1]+['NA'] 
					self.rage.progress.notate('..Options: '+",".join([m.split('~')[-1] for m in missed])+' are combined due to low numbers ('+",".join([str(cdic[m]) for m in  missed])+')...')

									

				continue 
			else:				

				self.variable_options[a] = [[a],[vals,np.mean(vals)]]
				pred = a.split('=')[0] 
				for i,s in enumerate(self.samples): 
					if a in self.predictors: 	s_preds[i].append([s.attributes[pred],True])
					else:				s_covars[pred][i] = (s.attributes[pred],s.attributes[pred] != 'NA')


		self.sample_predictors, self.sample_covariates = s_preds, s_covars

	
	
	def set_individual_predictors(self): 

		for v in [v for v in self.variable_options if self.variable_options[v][1] == []]:


			if  self.variable_options[v][0][0].split('~')[-1] == 'YES': self.variable_options[v][0] = self.variable_options[v][0][1::]+[self.variable_options[v][0][0]]


			self.variable_options[v] = [self.variable_options[v][0][1::],[self.variable_options[v][0][0]]] 
		self.seg = {} 
		self.variable_key = {} 
		for v,(p,f) in self.variable_options.items():
			if p != [v]:
				self.variable_key[v] =  {a:b for (a,b) in [(p[i],[1 if j==i else 0 for j in range(len(p))]) for i in range(len(p))]+[(k,[0 for n in p]) for k in f]}
				self.seg[v] = {p[j]: [i for i in range(len(self.samples)) if self.samples[i].attributes[v] == p[j]] for j in range(len(p))}
			else:
				self.variable_key[v] = {a:[b] for (a,b) in [(a,a) for a in f[0]]+[('NA',f[1])]}
				self.seg[v] = {'HI': [i for i,s in enumerate(self.samples) if (s.attributes[v] > f[1] and s.attributes[v] != 'NA')]}
				self.seg[v]['LO'] = 	[i for i,s in enumerate(self.samples) if (s.attributes[v] < f[1] and s.attributes[v] != 'NA')]


	

	#def set_sample_variables(self,model='OLS',combine=False):
	def set_sample_variables2(self,predictors, covariates = [], model='OLS',combine=False):



		self.predictors, self.covariates, self.variables = predictors, covariates, predictors + covariates



		self.set_sample_values() 
		self.set_individual_predictors() 



		print self.variable_options
		print "" 
		print self.variable_key 
		print ""
		print self.predictors 
		print self.covariates 

		# sample variable dictionary 
		# variable options dict - double list - pass/fail if bin else match, range 
		# print self.variable_key - [0 0 1] coding for each option [ 0 0 if it is a fail ] and match if continuous - list len 1 
		# list of preds/covars  
		
		return  RegVariables(self.samples.get_attributes(self.variables),self.variable_options,self.variable_key,self.predictors,self.covariates)
		 












	def set_sample_variables(self,predictors, covariates = [],MINSIZE=2,MAXGROUPS=5, model='OLS'):


		self.predictors, self.covarates, self.variables = predictors, covariates, predictors + covariates


		self.variable_options = {} 
		self.variable_key = {} 
		self.value_key = {} 
		self.sample_variables = dd(lambda: dd(list)) 

		S, s_preds, s_covars = self.samples, dd(list), dd(lambda: {}) 

		for V in predictors+covariates: 


			PREDICTOR = (V in predictors) 

			self.rage.progress.start_minor('Setting Sample Values For: '+V+'..',False)
			
			 
					

			s_vars =  [[s.attributes[v.split('=')[0]] for v in V.split('*')] for s in self.samples] 
			s_types = [self.samples.attribute_class[v.split('=')[0]] for v in V.split('*')] 

			if list(set(s_types)) == ['binary']: 
				self.value_key[V], opts, cnts  = ["_*_".join(sv) for sv in s_vars], list(set(["_*_".join(sv) for sv in s_vars])) ,sorted(cc(["_*_".join(sv) for sv in s_vars]).items(),key=lambda X: X[1], reverse=True)  
				p_opts = sorted([c[0] for i,c in enumerate(cnts) if (i==0 or (i<MAXGROUPS and c[1]>MINSIZE))],reverse=True) 
				f_opts =  [opt for opt in opts if opt not in p_opts] + [p_opts.pop() for x in range(1) if len(p_opts) == len(opts)] 
				self.variable_options[V] = [p_opts, f_opts] 
				self.variable_key[V] = {val: [1 if p_opts[i] == val else 0 for i in range(len(p_opts))] for val in list(set(self.value_key[V]))}

				 
			elif list(set(s_types)) == ['continuous']:

				self.value_key[V] = [np.prod(sv) for sv in s_vars]
				self.variable_options[V] = [[V],[self.value_key[V],np.mean(self.value_key[V])] ]
				self.variable_key[V] = {sv: [sv] for sv in self.value_key[V]}

			elif s_types == ['continuous','binary'] or s_types == ['binary','continuous']: 

				b_mod = {x: -0.5+i if i <2 else 0 for i,x in enumerate([c[0] for c in sorted(cc([sv[s_types.index('binary')] for sv in s_vars]).items(), key = lambda X: X[1],reverse=True)])}
				self.value_key[V] = [np.prod([b_mod[x] if s_types.index('binary') == i else x for i,x in enumerate(sv)]) for sv in s_vars]
				self.variable_options[V] = [[V],[self.value_key[V],np.mean(self.value_key[V])] ]
				self.variable_key[V] = {sv: [sv] for sv in self.value_key[V]}

		return  RegVariables(self.value_key,self.variable_options,self.variable_key,predictors,covariates)
		

















































			
		


class RegVariables:
	def __init__(self,sample_attributes,variable_options,variable_key,predictors,covariates):


		self.predictors, self.covariates = predictors,covariates
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
		shuffle_variables = [] 
		s_parents, s_opts, s_lists = [],[],[[] for s in range(len(self.sample_vals))]	
		for j,v in enumerate(self.names):
			if v == 'intercept' or v in variables: 
				if v in permute: 
					shuffle_variables.append((j,v)) 
				else:
					s_parents.append(v) 
					s_opts.append(self.options[v]) 
					for i in range(len(self.sample_vals)): s_lists[i].append(self.key[v][self.sample_vals[i][j]])
		

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
			for i in range(len(self.sample_vals)): s_lists[i].append(self.key[v][sv[i]])

		return RegVariable(self.PREDICTOR,self.COVARIATE,permute).create(s_parents,s_opts,s_lists) 




	def __str__(self):
        	return str(self.names)
			
			

class RegVariable:
	def __init__(self,PREDICTOR,COVARIATE,permute): 
		self.PREDICTOR,self.COVARIATE, self.permuted, self.parent, self.children =  PREDICTOR, COVARIATE, permute,{} , dd(list) 


	def create(self,parents,options,lists):

		self.parents, self.options, self.lists =  parents, options , lists
		self.data = [[a for b in x for a in b] for x in lists]
	 	self.array = np.array([np.array(x) for x in self.data])

		for i in range(len(self.options)):
			for opt in self.options[i]: 
				self.parent[opt] = self.parents[i] 
				self.children[self.parents[i]].append(opt) 
			

		self.names = [a for b in self.options for a in b] 

		

		self.predictor_idx, self.covariate_idx = [] , [] 

		for i,n in enumerate(self.names): 

			self.PREDICTOR[n] = self.PREDICTOR[self.parent[n]] 
			self.COVARIATE[n] = self.COVARIATE[self.parent[n]] 
			if n != 'intercept' and self.PREDICTOR[n]:   self.predictor_idx.append(i) 
			elif n != 'intercept' and self.COVARIATE[n]: self.covariate_idx.append(i) 

		return self

	def isolate_covariate(self,covariate):

		s = copy.deepcopy(self)
		


	def __str__(self):
        	return str(self.names)







		
		
			
class Member2:
	def __init__(self,name,idx):
		self.name = name 
		self.idx  = idx 
		self.cnts = cnt_key(None) 
		self.attributes = attribute_key(None)
#		self.cnts = ddd(None) 
#		self.attributes = ddd(None)
		self.preds =  dd(float)  
		self.regVariables = dd(list) 
		self.label = None 
		self.notes = dd(bool)
		self.null = False  
		self.hide = False
		self.norms = dd(lambda: cnt_key(None)) 

 
	def update(self):
		self.cnt_total = sum(self.cnts.values()) 	
		self.len = len(self.cnts) 
		return self.name 


	def add_line_cnts(self,cnts):
		for idx in [i for i in range(len(cnts)) if cnts[i] > 0]: self.cnts[idx] = cnts[idx]  
		self.cnt_total = sum(self.cnts.values()) 	
		self.len = len(self.cnts) 
		return self






class Members2(MutableSequence):

	"""A container for manipulating lists of hosts"""
	def __init__(self, data=None):
		"""Initialize the class"""
		super(Members, self).__init__()
		if (data is not None):
			self._list = [Member(data[i],i) for i in range(len(data))]
			self.len = len(self._list) 
			self.removed, self.attributes,self.attribute_type,self.attribute_class = [], [], [] , {} 
		else:
	    		self._list = list()
			self.len = 0 
			self.removed, self.attributes,self.attribute_type,self.attribute_class = [], [], [] , {} 

		self.keys = {} 
		self.notes = {} 



	def collate(self,label): 
		self.label = label
		self.names = [x.update() for x in self._list]
		self.lookup = {member.name : member.idx for member in self._list}
		self.len = len(self._list)
#		print len(self._list)
#		print len(self.__list) 
#		self.__list = self.__list[0:20] 


	def add_attribute(self,opt,optType):
		self.attributes.append(opt) 
		self.attribute_type.append(optType) 
		self.attribute_class[opt] = optType

	def get_attributes(self,attributes):

		s_key = dd(lambda: {}) 
		self.key  = {} 
		for a in attributes: 
			s_key[a] = [s.attributes[a] for i,s in enumerate(self._list)]
		return s_key


	def __repr__(self):
		return "<{0} {1}>".format(self.__class__.__name__, self._list)

	def __len__(self):
		"""List length"""
		return len(self._list)

	def __getitem__(self, ii):
		"""Get a list item"""
		return self._list[ii]

	def __delitem__(self, ii):
		"""Delete an item"""
		del self._list[ii]

	def __setitem__(self, ii, val):
		# optional: self._acl_check(val)
		return self._list[ii]

	def __str__(self):
		return str(self._list)

	def insert(self, ii, val):
		# optional: self._acl_check(val)
		self._list.insert(ii, val)

	def append(self, member):
		self.insert(len(self._list), member)

			
	def segregate(self,attribute,preds=[]):

		

		if attribute in self.attribute_class and self.attribute_class[attribute] != 'binary': 
			slist= sorted([(s.attributes[attribute],i) for i,s in enumerate(self._list)])
			slist= [g for g in slist if g[0] != 'NA']
			gIdx = len(slist)/2 
			g1 = slist[0:gIdx]


			while slist[gIdx][0] == g1[-1][0]: 
				g1.append(slist[gIdx])
				gIdx+=1 
				if gIdx == len(slist):
					seg =  {'HI': [g[1] for g in slist if g[0] !='NA'], 'LO': [g[1] for g in slist if g[0] == 'NA']} 
					return seg,{sa: len(seg[sa]) for sa in seg.keys()}

			g2 = slist[gIdx::]  
			seg =  {'HI': [g[1] for g in g2], 'LO': [g[1] for g in g1]} 

		else:
			seg = dd(list) 
			self.attribute_class[attribute] = 'binary'
			for i,s in enumerate(self._list):
				seg[s.attributes[attribute]].append(i) 
		return seg,{sa: len(seg[sa]) for sa in seg.keys()}


	def prune_cnts(self,min_val):

		self._list = [f for f in self._list if len(f.cnts) > min_val]
		self.len = len(self._list) 
		return self	
 


	def transpose(self): 
		self._list = [[self._list[i][j] for i in range(len(self._list))] for j in range(len(self._list[0]))]
		self.n_len,self.m_len,self.n_rng,self.m_rng = self.m_len,self.n_len,self.m_rng,self.n_rng
		return self

	def transform(self,transform):
		hide = False
		if   transform == 'log':	self._list = [[log(self._list[i][j]+1,2) for j in range(self.n_len)] for i in range(self.m_len)]
		elif transform == 'hide':  hide = True
		else:
			print transform,'wtf'
			sys.exit()
		return hide,self

	
	def center(self):
		self._list = [[self._list[i][j] - sum(self._list[i])/float(self.n_len) for j in range(self.n_len)] for i in range(self.m_len)]
		return self


	def filter(self,attribute_list,attribute_prune_list=[]):

		unknown_attributes = [p.split('=')[0] for p in attribute_prune_list if p.split('=')[0].split('!')[0] not in self.attributes]+[a.split('=')[0] for a in attribute_list if a.split('=')[0] not in self.attributes]

		if len(unknown_attributes)>0: 	predictor_error('Unknown Predictors: '+",".join(unknown_attributes))


	


		self.removed = [] 
		for a in attribute_list:

			
			if len(self._list[0].name.split('~'))>2 and self._list[0].name.split("~")[-2] == a: 
				for s in self._list:
					s.attributes[a] = s.name.split('~')[-1] 


			else:

				if self.attribute_class[a.split('=')[0]] == 'binary':	self.removed.extend([s.idx for s in self._list if s.attributes[a].split('~')[-1] == 'NA'])	
				else:							self.removed.extend([s.idx for s in self._list if s.attributes == 'NA'])	




		for p in attribute_prune_list:


			pA,pR,pX = p.split('=')[0], p.split('=')[0].split('!')[0],p.split('=')[-1].split(',') 

			if pR != pA: 
				self.removed.extend([s.idx for s in self._list if s.attributes[pR].split('~')[-1] in pX])
			else:
				self.removed.extend([s.idx for s in self._list if s.attributes[pA].split('~')[-1] not in pX])

		self.removed = sorted(list(set(self.removed)))
		
		if len(self._list) - len(self.removed) < 2: 
			predictor_error('Incongrous Covariates/Prune Requests: All samples removed') 

		self.swap = {}  
		self._list = [s for s in self._list if s.idx not in self.removed]
		for i in range(len(self._list)):	
			self.swap[self._list[i].idx] = i 
			self._list[i].idx = i 
			
		self.len = len(self._list) 
		return self
	



	def colorize_names(self,key):


		for i,s in enumerate(self._list):


			
			self._list[i].label 


			if s.name.split('=')[-1] in key: 
			
				self._list[i].label = key[s.name.split('=')[-1]]

				#s.notes['color'] = key[s.name.split('=')[-1]]
			elif s.name.split('~')[-1] in key:
				
				#s.notes['color'] = key[s.name.split('~')[-1]]
				self._list[i].label = key[s.name.split('~')[-1]]

			elif "~".join(s.name.split('~')[1::]) in key:
				self._list[i].label = key["~".join(s.name.split('~')[1::])]

			else: 
				s.notes['color'] = 'silver'














	def create_plot_labels(self,args,min_group_size=MIN_GRP_SIZE):

		self.notes['colors'], label_notes, self.colorize_key, binary_colors, continuous_colors, binary_markers = dd(list), {}, {}, {},{}, {} 
		colors, missed, marker, size, left_key, right_key, edge_key, size_labels, marker_labels = [], [], None, None, None,None,None, None, None 
		if 'color'  in args:    colors  = args.color 
		if 'marker' in args:    marker  = args.marker 
		if 'size'   in args:    size    = args.size 	
		if args.option == 'regression':
			
			p_bin   = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.predictors] for a in b] if self.attribute_class[p] == 'binary']
			p_cont  = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.predictors] for a in b] if self.attribute_class[p] != 'binary']
			c_bin   = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.covariates] for a in b] if self.attribute_class[p] == 'binary']
			c_cont  = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.covariates] for a in b] if self.attribute_class[p] != 'binary']
			variables = p_bin + p_cont + c_bin + c_cont 



#			preds = [p.split('=')[0] for p in args.predictors if self.attribute_class[p.split('=')[0]] == 'binary'] + [p.split('=')[0] for p in args.predictors if self.attribute_class[p.split('=')[0]] != 'binary'] 
#			covars = [p.split('=')[0] for p in args.covariates if self.attribute_class[p.split('=')[0]] == 'binary'] + [p.split('=')[0] for p in args.covariates if self.attribute_class[p.split('=')[0]] != 'binary'] 
			for r in variables: 
				if (r in colors) or (r == marker) or (r == size): continue 
				if len(colors) == 0: colors.append(r) 
				elif self.attribute_class[r] == 'binary': 
					if marker == None:  marker = r 
					elif len(colors) < 3: colors.append(r) 
					else: 		     missed.append(r) 
				else: 
					if size == None:     	size = r 
					elif len(colors) < 3:  	colors.append(r)       
					else:			missed.append(r) 
		for i,a in enumerate(colors):
			if self.attribute_class[a] == 'binary': 	binary_colors[a] = self.binary_labels(a,'color',rank=(i+1,len(colors))) 
			else:						continuous_colors[a] = self.continuous_labels(a,'color')


		if len(binary_colors) > 0: 
			left_key, right_key  = binary_colors.keys()[0],  binary_colors.keys()[-1] 
			if len(continuous_colors)> 0: 	edge_key  = continuous_colors.keys()[0]
		elif len(continuous_colors) > 0:
			left_key, right_key  = continuous_colors.keys()[0],  continuous_colors.keys()[-1] 

		if marker != None: 	binary_markers = self.binary_labels(marker,'marker') 
		if size != None: 
			if self.attribute_class[size] == 'binary': 
				print 'binary size'
				sys.exit() 
			else:
				s_data = np.array([s.attributes[size] for s in self._list]).reshape(-1,1)
				scaler = MinMaxScaler().fit_transform(s_data)
				size_labels = [(scaler[j][0]*7)+3 for j in range(len(scaler))]
		if left_key: 
			if left_key == right_key: self.notes['labels'] = {'color': left_key,'shape': marker, 'size': size} 
			else: 			  self.notes['labels'] = {'left_color': left_key, 'right_color': right_key, 'shape': marker, 'size': size}



		for i,s in enumerate(self._list):
			p = PlotLabel()


			if left_key:				
				try: p.set_left_color(binary_colors[left_key][1][i], binary_colors[left_key][0][i])
				except KeyError: p.set_left_color(continuous_colors[left_key][1][i], continuous_colors[left_key][0][i])
			if right_key and right_key != left_key:	
				try: p.set_right_color(binary_colors[right_key][1][i], binary_colors[right_key][0][i])
				except KeyError: p.set_right_color(continuous_colors[right_key][1][i], continuous_colors[right_key][0][i])
			if edge_key:				p.set_edge_color(continuous_colors[edge_key][1][i], continuous_colors[edge_key][0][i]) #  edge_labels[i], edge_colors[i])
			if size_labels: 			p.set_size_label(size_labels[i]) 
			if marker: 				p.set_marker_label(binary_markers[1][i], binary_markers[0][i]) #  marker_labels[i],markers[i]) 
			self._list[i].label = p.check() 
			self.colorize_key[p.id] =  self._list[i].label  
			





	def continuous_labels(self,attribute,label='color'):

		sIds = [s.attributes[attribute] for s in self._list if s.attributes[attribute] != 'NA']
		sAvg = np.mean(sIds)
		sImpute = [s.attributes[attribute] if s.attributes[attribute] != 'NA' else sAvg for s in self._list]
		if label == 'color':
			norm = plt.Normalize(min(sIds),max(sIds))
			return [plt.cm.hot(norm(x)) for x in sImpute], [attribute for x in sImpute] 








	def binary_labels(self,attribute,label='color',rank=None,min_group_size=MIN_GRP_SIZE,max_group_members=MAX_GRP_MEMBERS):

		
		sAll,sValid = [s.attributes[attribute] for s in self._list], [s.attributes[attribute] for s in self._list if s.attributes[attribute] != 'NA']
		min_group_size, max_group_members = 5, 10	

		if len(set(sValid)) > 0:	
			try: minSize = min(min_group_size,sorted(cc(sValid).values())[-2])
			except IndexError: minSize = 5
			sGrps, sIdx, sMissing = list(set(sValid)) , [], [] 
			if len(sGrps) > max_group_members:  sGrps = sorted([a for (a,b) in cc(sValid).items() if b > minSize],reverse=True)
			if len(sGrps) > max_group_members:  sGrps = sGrps[0:max_group_members]
			for s in self._list:
				if s.attributes[attribute] in sGrps:	
					sIdx.append(sGrps.index(s.attributes[attribute]))
				else:
					sIdx.append(len(sGrps))
					sMissing.append(s.attributes[attribute])
			if label == 'color': 
				label_list, missing_label = get_color_list(len(sGrps),len(sMissing)!=0,rank)
			elif label == 'marker': 
				label_list, missing_label = get_marker_list(len(sGrps),len(sMissing)!=0,rank)

			if len(sMissing) > 0:
				sGrps.append(attribute+'='+",".join([sM.split('~')[-1] for sM in list(set(sMissing))]))
				label_list.append(missing_label) 
		
			label_vals, label_labels = [label_list[i] for i in sIdx], [sGrps[i] for i in sIdx]
			return label_vals, label_labels 
				






























































'''	
	def colorize(self,attribute_list,min_group_size=MIN_GRP_SIZE):



		self.notes['colors'] = dd(list) 
		binary_colors = {} 
		continuous_colors = {} 
		for a in attribute_list:
			if self.attribute_class[a] == 'binary':	binary_colors[a] = self.binary_labels(a,'color',min_group_size) 
			else:					continuous_colors[a] = self.continuous_labels(a,'color')

		if len(binary_colors) > 0: 
			left_key, right_key  = binary_colors.keys()[0],  binary_colors.keys()[-1] 
			left_colors, left_labels = binary_colors[left_key]	
			right_colors, right_labels = binary_colors[right_key]
			primary_colors, edge_colors = [], []  
			for i,s in enumerate(self._list):
				primary_colors.append([left_colors[i],right_colors[i]])
				
			if len(continuous_colors)> 0: 
				edge_key  = continuous_colors.keys()[0]
				edge_colors, edge_labels = continuous_colors[edge_key]	
		
		for i,s in enumerate(self._list): 
			print [s.attributes[a] for a in attribute_list], primary_colors[i], edge_colors[i]

		return primary_colors,edge_colors
'''



#def method(**kwargs):
#  print kwargs

#keywords = {'keyword1': 'foo', 'keyword2': 'bar'}
#method(keyword1='foo', keyword2='bar')
#method(**keywords)









class PlotLabel2:
	def __init__(self):


		self.left_color  = None
		self.right_color = None
		self.edge_color  = None
		self.marker      = None 
		self.size        = 7
		self.label_id    = [] 
		self.labels = [] 

	def set_size_label(self,size_label): 

		self.size = size_label

	def set_marker_label(self,label,marker): 

		self.marker = marker 
		self.marker_label = label
		self.label_id.append([label,marker]) 

	def set_left_color(self,label,color):

		self.left_color = color
		self.left_color_label = label 
		self.label_id.append([label,color]) 

	def set_right_color(self,label,color):

		self.right_color = color
		self.right_color_label = label 
		self.label_id.append([label,color]) 
		
	def set_edge_color(self,label,color):
		self.edge_color = color
		self.edge_color_label = label
		self.label_id.append([label,color]) 



	def check(self,msize=10):


		self.id = '|'.join([L[0] for L in self.label_id]) 


		if self.left_color and self.right_color:	self.vals = {'color': self.left_color, 'fillstyle': 'left', 'marker': 'o', 'markerfacecoloralt': self.right_color, 'mew': 0.0}
		elif self.left_color:				self.vals = {'color': self.left_color, 'marker': 'o', 'mew': 0.0} #, 'markersize': self.size}
		else: 						self.vals = {'mew': 0.0} 
	

		if self.marker: self.vals['marker'] = self.marker 		
		self.proxy =  Line([0],[0],linestyle='none',markersize=msize,**self.vals)

		self.vals['markersize'] = self.size 


		return self























































	



























































	
