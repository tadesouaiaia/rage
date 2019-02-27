import sys
from collections import defaultdict as dd
from collections import Counter as cc
from collections import MutableSequence
import numpy as np
from math import log
import mmap 

def RageError(msg,error):
	raise error(msg) 







class cnt_key(dd):
    	def __missing__(self, key):
        #value = self.default_factory(key)
        	#self[key] = 'NA'
        	return 0.0

class attribute_key(dd):
    	def __missing__(self, key):
        #value = self.default_factory(key)
        	self[key] = 'NA'
        	return 'NA'


class array_check:
	def __init__(self,cstr):
		self.error_type = error_type
		print 'yo',cstr
	



class read: 
	def __init__(self,rage):

		self.rage = rage 
		rage.progress.start_major('Reading Input') 	
		
		if len(rage.args.counts) > 0:
			self.read_primary_counts(rage.args.counts[0]) 
 			for c,cnt_file in enumerate(rage.args.counts[1::]): 
				self.read_secondary_counts(cnt_file,'cnts-'+str(2+c)) 

		if 'sampleKey' in vars(rage.args) and rage.args.sampleKey != None: self.add_sample_attributes(rage.args.sampleKey)





	
	def read_primary_counts(self,cnt_file): 
		self.rage.progress.start_minor('Reading Primary Counts',5000,True) 	
		f_line = cnt_file.readline().split()[1::]
		self.features, kv = Members(), 0 
		if self.rage.args.test:
			self.test_idx = [i for i in range(len(f_line)) if i % 20 == 0]
			self.samples = Members([f_line[i] for i in self.test_idx]) 
			for j,line in enumerate(cnt_file): 
				line = line.split() 
				fr,cnts = line[0],[float(x) for x in line[1::]]
				test_cnts = [cnts[i] for i in self.test_idx] 
				if len([x for x in test_cnts if x >0]) > 10: 
					self.features.append(Member(fr,kv).add_line_cnts(test_cnts))
					for i,c in enumerate(test_cnts): 
						if c > 0: self.samples[i].cnts[kv] = c
					kv+=1
				if kv > self.samples.len: break 
		else: 
			self.samples = Members(f_line) 
			for j,line in enumerate(cnt_file): 
				self.rage.progress.mark() 
				line = line.split() 
				fr,cnts = line[0],[float(x) for x in line[1::]]
				if sum(cnts) > 0: 
					self.features.append(Member(fr,kv).add_line_cnts(cnts))
					for i,c in enumerate(cnts): 
						if c > 0: self.samples[i].cnts[kv] = c
					kv+=1
		self.samples.collate() 
		self.features.collate()


	def read_secondary_counts(self,cnt_file,c_idx): 
		self.rage.progress.start_minor('Reading Secondary Counts',5000,True) 	
		f_line = cnt_file.readline().split()[1::]
		kv = 0 
		self.test_idx = [i for i in range(len(f_line)) if f_line[i] in self.samples]
		for j,line in enumerate(cnt_file): 
				line = line.split() 
				fr,cnts = line[0],[float(x) for x in line[1::]]
				test_cnts = [cnts[i] for i in self.test_idx] 
				if fr in self.features.names:  
					for i,c in enumerate(test_cnts): 
						f_idx = self.features.lookup[fr] 
						if c > 0: 
							self.features[f_idx].notes[c_idx][i] = c 
							self.samples[i].notes[c_idx][kv] = c 
				
					kv+=1 
					if kv == self.features.len: break 



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
					if s_name in self.samples.lookup: self.samples[self.samples.lookup[s_name]].attributes[opt] = s_val 	
				self.samples.add_attribute(opt,'binary') 






	def data_matrix(self,transform,scaling=None):
		if transform == 'log': 

			dMatrix =  np.matrix([[0 if f.idx not in s.cnts else log(1.0+s.cnts[f.idx]) for f in self.features] for s in self.samples])

                     	return dMatrix - dMatrix.mean(axis=0)








		 
		
		
			
class Member:
	def __init__(self,name,idx):
		self.name = name 
		self.idx  = idx 
		self.cnts = cnt_key(None) 
		self.attributes = attribute_key(None)
#		self.cnts = ddd(None) 
#		self.attributes = ddd(None)


		self.notes = dd(bool) 
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






class Members(MutableSequence):

	"""A container for manipulating lists of hosts"""
	def __init__(self, data=None):
		"""Initialize the class"""
		super(Members, self).__init__()
		if (data is not None):
			self._list = [Member(data[i],i) for i in range(len(data))]
			self.len = len(self._list) 
			self.attributes,self.attribute_type = [], [] 
		else:
	    		self._list = list()
			self.len = 0 
			self.attributes,self.attribute_type = [], [] 


	def collate(self): 

		self.names = [x.update() for x in self._list]
		self.lookup = {member.name : member.idx for member in self._list}
		self.len = len(self._list) 

	def add_attribute(self,opt,optType):
		self.attributes.append(opt) 
		self.attribute_type.append(optType) 

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

			
	def segregate(self,m_idx,n_idx):
		return  [[self._list[i][j] for j in n_idx] for i in m_idx]
		#self._list = [[self._list[i][j] for j in m_idx] for i in n_idx]
		#return self

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

	

