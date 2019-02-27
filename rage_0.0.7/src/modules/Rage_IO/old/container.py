import sys
from collections import defaultdict as dd
from collections import Counter as cc
from collections import MutableSequence
import pickle as pk 
import numpy as np
from math import log
from sklearn import linear_model
from sklearn import feature_selection


def RageError(msg,error):
	raise error(msg) 







class ddd(dd):
    	def __missing__(self, key):
        #value = self.default_factory(key)
        	self[key] = 'NA'
        	return 'NA'


class array_check:
	def __init__(self,cstr):
		self.error_type = error_type
		print 'yo',cstr
	
		
class Ct(MutableSequence):


	"""A container for manipulating lists of hosts"""
	def __init__(self, data=None):
		"""Initialize the class"""
		super(Ct, self).__init__()
		if (data is not None):

			d,self.annotation,self.size = [x for x in data],None,0


			while type(d) == list:
				if len(d) == 0:
					break
				else:
					d = d[0] 
					self.size+=1


			if self.size == 3: 			
				try:			self._list = [[[float(v) for v in w] for w in V] for V in data]
				except ValueError:	self._list = list(data) 
				self.k_len,self.n_len,self.m_len = len(data),len(data[0][0]),len(data[0])
				self.k_rng,self.n_rng,self.m_rng = range(self.k_len),range(self.n_len),range(self.m_len)


			elif self.size == 2:
 
				try:			self._list = [[float(v) for v in V] for V in data]
				except ValueError:	self._list = list(data) 
				self.k_len,self.n_len,self.m_len = 1,len(data[0]),len(data) 
				self.k_rng,self.n_rng,self.m_rng = range(self.k_len),range(self.n_len),range(self.m_len)


			elif self.size == 1:
				try:			self._list = [float(v) for v in data]
				except ValueError:	self._list = list(data) 
				self.k_len,self.n_len,self.m_len = 1,None,len(data) 
				self.k_rng,self.n_rng,self.m_rng = [],[],range(self.m_len)

			else:
				print self.size,'wtf,'


				
		else:
	    		self._list = list()
			self.n_len,self.m_len = None,0 
			self.n_rng,self.m_rng = [],[]

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

	def append(self, val):
		self.insert(len(self._list), val)

	def match(self, d_list,verbose=False):
		if len(self._list) != len(d_list): return False	
		for i in range(len(self._list)):
			if self._list[i] != d_list[i]: return False
		if self.annotation and self.annotation != d_list.annotation: return False
			
		
			



		return True
			
	def annotate(self,name,vals,axis,members=['self']):

                self.annotation,self.ID,self.name,self.axis,self.mem = ddd(),name,axis+"."+name,axis,members
                if not self.mem or self.mem == 'input': self.mem = ['self']

                if name == axis: self.type = 'text'


		else:
			try: 	
				floats,ints = [float(v) if v!= 'NA' else 'NA' for v in vals],[int(float(v)) if v!='NA' else 'NA' for v in vals]
				diffs = sum([floats[i]-ints[i] for i in range(len(floats)) if floats[i]!='NA'])
				if diffs == 0:  vals,self.type = ints,'discrete'
				else:		vals,self.type = floats,'continuous'
			except TypeError:  self.type = 'array'
			except ValueError: self.type = 'binary'


		if self.type == 'array':

			for i in range(len(self._list)):	
				try: vals[i] = [float(v) for v in vals[i]]
				except ValueError: vals[i] = vals[i]
				self.annotation[self._list[i]] = vals[i]	
                	self.opts = list(set(sorted([a for b in self.annotation.values() for a in b])))
		else:
			for i in range(len(self._list)):	self.annotation[self._list[i]] = vals[i]	
			self.opts = list(set(sorted([v for v in self.annotation.values() if v!='NA'])))

		return self
	
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



