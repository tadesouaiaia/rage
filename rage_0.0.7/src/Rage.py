#!/usr/bin/env pythonw

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
from math import log
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import math
from random import shuffle
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import MinMaxScaler
from math import exp
from math import factorial 
from modules.Rage_IO import rage_inputs, rage_outputs, rage_progress



	



class Rage:
        def __init__(self,args,command_line):

		self.args, self.p_add = args, "."+args.command


		if len(args.counts) + len(args.condensedcnts) > 0: 

			self.args.prefix = self.create_count_prefix(args.prefix,args.counts,args.condensedcnts,args.test) 

			for a in vars(self.args):
				if (type(vars(args)[a]) == list) and (a not in ['counts','condensedcnts']) and (len(vars(args)[a])>0):	self.p_add+='_'+a+'__'+"_".join(vars(args)[a])

			self.progress = rage_progress.dot_log(self.args,self.args.prefix+self.p_add,command_line)
			if len(args.counts) ==1:		
				self.data = rage_inputs.read(self) 
				if len(args.condensedcnts) == 1: self.condensed_data = rage_inputs.read(self,data_type='CONDENSED') 
			elif len(args.counts) > 1:		self.datasets = [rage_inputs.read(self,i) for i,c in enumerate(args.counts)]
			elif len(args.condensedcnts) > 1:	self.datasets = [rage_inputs.read(self,data_idx=i,data_type='CONDENSED') for i,c in enumerate(args.condensedcnts)]
			elif len(args.condensedcnts) == 1:	self.data = rage_inputs.read(self,data_type='CONDENSED') 

			if args.option == 'norm': 
				from modules.RageNorm import Norm
				norm = Norm(self).run() 
			elif args.option == 'summarize': 
				from modules.RageSummary import Summary
				summary = Summary(self).run()  	
			elif args.option == 'regression': 		
				from modules.RageRegression import Regression
				regression = Regression(self).run() 
			elif args.option == 'transform': 		
				from modules.RageTransform import Transform
				transform = Transform(self).run() 
			elif args.option == 'simulate':
				from modules.RageSimulate import Simulate
				simulation = Simulate(self).run() 
			elif args.option == 'classify': 		
				from modules.RageClassify import Classify
				classification = Classify(self).run() 
			else:
				self.RageInputError('Counts File Required (--counts)') 

			self.progress.end() 
						
				


	def create_count_prefix(self,prefix,cnts,condensedcnts,TEST):
		if prefix: return prefix 
		else: 
			if len(cnts) == 1:		filepath = cnts[0].name 
			elif len(condensedcnts) == 1: 	filepath = condensedcnts[0].name 
			elif len(cnts) > 1: 		filepath = "_".join([x.name.split('/')[-1].split('.cnts')[0] for x in cnts])
			else: 				filepath = 'MULTI_RUN'
			if TEST: 
				return filepath.split('/')[-1].split('.cnts')[0]+'.testout'
			else:
				return filepath.split('/')[-1].split('.cnts')[0]+'.rageout'









	def RageInputError(self,msg):
		sys.stderr.write('Rage Input Error:\n') 
		sys.stderr.write('  Option:   '+self.args.option+'\n')
		sys.stderr.write('  Command:  '+self.args.command+'\n')
		sys.stderr.write('  Error:    '+msg+'\n')
		sys.exit()  
