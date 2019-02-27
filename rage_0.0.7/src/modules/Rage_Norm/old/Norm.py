#!/usr/bin/env python

import os
import sys
from collections import defaultdict as dd
import random
import math
import os
import sys
from math import log
import random
import pickle
#import Types as Type
from types.dataset import Dataset

def mean(mylist):
	return sum(mylist)/float(len(mylist))

def swap(mylist):
	return [[mylist[i][j] for i in range(len(mylist))] for j in range(len(mylist[0]))]

class Norm:
	def __init__(self,args,choices=None):
                self.args, self.cmd,self.opts, self.results = args, args.option+'.'+args.command,"|".join(choices), []		


		if self.args.log == True: 	self.transform = 'log'
		else:			 	self.transform = None

		self.output = sys.stderr



	def talk(self,command):
		if not self.args.silent:
			self.output.write("Rage Norm: ")
			self.output.write(command)


	def run(self,data):	
		
		if self.args.command == 'downsample' or self.args.command == 'multi':
			self.talk("Begin Downsampling\n")
			self.downsample()	

		if self.args.command == 'rank' or self.args.command == 'quantile' or self.args.command == 'multi'  or self.args.command == 'top':
			self.rank_normalize()
			if self.args.command == 'quantile' or self.args.command == 'multi':
				self.produce_quantiles()
			
			if self.args.command == 'top':
				self.produce_tops() 


		if self.args.command == 'rpm' or self.args.command == 'multi':
			self.rpm_normalize()

		

		if self.args.command == 'condensation':
			self.condense()

		return self.results


	def downsample(self):

		D,nKey,self.args.size = self.args.dataset,{},self.args.size[0]

		nVals = [[self.args.dataset.v[j][i] for j in range(len(self.args.dataset.m))] for i in range(len(self.args.dataset.n))]

		nTotals = sorted([(int(sum([self.args.dataset.v[j][i] for j in range(len(self.args.dataset.m))])),self.args.dataset.n[i]) for i in range(len(self.args.dataset.n))])
		if self.args.size == 0: 
			self.args.size = int(nTotals[int(len(nTotals)/10.0)][0])
			
			nLen = float(len(nTotals))
			perc5,perc10,perc20,perc25,perc75,perc80,perc90,perc95 = int(nLen*0.05),int(nLen*0.10),int(nLen*0.20),int(nLen*0.25),int(nLen*0.75),int(nLen*0.8),int(nLen*0.9),int(nLen*0.95)


			nFloats = [float(N[0]) for N in nTotals][perc25:perc75]

			print nFloats[0:20]

			nMean = sum(nFloats)/len(nFloats)
			nStd  = (sum([(nf-nMean)*(nf-nMean) for nf in nFloats])/(len(nFloats)-1.0))  

			print nMean,nStd
			
			


			self.talk("No minimum size supplied - will use 10th percentile: "+str(self.args.size)+"\n")
		else:
			self.talk("Using supplied minimum downsample value: "+str(self.args.size)+"\n")
		

		nSmall = [x for x in nTotals if x[0]*5 < self.args.size]
		if len(nSmall)>0:
			nList = ",".join([x[1] for x in nSmall])
			self.talk("Warning: Outlier samples with low counts will be heavilty upsampled: "+nList+"\n")



		for i,V in enumerate(nVals):
			idx, iRand,nKey[D.n[i]], k = [[(0,0),'INIT']], [], dd(int), 0
			for j,v in enumerate(V):
				if v >= 1: idx.append([(idx[-1][0][1],idx[-1][0][1]+int(v)),D.m[j]])
			iRange,iLen = range(idx[-1][0][1]+1),len(range(idx[-1][0][1]+1))



			while len(iRand) < self.args.size: iRand += random.sample(iRange,min(iLen,self.args.size-len(iRand)))
			for s in  sorted(iRand):
				while s > idx[k][0][1]: k+=1
				nKey[D.n[i]][idx[k][1]]+=1
		vals = [[nKey[s][D.m[j]] for s in D.n] for j in range(len(D.m))]

		if self.args.log: name = 'ds.'+str(self.args.size)+'.log.cnts'
		else: 	          name = 'ds.'+str(self.args.size)+'.cnts'

              	self.results.append(Dataset(name,False,self.transform).create(D.M,D.N,'dowsampled-cnts',D.m,D.n,vals))
		self.talk("Downsampling Complete\n")

		
	def condense(self):
		D,key,kName = self.args.dataset,self.args.n_keys['groups'],self.args.n_keys['groups'].names[0]
		kData,kType,kOpts = [[k[2],1,'NO'] for k in key.data[kName]],key.types[kName],list(set([k[2] for k in key.data[kName] if k[0]]))
		
		
		if not key.default:	groups = [[i for i in range(len(kData)) if kData[i][0]==kOpts[j]] for j in range(len(kOpts))]

		else:
			print key.default
			sys.exit()

		for size in self.args.size:
				
			size,grps = max(size,min(min([int(len(g)/2.0) for g in groups]),2)),[g for g in groups]

			for i,group in enumerate(grps):
				random.shuffle(group)
				grps[i] = [group[j:j+size] for j in range(0,len(group),size)]


			n =          [x for x in D.n]+[a for b in [[kOpts[j]+'@'+str(i+1) for i in range(len(grps[j]))] for j in range(len(kOpts))] for a in b]
			gk = kData + [a for b in [[[kOpts[j],len(g),'YES'] for g in grps[j]] for j in range(len(kOpts))] for a in b]
			gv=[v+[a for b in [[sum([v[x] for x in g]) for g in grps[i]] for i in range(len(grps))] for a in b] for j,v in enumerate(D.v)]
			gP,gN = 'condense'+str(size)+'.',D.N+'_c'+str(size)
			#gP,gN = ".".join([D.name.split(".")[0],'condense',kName,str(size),".".join(D.name.split(".")[1:-1])]),D.N+'.'+str(size)+'.condense'
			if self.args.log: gv = [[log(x+1) for x  in X] for X in gv]


			self.results.append(Dataset(gP+D.V).create(D.M,gN,D.V,D.m,n,gv))
			self.results[-1].parents['N'].append(D.N)
              		self.results.append(Dataset(gP+'key',False,'hide').create(gN,'info','anno',n,[kName,'size','condensed'],gk,gN))



	def rank_normalize(self):

		

		self.sort_tuples = [sorted([(self.args.dataset.v[j][i],self.args.dataset.m[j]) for j in range(len(self.args.dataset.m))]) for i in range(len(self.args.dataset.n))]
		self.rank_lists = [[] for d in self.args.dataset.m]

		self.rank_map,self.range_map,self.quantile_map = dd(lambda: dd(float)) , dd(lambda: dd(float)), dd(lambda: dd(float))

		for i,S in enumerate(self.sort_tuples):
			m,n,scr,rank_tuple,rank_avg,group=0,0,S[0][0],[],[],[]
			while True:
				while n < len(S) and S[n][0] == scr: 
					group.append(S[n][1])
					n+=1
				if scr == 0: rank_tuple.extend([(g,0.0,0.0,0.0) for g in group])
				else:	     rank_tuple.extend([(g,m,n-1,scr) for g in group])


				if n == len(S): break 
				m,scr,group=n,S[n][0],[]

			sample = self.args.dataset.n[i]
			for k,ra in enumerate(sorted(rank_tuple)):
				self.rank_map[ra[0]][sample] = (ra[2]-ra[1])/2.0
				self.range_map[ra[0]][sample] = (ra[1],ra[2])
				self.rank_lists[k].append(ra[3])
				

		
		rank_sums = sorted([sum(RL)/float(len(RL)) for RL in self.rank_lists])
		
		for gene in self.range_map.keys():
			for sample,(a,b) in self.range_map[gene].items():
				if b == 0:	self.quantile_map[gene][sample] = 0 
				else:		self.quantile_map[gene][sample] = sum([rank_sums[x] for x in range(a,b+1)])/(float(b+1)-a)

		
		self.quantile_vals = [[self.quantile_map[m][n] for n in self.args.dataset.n] for m in self.args.dataset.m]
		self.rank_vals     = [[self.rank_map[m][n] for n in self.args.dataset.n] for m in self.args.dataset.m]
		
	def produce_quantiles(self):

		if self.args.log:	name = 'quantile.log.cnts'
		else:			name = 'quantile.raw.cnts'

              	self.results.append(Dataset(name,False,self.transform).create(self.args.dataset.M,self.args.dataset.N,'quantile-cnts',self.args.dataset.m,self.args.dataset.n,self.quantile_vals))
		self.talk("Quantile Normalization Complete\n")

		


	





	


	def rpm_normalize(self):


		D = self.args.dataset
		nVals = [[self.args.dataset.v[j][i] for j in range(len(self.args.dataset.m))] for i in range(len(self.args.dataset.n))]
		nTotals = [float(sum([self.args.dataset.v[j][i] for j in range(len(self.args.dataset.m))])) for i in range(len(self.args.dataset.n))]

		constant = sum(nTotals)/float(len(nTotals))


		self.rpm_vals =  [[(self.args.dataset.v[j][i]*constant)/(nTotals[i]+0.01) for i in range(len(self.args.dataset.n))] for j in range(len(self.args.dataset.m))]


		if self.args.log:	name = 'rpm.log.cnts'
		else:			name = 'rpm.raw.cnts'

              	self.results.append(Dataset(name,False,self.transform).create(self.args.dataset.M,self.args.dataset.N,'rpm-cnts',self.args.dataset.m,self.args.dataset.n,self.rpm_vals))
		self.talk("RPM Normalization Complete\n")















	def size_factor_normalize(self):
		log_avgs = [sum([math.log(v+1) for v in  self.data.feature_vals[j]])/len(self.data.samples) for j in range(len(self.data.features))]
		size_factors = [math.exp(sorted([log(s[j]+1) - log_avgs[j] for j in range(len(log_avgs))])[len(log_avgs)/2]) for s in self.sample_cnts]
		
		fname = 'sizefact'
		self.write(fname)
		for j in range(len(self.data.features)):
			feature,cnts = self.data.features[j],[self.data.feature_vals[j][i]*size_factors[i] for i in range(len(self.data.samples))]
			self.write(fname,feature,cnts)



