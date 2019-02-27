#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 
import numpy as np

fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}



class GeneVals:
        #def __init__(self,cnt_key,color_key, options=None):
        def __init__(self,gene_key,sample_key,gene_lookup,sample_lookup,VERBOSE=False):

		self.gene_lookup = gene_lookup
		self.sample_lookup = sample_lookup
		self.remaining_genes = self.gene_lookup.keys() 
		self.slices = [200,100, 50,25,10,5,3,2,1] 
		self.mins =   [20, 10,  6, 4, 3,2,2,1,0] 
		self.genes   = gene_key
		self.samples = sample_key
		self.added_genes = [dd(int) for s in self.samples]
		self.sampleRange = range(len(self.added_genes))
		self.si, self.sj, self.iter = 0,0,0
		self.VERBOSE = VERBOSE
		

	def add_init(self,INIT_ADD):
		k = 0 
		t_added = []

		if self.VERBOSE: sys.stderr.write('adding initial high covered genes....') 
 
		self.top = sorted(self.genes['len'].items(),key = lambda X:X[1], reverse=True)
		while k < INIT_ADD:
			g_idx,gL = self.top[k] 
			self.genes['len'].pop(g_idx)
			for i in self.genes['samples'][g_idx]:
				self.added_genes[i]['obs'] += 1
				self.added_genes[i]['ranks'] += self.samples[i][g_idx]['rank']	
			t_added.append(g_idx) 
			k+=1

		if self.VERBOSE: sys.stderr.write('done\n') 


		self.summarize_progress(t_added) 


	def rebalance(self,RB_GOAL=10):

		if self.VERBOSE: sys.stderr.write('rebalancing geness....') 
		self.iter+=1
	 
		s_obs =  sorted([(self.added_genes[i]['obs'],i) for i in self.sampleRange])
		s_low, s_hi = s_obs[0:self.slices[0]], s_obs[-1::-1][0:self.slices[0]]


		top_remain = sorted(self.genes['len'].items(),key = lambda X:X[1], reverse=True) 	
		si,sj = 0,0 
		k,n,m=0,0,0


		
		t_obs = [self.samples[i][top_remain[0][0]]['exist'] for c,i in s_obs]
		TSUM =  sum([self.samples[i][top_remain[0][0]]['exist'] for c,i in s_obs])


		t_data,t_added  = [],[] 
		while sum(t_obs) > TSUM/4.0 and len(t_added) <= RB_GOAL: 

			t = top_remain[k][0]  
			t_obs =  [self.samples[i][t]['exist'] for c,i in s_obs]	
			sL,sH = t_obs[0:self.slices[si]]  , t_obs[-1::-1][0:self.slices[sj]]
			if sum(sL) > sum(sH): 
				for i in self.genes['samples'][t]:
					self.added_genes[i]['obs'] += 1
					self.added_genes[i]['ranks'] += self.samples[i][t]['rank']	
				t_added.append(t)
				self.genes['len'].pop(t)
			elif sum(sL) > self.mins[si]: t_data.append([sum(sL),t,t_obs,sL,sH])
			k+=1	
		BACKTRACE=False
		while len(t_added) < RB_GOAL:
			k=0
			t_data.sort(reverse=True) 

			if BACKTRACE or si +2 == len(self.slices):
				BACKTRACE=True
				if si > 0: si -= 1
				elif sj + 1 < len(self.slices): sj+=1
			else: 
				si +=1
				sj +=1 


			if len(t_added) > RB_GOAL or (si==0 and sj+1 == len(self.slices)): break 

			while k < len(t_data):
			
				scr,t,t_obs,sL,sH = t_data[k] 
				tL,tH = sL[0:self.slices[si]],sH[0:self.slices[sj]]
			
	
				if sum(tL) > sum(tH) and sum(tH) < len(tH):
					scr,t,t_obs,sL,sH = t_data.pop(k) 
					for i in self.genes['samples'][t]:
						self.added_genes[i]['obs'] += 1
						self.added_genes[i]['ranks'] += self.samples[i][t]['rank']	
					t_added.append(t)
					self.genes['len'].pop(t)
				else:	k += 1		
			
			 
		s_obs =  sorted([(self.added_genes[i]['obs'],i) for i in self.sampleRange])
		s_low, s_hi = s_obs[0:self.slices[0]], s_obs[-1::-1][0:self.slices[0]]



		if len(t_added) == 0 and len(t_data)>0: 
			scr,t,t_obs,sL,sH = t_data.pop(0) 
			for i in self.genes['samples'][t]:
				self.added_genes[i]['obs'] += 1
				self.added_genes[i]['ranks'] += self.samples[i][t]['rank']	
			t_added.append(t)
			self.genes['len'].pop(t)

		if self.VERBOSE: sys.stderr.write('done\n') 

		self.summarize_progress(t_added) 


		return

	def summarize_progress(self,t_added):
		rLen = len(self.remaining_genes)
		addLen    = len(t_added) 
		w = sys.stdout
		if self.iter == 0: 
			w.write('%-40s %10s %10s %10s %10s %10s %10s %10s %10s %10s %30s\n' % ('---','ITER','GENES','ADDED','REMAINING','MEANOBS','MIN_OBS','MIN_ID','MAX_OBS','MAX_ID','OBS10,25,50,75,90'))


		sample_stats = sorted([(self.added_genes[i]['obs'],self.added_genes[i]['ranks'],i) for i in self.sampleRange])


		minObs,minRank,minId = sample_stats[0] 
		maxObs,maxRank,maxId = sample_stats[-1] 
		obs_stat = [ss[0] for ss in sample_stats] 		

		obsMean = round(np.mean(obs_stat),2) 
		percs = " ".join([str(int(round(np.percentile(obs_stat,x)))) for x in [10,25,50,75,90]])

		for t in t_added:
			w.write('%-40s %10s %10s %10s %10s %10s %10s %10s %10s %10s %30s\n' %(self.gene_lookup[t],self.iter,rLen,addLen,rLen-addLen,obsMean,minObs,self.sample_lookup[minId],maxObs,self.sample_lookup[maxId],percs))

		self.remaining_genes = [gi for gi in self.remaining_genes if gi not in t_added]
		return







def scan(cnts_file,options):
	VERBOSE=True
	gene_lookup = {} 
	if VERBOSE: sys.stderr.write('READING FILE....') 
	for line in open(cnts_file):
		line = line.split()
		if line[0] == '---':
            		names,totals,gene_idx = line[1::],[[] for i in range(len(line)-1)],0
			sample_lookup = {n: names[n] for n in range(len(names))}


        	else:
			
            		gene,cnts = line[0],[float(x) for x in line[1::]]

			gene_lookup[gene_idx] = gene
 
			for i in range(len(cnts)): 
				if cnts[i] > 0: totals[i].append((cnts[i],gene_idx))
			gene_idx +=1 


	if VERBOSE: sys.stderr.write('COMPLETE\n') 
	if VERBOSE: sys.stderr.write('SETTING UP INDEXES....') 

	sample_key = [dd(lambda: dd(int)) for i in range(len(totals))]
	gene_key =   dd(lambda: dd(list)) 
	for i,t in enumerate(totals): 

		ts = [x for x in sorted(t,reverse=True)] 
		sample_key[i]['obs'] = len(ts) 
		for k,(gene_val,gene_idx) in enumerate(ts):
			sample_key[i][gene_idx]['exist'] = 1
			sample_key[i][gene_idx]['rank'] = k+1
			sample_key[i][gene_idx]['val'] = gene_val
			gene_key['ranks'][gene_idx].append(k+1) 
			gene_key['samples'][gene_idx].append(i) 
	
	for gene_idx,ranks in gene_key['ranks'].items():

		gene_key['len'][gene_idx]   = len(ranks) 
		gene_key['avg'][gene_idx]   = np.mean(ranks) 
		gene_key['ranks'][gene_idx] = sorted(ranks) 
		gene_key['obs'][gene_idx] = len(ranks) 

	
	if VERBOSE: sys.stderr.write('COMPLETE\n') 
	INIT_ADD = 500
	
	genevals = GeneVals(gene_key,sample_key,gene_lookup,sample_lookup,VERBOSE) 		
	genevals.add_init(INIT_ADD) 

	for rb in range(400):
		genevals.rebalance(10)




			 
















if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"

	parser = OptionParser(usage=usage)


	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	(options, args) = parser.parse_args()



	scan(args[0],options)













