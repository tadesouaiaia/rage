#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as multitest 


import statsmodels.sandbox.stats.multicomp as mpt 
#import statsmodels.sandbox.stats.multicomp.fdrcorrection0 as fdr

#statsmodels.sandbox.stats.multicomp.multipletests(p,alpha=0.05,method='fdr_bh')
#statsmodels.sandbox.stats.multicomp.fdrcorrection0(p,alpha=0.05)


class out_line:
        def __init__(self,line,missing_grp='FULLGROUP~ES'):
		self.grps = [] 
		self.missing_grp = missing_grp
		for i in range(12,len(line),5): 
			self.grps.append(line[i]) 




	def add_line(self,line): 
		self.key = {} 
		self.gene = line[0] 
		self.RS, self.CV, self.obs, self.len = [float(x) for x in line[1:5]]
		self.OR = self.obs / self.len 
		self.parent = line[5] 
		self.max_grp = line[6] 
		self.maxMean, self.maxObs, self.maxChi, self.maxP = [float(x) for x in line[7:11]]
		self.mp = {} 
		self.np = {} 
		self.bool = {} 
		psums = [] 
		for i in range(12,len(line),5): 
			grp,gbool,nullP,modelP = line[i],line[i+1],float(line[i+2]),float(line[i+3])
			self.key[grp] = [gbool,nullP,modelP]
			self.np[grp], self.mp[grp] = nullP, modelP
			self.bool[grp] = gbool 
			psums.append((nullP+modelP,grp))

		psums.sort() 
		nps =  [self.np[g] for g in self.grps] 
		mps =  [self.mp[g] for g in self.grps] 
		mybools = [self.bool[g] for g in self.grps]
		gbools = list(set([self.bool[g] for g in self.grps]))
		NP,MP = np.mean(nps), np.mean(mps) 
		self.over, self.under = None, None 
		if self.max_grp in self.grps: self.over = [self.max_grp, self.OR, self.maxMean,self.maxObs,self.maxChi,self.np[self.max_grp], self.mp[self.max_grp]] 
		if len(gbools) == 1 and gbools[0] == 'False': 
			if self.max_grp not in self.grps:  
				self.over = [self.max_grp, self.OR, self.maxMean,self.maxObs,self.maxChi,NP, MP] 
				self.under = [psums[0][-1], 'NA', 'NA', 'NA', 'NA', self.np[psums[0][-1]], self.mp[psums[0][-1]]]

 
		elif len(gbools) == 1 and gbools[0] == 'True': 
			if self.max_grp in self.grps: 
				self.under = [self.missing_grp, 'NA', 'NA','NA','NA',NP, MP] 
				
		else: 
			
			gN,gI = 0,1  
			if mybools[0] != False: 
				gN,gI = 1, 0  




			self.under = [self.grps[gI], 'NA', 'NA', 'NA', 'NA', self.np[self.grps[gI]], self.mp[self.grps[gI]]]
			self.over = [self.grps[gN], 'NA', 'NA', 'NA', 'NA', self.np[self.grps[gN]], self.mp[self.grps[gN]]]
			
		if self.over == None: 
			if len(gbools) == 1 and gbools[0] == 'True': 
				self.over = [psums[-1][-1], 'NA', 'NA', 'NA', 'NA', self.np[psums[-1][-1]], self.mp[psums[-1][-1]]]
			 
			
		if self.under == None: 
			self.under = [psums[0][-1], 'NA', 'NA', 'NA', 'NA', self.np[psums[0][-1]], self.mp[psums[0][-1]]]

		return 





#def parse_out_file(line): 

#			--- RS CV obs len parent maxG maxMeans maxObs maxChi maxP | params

#			['ENSG00000000971;chr1;CFH', '0.007', '3.537', '276', '2455', 'FULLGROUP', 'FULLGROUP~ES', '0.34', '0.19', '2.6e-05', '1.3e-03', '|', 'FULLGROUP~AB', 'False', '2.51e-03', '2.51e-03', '|', 'FULLGROUP~EB', 'False', '3.00e-05', '3.00e-05']



def run_script(data_file,options):
	k=0
	genes, p_vals, gene_data = [], [] , []
	overs, unders,obs = [], [] , {} 
	pm_under, pn_under, pm_over, pn_over = [], [], [], [] 
	 
	for line in open(data_file): 
		line = line.split() 
		if line[0] == '---': header = line 
		else:
			if k == 0: 
				dex = out_line(line) 
 		
			if '|' in line[0] or ';RP' in line[0]: continue 
			dex.add_line(line) 
			genes.append(dex.gene)
			obs[dex.gene] = dex.OR  
			overs.append(dex.over) 
			unders.append(dex.under) 
			pn_under.append(dex.under[-2])
			pm_under.append(dex.under[-1])
			pn_over.append(dex.over[-2])
			pm_over.append(dex.over[-1])
			k+=1
#			genes.append(line[0]) 
#			gene_data.append(line[6]) 
#			p = float(line[-1]) 
#			p = float(line[10]) 
#			p_vals.append(p) 
			#k+=1
		#if k > 10: break 
		
	over_fdr_bool, over_fdr_res =  mpt.fdrcorrection0(pm_over,alpha=0.05) 
	under_fdr_bool, under_fdr_res =  mpt.fdrcorrection0(pm_under,alpha=0.05) 

#	over_fdr_bool, over_fdr_res =  mpt.fdrcorrection0(pn_over,alpha=0.05) 
#	under_fdr_bool, under_fdr_res =  mpt.fdrcorrection0(pn_under,alpha=0.05) 
	
	#print mpt.multipletests(p_vals,alpha=0.05,method='fdr_bh') 
	#print 'hu' 	
#	for g,b,f in zip(genes,over_fdr_bool,over_fdr_res): 
#		if b == True:
#			print g,b,f

#	print len(genes), len(under_fdr_bool), len(under_fdr_res) 
#	sys.exit() 

	for g,b,f,d in zip(genes,over_fdr_bool,over_fdr_res,overs): 

		#if b == True:
		print g,obs[g],'over',b,f,d[0]


	for g,b,f,d in zip(genes,under_fdr_bool,under_fdr_res,unders): 
		#if b == True:
		print g,obs[g],'under',b,f,d[0] 










if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args[0],options)	














