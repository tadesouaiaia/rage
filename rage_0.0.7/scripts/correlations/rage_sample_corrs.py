#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 
import numpy as np 
import warnings
warnings.filterwarnings("ignore")





def scan(cnts_file,options):
	genes = [] 
	len_key = []  
	rate_key = []  
	for line in open(cnts_file):
		line = line.split()
		if line[0] == '---':
			names,sample_cnts = line[1::],[[] for i in range(len(line)-1)]

		else:
			

			genes.append(line[0]) 
			cnts = [float(x) for x in line[1::]]


	    		for j in range(len(cnts)): 
				sample_cnts[j].append(cnts[j]) 

	tGene = float(len(genes)) 
	for i in range(len(sample_cnts)): 
		len_key.append(len([x for x in sample_cnts[i] if x > 0]))
		rate_key.append(len_key[i]/tGene) 


	print '--- --- obs1 obs2 mix_obs exp_mix |','dist','|','Rall Sall Rpa Spa','|','Rkeep Skeep Rpk Spk'
	for i in range(len(sample_cnts)): 
		sI = names[i] 
		cI = sample_cnts[i] 
		iL,iR = len_key[i] , rate_key[i] 
		cIL = [log(1+x,2) for x in cI]



		for j in range(i+1,len(sample_cnts)): 
			sJ = names[j] 
			cJ = sample_cnts[j] 
			jL, jR = len_key[j] , rate_key[j] 
			cJL =  [log(1+x,2) for x in cJ] 

    			iV = [cI[z] for z in range(len(cI)) if cI[z] != 0 and cJ[z] != 0] 
    			jV = [cJ[z] for z in range(len(cI)) if cI[z] != 0 and cJ[z] != 0] 


			mixL = len(iV) 

			cIL,cJL = [log(1+x,2) for x in cI], [log(1+x,2) for x in cJ] 
			bIL,bJL = [log(x,2) for x in iV], [log(x,2) for x in jV] 
			R_ALL = stats.pearsonr(cIL,cJL) 
			S_ALL = stats.spearmanr(cIL,cJL) 

			R_KEEP = stats.pearsonr(bIL,bJL) 
			S_KEEP = stats.spearmanr(bIL,bJL) 
			match_exp = iR*jR*tGene


			
			dist =  round(np.linalg.norm( np.array(cIL) - np.array(cJL)),4) 
			print sI,sJ,iL,jL,mixL,int(match_exp),'|',dist, 

			print '|', round(R_ALL[0],3), round(S_ALL[0],3), round(R_ALL[1],5), round(S_ALL[1],5),'|', round(R_KEEP[0],3), round(S_KEEP[0],3), round(R_KEEP[1],5), round(S_KEEP[1],5)
			











if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)


    parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
    parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

    (options, args) = parser.parse_args()



    scan(args[0],options)













