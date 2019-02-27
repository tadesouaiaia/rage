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
	for line in open(cnts_file):
		line = line.split()
		if line[0] == '---':
			names,sample_cnts = line[1::],[[] for i in range(len(line)-1)]

		else:
			

			genes.append(line[0]) 
			cnts = [float(x) for x in line[1::]]


	    		for j in range(len(cnts)): 
				sample_cnts[j].append(cnts[j]) 


	sample_obs = sorted([(len([x for x in sample_cnts[i] if x>0]),i) for i in range(len(sample_cnts))])

	if options.end == None: 
		options.end = len(sample_obs) 
	for i,(obsA,idxA) in enumerate(sample_obs):

		if i < options.start: continue 
		if i >= options.end:  break 
		cntA,nameA = sample_cnts[idxA] , names[idxA]
		
		for j in range(i+1,len(sample_obs)):
			obsB,idxB = sample_obs[j] 
			cntB,nameB = sample_cnts[idxB] , names[idxB]
			
			valid_idx = [x for x in range(len(cntA)) if cntA[x] > 0 or cntB[x] > 0] 

			valA,valB = [cntA[x] for x in valid_idx], [cntB[x] for x in valid_idx] 

			Sv = stats.spearmanr(valA,valB)[0] 
#			if Sv < 0.5 or np.isnan(Sv): continue 

			R = stats.pearsonr(cntA,cntB)[0] 
			Rv = stats.pearsonr(valA,valB)[0] 
			S = stats.spearmanr(cntA,cntB)[0] 


			print nameA,nameB,i,j,obsA,obsB,'|',round(R,4), round(Rv,4),'|',round(S,4), round(Sv,4)












if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)


    parser.add_option("-s", "--start", default = 0, type=int, help="Output Filename Prefix")
    parser.add_option("-e", "--end", default = None, type=int, help="Output Filename Prefix")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
    parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

    (options, args) = parser.parse_args()



    scan(args[0],options)













