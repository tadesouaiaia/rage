#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 






def scan(cnts_file,options):
	genes = []
	gene_cnts= []  
	for line in open(cnts_file):
		line = line.split()
		if line[0] == '---':
			names,sample_cnts = line[1::],[[] for i in range(len(line)-1)]

		else:
			

			genes.append(line[0]) 
			gene_cnts.append([float(x) for x in line[1::]])
			

#	    		for j in range(len(cnts)): 
#				sample_cnts[j].append(cnts[j]) 

	print '--- --- | shared missed | pv1 pv2 pv3 pv4 |','Rx,Sx,Rv,Sv'
	for i in range(len(gene_cnts)):

		gI,cI = genes[i],gene_cnts[i] 
		iL = [log(x+1,2) for x in cI] 
		for j in range(i+1,len(gene_cnts)):

			gJ,cJ = genes[j],gene_cnts[j]
			jL = [log(x+1,2) for x in cJ] 


			gv = [k for k in range(len(cI)) if (cI[k]+cJ[k]) > 0 ] 

			gshare = [k for k in range(len(cI)) if cI[k] > 0 and cJ[k] > 0] 

			shared = len(gshare) 
			missed = len(gv) - shared
			viL,vjL = [iL[k] for k in gv],[jL[k] for k in gv] 

			Rx = stats.pearsonr(iL,jL)
			Sx = stats.spearmanr(iL,jL) 
			Rv = stats.pearsonr(viL,vjL)
			Sv = stats.spearmanr(viL,vjL) 

			
			print gI, gJ, '|',shared,missed,'|',Rx[1],Sx[1], Rv[1],Sv[1],'|',
			print round(Rx[0],3), round(Sx[0],3),
			print round(Rv[0],3), round(Sv[0],3)













if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)


    parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
    parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

    (options, args) = parser.parse_args()



    scan(args[0],options)













