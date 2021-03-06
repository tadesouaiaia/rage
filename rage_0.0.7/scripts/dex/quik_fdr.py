#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as mlt 
import statsmodels.sandbox.stats.multicomp as mpt 
#import statsmodels.sandbox.stats.multicomp.fdrcorrection0 as fdr

#statsmodels.sandbox.stats.multicomp.multipletests(p,alpha=0.05,method='fdr_bh')
#statsmodels.sandbox.stats.multicomp.fdrcorrection0(p,alpha=0.05)


#def parse_out_file(line): 

#			--- RS CV obs len parent maxG maxMeans maxObs maxChi maxP | params

#			['ENSG00000000971;chr1;CFH', '0.007', '3.537', '276', '2455', 'FULLGROUP', 'FULLGROUP~ES', '0.34', '0.19', '2.6e-05', '1.3e-03', '|', 'FULLGROUP~AB', 'False', '2.51e-03', '2.51e-03', '|', 'FULLGROUP~EB', 'False', '3.00e-05', '3.00e-05']



def run_script(data_file,options):
	k=0
	genes, p_vals, gene_data = [], [] , []
	overs, unders,obs = [], [] , {} 
	pm_under, pn_under, pm_over, pn_over = [], [], [], [] 
	lines = [] 
	pvs = [] 
	for line in open(data_file): 
		line = line.split() 
		lines.append(line) 
		if line[0] != '---': 
			pvs.append(float(line[-1])) 
	

	fdr_bools, fdr_pvs =  mlt.fdrcorrection(pvs) 
	for j,line in enumerate(lines): 
		if line[0] == '---': 
			line.extend(['fdr_bool','fdr_val'])
		else:
			line.extend([fdr_bools[j-1],fdr_pvs[j-1]])
		print " ".join([str(s) for s in line]) 	













if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args[0],options)	














