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
from random import shuffle 
#import statsmodels.sandbox.stats.multicomp.fdrcorrection0 as fdr

#statsmodels.sandbox.stats.multicomp.multipletests(p,alpha=0.05,method='fdr_bh')
#statsmodels.sandbox.stats.multicomp.fdrcorrection0(p,alpha=0.05)


#def parse_out_file(line): 

#			--- RS CV obs len parent maxG maxMeans maxObs maxChi maxP | params

#			['ENSG00000000971;chr1;CFH', '0.007', '3.537', '276', '2455', 'FULLGROUP', 'FULLGROUP~ES', '0.34', '0.19', '2.6e-05', '1.3e-03', '|', 'FULLGROUP~AB', 'False', '2.51e-03', '2.51e-03', '|', 'FULLGROUP~EB', 'False', '3.00e-05', '3.00e-05']



def run_script(data_file,options):

	grp_key = dd(list) 
	h_key = dd(lambda: {})  	
	samples = [] 
	for line in open(options.key):
		line = line.split() 
		if line[0] == '---': headers = line
		else:
			grp_key[line[5]].append(line[0]) 
			
			for i in range(1,len(line)):
				h_key[line[0]][headers[i]] = line[i] 
			samples.append(line[0]) 

	gene_key = dd(lambda: dd(list)) 
	genes = [] 
	X_samples, Y_samples = [], [] 
	X_names, Y_names = [], [] 
	for line in open(data_file):
		line  = line.split() 
		if line[0] == '---': 
			for i,s in enumerate(line):
				if s == '---': continue 
				elif h_key[s]['S_GROUP'] == 'X': 
					X_samples.append(i) 
					X_names.append(s) 
				elif h_key[s]['S_GROUP'] == 'Y': 
					Y_samples.append(i) 
					Y_names.append(s) 
				else:
					print 'wtf'

		else:
			gene,cnts = line[0],['---']+[int(x) for x in line[1::]] 
			genes.append(gene) 
			X_cnts = [cnts[z] for z in X_samples]
			Y_cnts = [cnts[z] for z in Y_samples]
			gene_key[gene]['X'] = X_cnts 
			gene_key[gene]['Y'] = Y_cnts 

	X_idxs = range(len(genes)) 
	Y_idxs = range(len(genes)) 
	shuffle(Y_idxs)

	print '---'," ".join(X_names)," ".join(Y_names) 
			 
	 
	for x,y in zip(X_idxs,Y_idxs):
		gX,gY = genes[x],genes[y] 
		cX,cY = gene_key[gX]['X'], gene_key[gY]['Y'] 
		xM,yM =  np.mean(cX),np.mean(cY) 
		if xM >= yM: 	    FC = xM/yM
		else:               FC = -1*yM/xM

		new_gene =  gX+'_'+gY.split(';')[-1]+'_'+str(x)+'_'+str(y)+'_'+str(FC)

		print new_gene," ".join([str(c) for c in cX])," ".join([str(c) for c in cY])

	sys.exit() 
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














