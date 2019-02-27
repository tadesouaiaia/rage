#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 






def run_script(data_file,options):
	k=0
	for line in open(data_file):
 		lp = line.strip() 
		line = line.split() 
		if line[0] == '---': 
			header = lp 
			gk = dd(lambda: dd(lambda: [0 for i in range(len(line)-1)])) 
		else:

			iso, cnts = line[0], [int(x) for x in line[1::]] 
			if len(iso.split('=')) != 2: continue 
			g,i = iso.split('=') 
			if 'T' not in i: continue 
			trans = ",".join(sorted([t.split('@')[0].split(':')[0] for t in i.split('|')])) 
			
			for j in range(len(cnts)): 
				gk[g][trans][j] += cnts[j] 

			k+=1
#			if k > 100: break 

	print header
	for g in gk: 
		if len(gk[g]) == 1: continue 
		for t in gk[g]: 
			t_cnts = gk[g][t] 
			if len([x for x in t_cnts if x > 0])/ float(len(t_cnts)) < 0.05: continue 
			if len([x for x in t_cnts if x > 5])/ float(len(t_cnts)) < 0.025: continue 
			
			print g+'='+t, ' '.join([str(s) for s in t_cnts]) 

			#print g+'='+t," ".join([str(s) for s in gk[g][t]])



if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args[0],options)	














