#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 


fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}




def scan(cnts_file,options):
	for line in open(cnts_file):
		line = line.split()
		if line[0] == '---':
            		names,totals = line[1::],[[] for i in range(len(line)-1)]

        	else:


            		gene,Tcnts = line[0],[float(x) for x in line[1::]]
			for i in range(len(Tcnts)):
				totals[i].append((Tcnts[i],gene)) 

	
	for j,T in enumerate(totals):
		ts = [x for x in sorted(T,reverse=True) if x[0] > 0]


		tlen = len(ts) 
		for i,t in enumerate(ts):
			
			print t[1],t[0],i,"|",tlen,names[j]
			











if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"

	parser = OptionParser(usage=usage)


	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	(options, args) = parser.parse_args()



	scan(args[0],options)













