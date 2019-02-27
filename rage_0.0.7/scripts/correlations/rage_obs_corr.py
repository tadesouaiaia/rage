#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 


fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}




def scan(cnts_file,options):
	cnts = []
	genes = [] 
	for line in open(cnts_file):
		line = line.split()
		if line[0] == '---':
            		names,totals = line[1::],[[] for i in range(len(line)-1)]

        	else:


            		gene,Tcnts = line[0],[float(x) for x in line[1::]]
			for i in range(len(Tcnts)):
				totals[i].append((Tcnts[i])) 
			cnts.append(Tcnts) 
			genes.append(gene) 

	sObs,sTot,slogTot = [],[],[]
	for j,T in enumerate(totals):
		ts = [x for x in T if x > 0] 
		sObs.append(len(ts)) 
		sTot.append(sum(ts))
		slogTot.append((log(sum(ts)+1.0,2)))
	skey = {} 
	skey['obs'],skey['tot'],skey['log_tot']  = sObs, sTot, slogTot
	for g,c in zip(genes,cnts):

		cL = [log(x+1,2) for x in c]
		rk = {} 
		lk = {}
		for k in ['obs','tot','log_tot']: 

			R,rpv = stats.pearsonr(skey[k],c) 
			S,spv = stats.spearmanr(skey[k],c) 
			
			rk[k] = [round(R,3),round(rpv,4),round(S,3),round(spv,4)]
			
			R,rpv = stats.pearsonr(skey[k],cL) 
			S,spv = stats.spearmanr(skey[k],cL) 
			lk[k] = [round(R,3),round(rpv,4),round(S,3),round(spv,4)]
		
		print g,'RAW',
		for k,v in rk.items(): 
			print k,v[0],v[2],v[1],v[3],"|",
		print 'LOG',
		for k,v in lk.items(): 
			print k,v[0],v[2],v[1],v[3],"|",


		print "" 







if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"

	parser = OptionParser(usage=usage)


	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	(options, args) = parser.parse_args()



	scan(args[0],options)













