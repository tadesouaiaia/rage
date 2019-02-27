#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 


fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}




def scan(cnts_file):
	info = dd(list) 
	samples = [] 
	for line in open(cnts_file):
		line = line.split()

		
		sample,desc = line[0],line[1].split('/') 
		for i in range(len(desc)):
			info[desc[i]].append(float(line[2+i]))
		samples.append(sample)
	
	totals = info['True'] 

	
	DS_TESTS = [25000,50000,100000]


	for i in range(len(samples)):

		t = info['True'][i]
		sI = info['Infer'][i]
		sIM = info['Infer-Missing'][i]
		sP  = info['Impute'][i]
		sPM = info['Impute-Missing'][i]

		sIO = sIM - t 
		sPO = sPM - t 

		sI_onlyPresent = sI - sIO 
		sP_onlyPresent = sP - sPO

		info['Infer-Relative'].append((sI-sIO))
		info['Impute-Relative'].append((sP-sPO))  

		


	totals = info['True'] 



	for k in info.keys():
		ds_vals = dd(list)
		up_samples = dd(list) 

		if k == 'True': continue
		if k.split('-')[-1] == 'Relative':
			X = [info[k][i] for i in range(len(samples))]
			RT = [info[k.split('-')[0]][i] for i in range(len(samples))]
			percs = [info[k][i]/RT[i] for i in range(len(samples))]
		elif k.split('-')[-1] == 'Missing':  
			X = [info[k][i] for i in range(len(samples))]
			percs = [totals[i]/X[i] for i in range(len(samples))]

		else:
			continue 
			print 'hmm',k
			continue


		for T in DS_TESTS:
			w=open('RAGE_DS-'+k+'.vals','w')
			if k.split('-')[-1] == 'Missing': 
				for s,t,x,p in zip(samples,totals,X,percs):
					ds_val = p*T
					ds_vals[T].append(ds_val) 
					if ds_val > t:
						up_samples[T].append(s) 


			elif k.split('-')[-1] == 'Relative': 
				for s,t,x,p in zip(samples,totals,X,percs):
					ds_val = p*T
					ds_vals[T].append(ds_val) 
					if ds_val > t:
						up_samples[T].append(s) 
					
				

			else:
				continue 
				print 'hmmm'
				continue 
				for s,t,x,p in zip(samples,totals,X,percs):
					if p > 1:	ds_val = T
					else:		ds_val = p*T
					
					ds_vals[T].append(ds_val) 
					if ds_val > t:
						up_samples[T].append(s) 
						
					




		for i in range(len(samples)):

			w.write('%s %d %5.2f %5.5f %s ' % (samples[i],totals[i],X[i],percs[i],k))

			for T in DS_TESTS:
				w.write(' %s %s %s' % (T,(samples[i] in up_samples[T]),int(ds_vals[T][i])))
			w.write('\n')











if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"

	parser = OptionParser(usage=usage)


	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	(options, args) = parser.parse_args()



	scan(args[0])













