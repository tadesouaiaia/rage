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
	INIT = True
	for line in open(cnts_file):
		line = line.split() 
		obs_key, avg_key = {},{} 
		for i in range(len(line)): 

			if line[i] == 'obs':
				obs_loc = i  
				obs_vals = line[i+1].split(',') 
				obs_names = line[i+2].split(',')
				for j in range(len(obs_names)): obs_key[obs_names[j]] = obs_vals[j]
			if line[i] == 'avg': 
				avg_vals = line[i+1].split(',') 
				avg_names = line[i+2].split(',')
				for j in range(len(avg_names)): avg_key[obs_names[j]] = avg_vals[j]
				avg_loc = i+3
			if line[i] == 'FC': 
				fc1,fc2 = line[i+1],line[i+2] 
			if line[i] == 'CHI': 
				chi1,chi2 = line[i+1],line[i+2] 

		if INIT: 
			names = sorted(obs_names)
			INIT = False 	
			print '---  ID MODEL FDR-VAL',"OBS:"," ".join(names),"avg"," ".join(names),'FC1 FC2 CHI1 CHI2'
			

		print " ".join(line[0:obs_loc]),'|',' '.join([obs_key[n] for n in names]),'|'," ".join([avg_key[n] for n in names]),fc1,fc2,chi1,chi2 
		 







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













