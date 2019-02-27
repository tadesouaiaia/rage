#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 
import numpy as np 

import warnings
warnings.filterwarnings("ignore")




SWAP_DICT={"EB1010": "EB891","EB1061": "EB1690","EB1069": "EB1202","EB1093": "EB1534","EB1202": "EB1069","EB1216": "EB1266","EB1266": "EB1216","EB1310": "EB1323","EB1323": "EB1310","EB1345": "EB1751","EB1373": "EB1379","EB1379": "EB1373","EB1418": "EB1730","EB1422": "EB459","EB1432": "EB1492","EB1447": "EB1473","EB1448": "EB1463","EB1463": "EB1448","EB1473": "EB1447","EB1492": "EB1432","EB1528": "EB1726","EB1534": "EB1093","EB1548": "EB1551","EB1551": "EB1548","EB1690": "EB1061","EB1726": "EB1528","EB1730": "EB1418","EB1751": "EB1345","EB459": "EB1422","EB891": "EB1010","EB899": "EB986","EB986": "EB899"}

def scan(cnts_file,options):
	genes = [] 
	for line in open(cnts_file):
		line = line.split()
		if line[0] == '---':
			names = line[1::]
			swap_names = [] 
			for n in names: 
				if n in SWAP_DICT: swap_names.append(SWAP_DICT[n])
				else:              swap_names.append(n)
			print line[0]," ".join(swap_names) 

		else:
			print " ".join(line) 			













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













