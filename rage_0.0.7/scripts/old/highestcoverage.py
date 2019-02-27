#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 






def run_script(data_file,options):

	k = 0

	count_list = [] 
	obs_dictionary = {}  

	for line in open(data_file): 
		linelist = line.split() 
		if linelist[0] == '---': continue 
		gene = linelist[0] 
		gene_counts = [int(c) for c in linelist[1::]] 
#		gene_cnts = linelist[1:
		number_of_samples = float(len(gene_counts)) 
		gene_obs = len([x for x in gene_counts if x > 0]) / number_of_samples 
		gene_data = [sum(gene_counts), gene] 

		if gene_obs > 0.5: continue 

		obs_dictionary[gene] = gene_obs 
		count_list.append(gene_data) 
#		print count_list
		k+=1 
#		if k == 10: break 

	count_list.sort(reverse=True) 

	print 'here are your top genes' 
	for i in range(options.top): 

		my_data = count_list[i] 
		my_counts = my_data[0] 
		my_name   = my_data[1] 
		my_obs    = obs_dictionary[my_name] 


		print 'rank', i+1, my_name, 'has', my_counts, 'total reads and', my_obs,'observation rate' 






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-c", "--counts", default = None, type='string', help="Output Filename Prefix")

	parser.add_option("-t", "--top", default = 3, type=int, help="top how many")
	(options, args) = parser.parse_args()



	run_script(args[0],options)	














