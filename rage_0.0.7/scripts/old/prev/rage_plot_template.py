#!/usr/bin/env python

import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc

import numpy as np
#from scipy.stats import spearmanr as spearmanr
import seaborn
from math import log
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
from matplotlib.lines import Line2D









def plot_example(): 

	ax1 = plt.subplot2grid((2,2), (0,0))
	ax2 = plt.subplot2grid((2,2), (0,1))
	ax3 = plt.subplot2grid((2,2), (1,0))
	ax4 = plt.subplot2grid((2,2), (1,1))

	days_of_week   = ['Monday','Tuesday','Wednesday','Thursday','Friday'] 
	anxiety_levels = [3,4,6,9,3] 

	ax1.plot([0,1,2,3,4],anxiety_levels) 
	ax1.set_xticks([0,1,2,3,4])
	ax1.set_xticklabels(days_of_week) 


	for i in range(len(anxiety_levels)):
		ax2.scatter(i,anxiety_levels[i],color='blue',marker='*',s=50)
		ax2.set_xticks([0,1,2,3,4])
		ax2.set_xticklabels(days_of_week) 

	for i in range(len(anxiety_levels)):
		ax3.bar(i,anxiety_levels[i],color='green') 
		ax3.set_xticks([0,1,2,3,4])
		ax3.set_xticklabels(days_of_week) 
	
	anxiety_fracs = [anxiety_level / float(sum(anxiety_levels)) for anxiety_level in anxiety_levels] 

	ax4.pie(anxiety_fracs, labels = days_of_week) 

	plt.suptitle('Multiple Plots') 

	
	plt.show() 





def scan_file(input_file,options):
	for line in open(input_file): 
		line = line.split() 
		print line





		


				






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	plot_example() 
#	scan_key(args[0],options)













