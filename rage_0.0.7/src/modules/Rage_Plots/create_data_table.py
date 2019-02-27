#!/usr/bin/env python

from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import Circle as Circ
from collections import defaultdict as dd 
from collections import Counter as cc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 

import matplotlib.patches as mpatches


def load_key(key_file): 

	key = dd(lambda: {}) 
	cell_key = dd(lambda: dd(list)) 
	cell_counts = dd(int) 
	for line in open(key_file): 
		line = line.split() 
		if line[0] == '---': headers = line
		else:
			for i in range(1,len(line)):
				try: 
					key[headers[i]][line[0]] = float(line[i])
				except ValueError: 
					key[headers[i]][line[0]] = line[i]
			
			for i in range(2,len(line)):

				cell_key[headers[i]][line[1]].append(line[i])
				
			cell_counts[line[1]]+=1
	return key,cell_key,cell_counts



if __name__ == '__main__':

	import sys
	import os
	from optparse import OptionParser

	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)

	parser.add_option('-n',  dest= "name", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('-d',  dest= "dir", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('-o',  dest= "output", type = 'str' , default = sys.stdout, help = "horizontal data")
	parser.add_option('--show',action='store_true', dest='show',default=False,help='center means?')
	parser.add_option('--plot',action='store_true', dest='plot',default=False,help='center means?')
	parser.add_option('--results',action='store_true', dest='results',default=False,help='center means?')
	parser.add_option('--phase',action='store_true', dest='phase',default=False,help='center means?')
	parser.add_option('--key',  dest= "key", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--title',  dest= "title", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--opts',  dest= "opts", type = 'str' , default = None, help = "horizontal data")
	(options, args) = parser.parse_args()


	key,cell_key,cell_counts = load_key(options.key) 

	cell_types = ['SVZ','SVZ+','IZ','IZ+','SP','SP+'] 
	
	table_columns = ['Cells','Avg Reads','Avg Genes','Avg GW','Patch-Clamp'] 
	
	cell_data = [] 
	for c in cell_types:
	
		cLen =  float(cell_counts[c])
		tR, gO, gW  = [float(x) for x in cell_key['TOTAL_READS'][c]], [float(x) for x in cell_key['GENES_OBSERVED'][c]], [float(x) for x in cell_key['GW'][c]]
		cell_data.append([cell_counts[c],round(sum(tR)/cLen,1), round(sum(gO)/cLen,1), round(sum(gW)/cLen,1),cc(cell_key['EFIZ'][c])['YES']])

		
 


	# Get some pastel shades for the colors
	colors = plt.cm.BuPu(np.linspace(0, 0.5, len(cell_types)))
	n_rows = len(cell_data)

	index = np.arange(len(table_columns)) + 0.3
#	bar_width = 0.4

	# Initialize the vertical-offset for the stacked bar chart.
#	y_offset = np.zeros(len(columns))

	# Plot bars and create text labels for the table
	#for row in range(n_rows):
	 #   y_offset = y_offset + data[row]
	    
	  #  cell_text.append(['%1.1f' % (x / 1000.0) for x in y_offset])
	# Reverse colors and text labels to display the last value at the top.
	colors = colors[::-1]
	#cell_text.reverse()

	#print cell_text

	#fig = plt.figure(1)

#Table - Main table
	ax = plt.subplot2grid((2,2), (0,0), colspan=2, rowspan=2)
	# Add a table at the bottom of the axes
	the_table = ax.table(cellText=cell_data,
			      rowLabels=cell_types,
			      rowColours=colors,
			      colLabels=table_columns,
  			      colWidths = [0.15, 0.15,0.15,0.15,0.15],
  			      #rowHeights = [0.15, 0.15,0.15,0.15,0.15],
			      loc='center')
			      #fontsize=55)
	the_table.auto_set_font_size(False)
	the_table.set_fontsize(40.5)
#	the_table.set_fontsize(18)
	the_table.scale(1.6,1.3)
	ax.axis("off") 
	plt.subplots_adjust(left=0.2,top=1.8, bottom=-0.8,right=0.9,wspace=0.01,hspace=0.01)
	# Adjust layout to make room for the table:
	#plt.subplots_adjust(left=0.2, bottom=0.2)

#	plt.ylabel("Loss in ${0}'s".format(value_increment))
#	plt.yticks(values * value_increment, ['%d' % val for val in values])
	plt.xticks([])
	#plt.title('Loss by Disaster')

	plt.show()



