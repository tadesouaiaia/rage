#!/usr/bin/env python

import os
import sys
import matplotlib
import seaborn as sns 
#matplotlib.use('Agg')
import pylab as plt
import numpy as np
import math
import matplotlib.lines as mlines
from scipy.optimize import curve_fit

from scipy.stats import pearsonr 
from scipy.stats import spearmanr
from collections import Counter as cc
from matplotlib import cm 
from collections import defaultdict as dd
from matplotlib.lines import Line2D as Line
from matplotlib import cm as cm
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Rectangle as Rect
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib._png import read_png
from PIL import Image
from math import exp
from math import sin
from math import pi
import matplotlib
#import seaborn
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
from matplotlib.lines import Line2D as Line
import seaborn
# Generate the image
import copy
import sys
import os
import pickle
from collections import defaultdict as dd 
from collections import defaultdict as dd 
from datetime import datetime
from math import fabs
from math import log
import numpy as np 
import pylab 
import scipy.stats as stats
import os
import sys
import numpy as np
from scipy.interpolate import spline
from scipy.interpolate import spline
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.collections as mcoll
from matplotlib.collections import LineCollection
#!/usr/bin/env python

import matplotlib
#import seaborn
import matplotlib.pyplot as plt
import seaborn as sns 
# Generate the image
import copy
import sys
import os
import pickle
from collections import defaultdict as dd 
from collections import defaultdict as dd 
from collections import Counter as cc 
from datetime import datetime
from math import fabs
#!/usr/bin/env python
from math import log
import numpy as np 
import pylab 
import scipy.stats as scistats
import os
import sys
import math 
import matplotlib.pyplot as plt
from scipy.stats import chi2
from collections import defaultdict as dd 
from random import random
from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import FancyArrowPatch as Arrow

#from matplotlib.colors import ListedColormap, BoundaryNorm

matplotlib.rcParams['xtick.labelsize'] = 14


blues,greens,reds,yellows,grays = ['blue','cyan','darkblue'],['green','limegreen','Olive'],['red','magenta','purple'],['yellow','Chartreuse','goldenrod'],['black','Brown','DarkGray']
axillary = ['DarkSalmon','DarkSlateGray','DeepPink','DarkTurquoise','ForestGreen','Crimson','Plum','Chocolate','FireBrick']
short_colors = ['r','b','g','y','k','R','B','G','R','Y']
mycolors = ['orange','darkblue','limegreen','purple','red','yellow','cyan','magenta','Brown','DarkGray']
colors = blues + greens + reds+ yellows + grays + axillary
VWCOLORS = mycolors+short_colors+colors
colors5 = ['blue','DarkGray','DarkGray','DarkGray','DarkGray']
colors10 = ['DarkGray','DarkGray','DarkGray','DarkGray','lime']
colors15 = ['DarkGray','orange','purple','green','DarkGray']
colors20 = ['DarkGray','DarkGray','DarkGray','cyan','DarkGray','red','DarkGray']
colors = colors5 + colors10+colors15+colors20
##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

blues,greens,reds,yellows,grays = ['blue','cyan','darkblue'],['green','limegreen','Olive'],['red','magenta','purple'],['yellow','Chartreuse','goldenrod'],['black','Brown','DarkGray']
axillary = ['DarkSalmon','DarkSlateGray','DeepPink','DarkTurquoise','ForestGreen','Crimson','Plum','Chocolate','FireBrick']
short_colors = ['r','b','g','y','k','R','B','G','R','Y']
circles,squares,triangles,others = ["o",".","8",'|'], ['s','p','D','_'],['>',"<","v","^"],['*','H','x',',']
flat_colors = [item for sublist in [[x[i] for x in [blues,greens,yellows,reds,grays]] for i in range(3)] for item in sublist]
flat_marks  = ['.']+[item for sublist in [[x[i] for x in [circles,triangles,squares,others]] for i in range(4)] for item in sublist]
flat_colors = flat_colors+axillary+flat_colors+flat_colors+short_colors
flat_colors = ['black']+flat_colors+flat_colors+flat_colors+flat_colors
flat_marks  = flat_marks+flat_marks+flat_marks
flat_sizes  = [30,50,20,20,30,20,40,30]

def select(x,lab='clr'):
	if not x: 
		return False
	elif len(x) == 1: 
		x = x[0]
		if type(x) == int: 
			if lab=='clr': return flat_colors[x]
			if lab=='mrk': return flat_marks[x]
			if lab=='size': return flat_sizes[x]
		return x
	else:
		print 'hmm'



def plot_error(msg):
	sys.stderr.write('plotError: '+msg+'\n')
	sys.exit()







class plot:
        #def __init__(self,cnt_key,color_key, options=None):
        def __init__(self,key,xLen=1,yLen=1, pheno=None,options=None):
		#self.options = options
		#self.color_key = color_key 

		self.key = key 
		self.pheno = pheno 
	        sns.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'cornflowerblue'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
                self.fig.set_facecolor('skyblue') 
                self.fig.patch.set_facecolor('skyblue')
                matplotlib.rcParams['savefig.facecolor'] = 'lightgrey'
                matplotlib.rcParams['ytick.labelsize'] = 7.5
                matplotlib.rcParams['xtick.labelsize'] = 7.5
                seaborn.set(rc={'axes.facecolor':'lightgrey', 'figure.facecolor':'lightgrey'})
                self.fig.patch.set_facecolor('lightgrey')
		self.xLen,self.yLen= xLen, yLen
		self.xLoc, self.yLoc =0,0






	def add_scatter(self,A,B):

		X=self.key[A]
		Y=self.key[B]

		ax = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		for x,y in zip(X,Y):
			ax.scatter(x,y) 

		R = pearsonr(X,Y) 


		ax.set_title('R= '+str(round(R[0],4)))
		ax.set_xlabel(A) 
		ax.set_ylabel(B)
		 
		self.yLoc += 1

		if self.yLoc == self.yLen:
			self.xLoc +=1
			self.yLoc =0


 




def load_cnts(stream,options,LOG=False,OBS = 0.22,pv=0.01):
	
	
	n=0
	key = dd(lambda: {}) 

	for line in open(options.key): 
		line = line.split() 	
		if line[0] == '---': header = line 
		else:
			for i in range(1,len(line)):
				try:
					key[header[i]][line[0]] = float(line[i])
				except ValueError:
					key[header[i]][line[0]] = line[i]

	for line in stream:
		line = line.split() 
		if line[0] == '---': 
			header = line 
			pg_idx = {opt: [i for i in range(1,len(line)) if (line[i] in key['pgGRP'].keys() and key['pgGRP'][line[i]] == opt)] for opt in list(set(key['pgGRP'].values()))}
			c_key = {} 


			for k in key.keys():
				if k != 'pgGRP': 
					c_key[k] = [key[k][line[i]] for i in range(1,len(line)) if line[i] in key[k].keys()]


			my_keys = c_key.keys()
			for a in range(len(my_keys)): 
				for b in range(a+1,len(my_keys)):
					k1,k2 = my_keys[a],my_keys[b] 
					R,pv = pearsonr(c_key[k1],c_key[k2])
					#if k1 == 'OBS':	
						#print k1,k2,R,pv 

			print '--- obsR','|','pgR','obsR','totR','VR_obs','VR_tot','|','pgPV','obsPV','totPV','VR_obsPV','VR_totPV','|','OBS_ORDER'

		else:
			gene,cnts = line[0],[log(float(x)+1.0,2) for x in line[1::]] 
			pg_cnts = {opt: [cnts[i-1] for i in v] for opt,v in pg_idx.items()}

			c_obs = [1 if c > 0 else 0 for c in cnts] 
			pVal,pObs,pTupes = [], [], []  
			for k,V in pg_cnts.items():
				pZ = [1 if v>0 else 0 for v in V]
				pVal.extend([float(k.split('pg')[0]) for v in V])
				pObs.extend([1 if v>0 else 0 for v in V])	
				pTupes.append((sum(pZ)/float(len(pZ)),k))

			pTupes.sort() 
			R_pg,PV_pg =  pearsonr(pVal,pObs) 
			R_obs,PV_obs =  pearsonr(c_key['OBS'],c_obs) 
			R_tot,PV_tot =  pearsonr(c_key['LOG_TOTAL'],c_obs) 
		
			R_cobs,PV_cobs =  pearsonr(c_key['OBS'],cnts) 
			R_ctot,PV_ctot =  pearsonr(c_key['LOG_TOTAL'],cnts) 

	

			print gene, round(sum(c_obs)/float(len(cnts)),3), "|",round(R_pg,4), round(R_obs,4),round(R_tot,4),round(R_cobs,4),round(R_ctot,4),"|",
			print round(PV_pg,4), round(PV_obs,4),round(PV_tot,4),round(PV_cobs,4),round(PV_ctot,4),
			print "|", " ".join([z[1] for z in pTupes])


			 
			
		





	sys.exit() 
	my_plot = plot(key,2,2)

	my_plot.add_scatter('OBS','RSQ')  
	my_plot.add_scatter('AVG','RSQ')
	my_plot.add_scatter('CV','RSQ')
	



	plt.show() 



#	AVG CV 0.635635
#	AVG OBS 0.884303
#	AVG RATE 0.88418
#	AVG RSQ 0.379864
#	CV OBS 0.839653
#	CV RATE 0.839553
#	CV RSQ 0.32863
#	OBS RATE 0.999858
#	OBS RSQ 0.407585
#	RATE RSQ 0.407627








if __name__ == '__main__':

	import sys
	import os
	from optparse import OptionParser

	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)



	parser.add_option('-n',  dest= "name", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--key',  dest= "key", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--bkey',  dest= "bkey", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('-c','--cnts',  dest= "cnts", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('-b','--bcnts',  dest= "bcnts", type = 'str' , default = None, help = "horizontal data")


	##############################################################################
	######################   STEP 1: LOAD FILE  ################################## 
	##############################################################################

	# colors = ['black','darkblue','limegreen','olive','purple','pink','Chartreuse','brown','cyan','black','Brown','DarkGray']
	(options, args) = parser.parse_args()

	if len(args) == 1:
		stream = open(args[0]) 
	else: 
		stream = sys.stdin

    	pheno  = load_cnts(stream,options)


	sys.exit()



	plt.suptitle(pheno,fontweight='bold',fontsize=30) 
	plt.subplots_adjust(left=0.0, bottom=0.06, right=0.95, top=0.90,wspace=-0.25,hspace=1.0)
	plt.show() 
