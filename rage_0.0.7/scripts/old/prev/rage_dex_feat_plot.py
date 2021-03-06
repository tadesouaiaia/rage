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
        def __init__(self,xLen=1,yLen=1, pheno=None,options=None):
		#self.options = options
		#self.color_key = color_key 

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


		self.xOffset = 0 

		self.cnts = [] 
		self.names = [] 
		self.genes = [] 
		self.xLen,self.yLen= xLen, yLen
		self.xLoc, self.yLoc =0,0



	def add_cont_gene(self,pheno,pval,samples,attributes,Yi,Yr):



		ax = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		ax2 = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,self.yLoc+1), rowspan = 1, colspan = 1)

		attD, rawD,resD = dd(list), dd(list), dd(list) 
		color_key = {'EB': 'red', 'ES': 'blue', 'AB': 'green', 'NA': 'k'} 		

		for i in range(len(samples)):

			s,a,yi,yr = samples[i], attributes[i], Yi[i], Yr[i]
			clr= 'k'

			key = 'NA' 

			if s[0:2] == 'EB': 
				clr = 'red'
				key = 'EB'
			elif s[0:2] == 'ES': 
				clr= 'cyan'
				key = 'ES' 
			elif s[0] == 'T':  
				clr = 'green'
				key = 'AB' 

			attD[key].append(a) 
			rawD[key].append(yi) 
			resD[key].append(yr) 
			

			ax.scatter(a,yi,color=clr,alpha=0.5) 
			ax2.scatter(a,yr,color=clr,alpha=0.5) 

		sR, iR, rR, pR, std_err = stats.linregress(attributes, Yi) 
		sS, iS, rS, pS, std_err = stats.linregress(attributes, Yr) 
		R_i,S_i = pearsonr(attributes,Yi), pearsonr(attributes,Yi) 
		R_r,S_r = pearsonr(attributes,Yr), pearsonr(attributes,Yr) 

		xEnd = ax.get_xlim()[1]*0.9
		for k in rawD.keys():

			if k == 'ES': continue  
			sR, iR, rR, pR, std_err = stats.linregress(attD[k],rawD[k]) 
			ax.plot([0,xEnd],[iR,iR+(sR*xEnd)],color=color_key[k])
			ax.text(xEnd,iR+(sR*xEnd),round(pR,4),color=color_key[k]) 	
			xEnd *= 1.2 

		ax.plot([0,xEnd],[iR,iR+(sR*xEnd)],color='k',linewidth=0.5)
	




		xEnd = ax2.get_xlim()[1]*0.9 
		for k in resD.keys(): 
			if k == 'ES': continue  
			sR, iR, rR, pR, std_err = stats.linregress(attD[k],resD[k]) 
			ax2.plot([0,xEnd],[iR,iR+(sR*xEnd)],color=color_key[k])
			ax2.text(xEnd,iR+(sR*xEnd),round(pR,4),color=color_key[k]) 	
			xEnd *= 1.2 

		sR, iR, rR, pR, std_err = stats.linregress(attributes, Yr) 
		ax.plot([0,xEnd],[iR,iR+(sR*xEnd)],color='k',linewidth=0.5)

		ax.set_title(pheno+' rawR= '+str(round(R_i[0],2)))
		#ax.set_xlabel(str(S_i[1]))
		ax2.set_title("resid r,pv= "+str(round(R_r[0],3))+'   '+str(pval)) 
		#ax.set_ylabel('Raw Gene Expression') 
		#ax2.set_ylabel('Resid Gene Expression') 
		#ax2.set_xlabel(str(S_r[1]))
		ax.set_yticks([]) 
		ax2.set_yticks([]) 
		self.xLoc += 1




	def setup_groups(self,attributes):
		opts = list(set(attributes)) 
		my_groups = {opts[j]: [i for i in range(len(attributes)) if attributes[i] == opts[j]] for j in range(len(opts))}
		groups = dd(list) 
		g_list = sorted([(len(k),k,g) for (g,k) in my_groups.items()],reverse=True) 
		for i,(kL,kG,g) in enumerate(g_list):
			
			if i < 3  or (i<5 and kL > 5): 
				groups[g] = kG 
			elif i < 10 and kL > 5: 
				groups[g] = kG 	
			else:
				groups['OTHER'].extend(kG) 

		return groups 

		



	def add_bin_gene(self,gene,pval,samples,attributes,Yi,Yr,key=None):
	


	

		groups = self.setup_groups(attributes)


		gY = {opt: [Yi[j] for j in K] for opt,K in groups.items()} 
		gR = {opt: [Yr[j] for j in K] for opt,K in groups.items()} 


		ax = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,0), rowspan = 1, colspan = 1)
		ax2 = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,1), rowspan = 1, colspan = 1)
		

		labels,items  = self.add_data(gY,ax) 
		labels,items  = self.add_data(gY,ax2) 
		if self.xLoc == 0: 
			ax.legend(items,labels,ncol=len(labels),loc='upper right',bbox_to_anchor=(1.5,1.8),fontsize=11,handletextpad=0.35,columnspacing=0.475)
		ax.axis('off') 
		ax2.axis('off') 
#		ax.set_yticks([]) 
#		ax2.set_yticks([]) 

		ax.set_title(gene.split(';')[-1]+' expression ',fontweight='bold') 
		ax2.set_title(gene.split(';')[-1]+' residuals   '+str(pval),fontweight='bold')
		ax.set_ylabel(str(pval)) 




		if key:
			yMe = 2
			for K  in sorted(key.keys()): 
				
				my_attributes = [key[K][s] for s in samples]
				my_attributes = ['MEGA' if x == '+' else x for x in my_attributes]
				my_groups = self.setup_groups(my_attributes) 	
				my_groups = {a: b for (a,b) in my_groups.items() if len(a.split('_')) == 1 and a not in ['OTHER','NA','CP']}
				#ax3 = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,2), rowspan = 1, colspan = 1)
				ax3 = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,yMe), rowspan = 1, colspan = 1)
				gY = {opt: [Yi[j] for j in K] for opt,K in my_groups.items()} 
				self.add_data(gY,ax3,True) 
				ax3.set_title(K)
				yMin = ax.get_ylim()[0] 
				mg = gY.keys() 
				for j,x in enumerate(ax3.get_xticks()):
					ax3.text(x,yMin-1,mg[j],fontweight='bold',fontsize=10)
				ax3.axis('off')  
				yMe += 1

		self.xLoc += 1




	def add_data(self,gkey,ax,KEY=False):

		opts = gkey.keys() 
		vals = gkey.values() 

		colors =['DarkSalmon','DarkSlateGray','DeepPink','DarkTurquoise','ForestGreen','Crimson','Plum','Chocolate','FireBrick']
		colors = VWCOLORS

		axillary = ['DarkSalmon','DarkSlateGray','DeepPink','DarkTurquoise','ForestGreen','Crimson','Plum','Chocolate','FireBrick']

		if KEY:
			if len(opts) == 2: colors = ['magenta','lime'] 
			elif len(opts) == 3:  colors = ['green','cyan','red'] 
			else:  	colors =['olive','DarkSalmon','DarkSlateGray','Plum','Chocolate','FireBrick']+VWCOLORS
        	try:
           		color_key = {opts[j]: colors[j] for j in range(len(opts))} 
		except IndexError:
            		print "okay this is error", opts, len(opts), len(colors)
            		sys.exit()

		

		vals = [gkey[opt] if opt in gkey else [] for opt in opts] 
		colors = [color_key[opt] if opt in color_key else 'k' for opt in opts] 
		
		labels,items=[],[]
		for v,c in color_key.items():
			labels.append(v) 
			items.append(Rect((0,0),1,1,fc=c))


		pos, width = [self.xOffset+0.10*x for x in range(len(opts))], 0.075

		self.bp = ax.boxplot(vals,positions=pos, widths=width, patch_artist=True, showmeans=True,whis=1)
		clevels = np.linspace(0., 3., len(opts)+1)

		xticks = [] 	
		means = sorted([(np.mean(v),opt) for (v,opt) in zip(vals,opts) if opt[0:2] == 'si'])
		medians = sorted([(np.percentile(v,60),opt) for (v,opt) in zip(vals,opts) if opt[0:2] == 'si'])
		sorts = [(sorted(v,reverse=True),opt) for (v,opt) in zip(vals,opts) if opt[0:2] == 'si']

		for i,opt in enumerate(opts):
			if opt == '': continue 
			clr = color_key[opt] 
			self.bp['boxes'][i].set_edgecolor(clr) 
			self.bp['boxes'][i].set_linewidth(1) 
			plt.setp(self.bp['medians'][i], color=clr, linewidth=3)
			plt.setp(self.bp['means'][i], marker='h',markersize=9,markerfacecolor=clr) 
			plt.setp(self.bp['caps'][(i*2)+1], color=clr,linewidth=1) 
			plt.setp(self.bp['caps'][(i*2)], color=clr,linewidth=1) 
			plt.setp(self.bp['whiskers'][i*2], color=clr,linewidth=1) 
			plt.setp(self.bp['whiskers'][1+(i*2)], color=clr,linewidth=1) 	
			plt.setp(self.bp['fliers'][i], markerfacecolor=clr, markeredgecolor = clr, marker='s',markersize=2.0)		


		clevels = np.linspace(0., 1., len(vals))
		xJitters = [np.random.normal(pos[i],0.01,len(vals[i])) for i in range(len(vals))]
		yJitters = [[vals[i][j]*np.random.normal(1,0.005,len(vals[i]))[j] for j in range(len(vals[i]))] for i in range(len(vals))]

		for xJ, yJ, clevel in zip(xJitters, yJitters, colors):
    			ax.scatter(xJ, yJ, c=clevel, alpha=0.7,s=4,zorder=9)


                for patch, clevel in zip(self.bp['boxes'], colors):
                        patch.set_facecolor(clevel) #cm.prism(clevel))
                        patch.set_alpha(0.2)

		
		ax.set_xticklabels([opt.split('~')[-1] for opt in opts],fontsize=10,rotation=-30) 
		ax.set_xlim([pos[0]-0.2,pos[-1]+0.2])

		return labels, items






















def read_key(key_file):
	my_key = dd(lambda: {}) 
	for line in open(key_file): 
		line = line.split() 
		if line[0] == '---':
			header = line
			#valid_opts = ['CELL_KIND', 'DETAIL_LAYER', 'CELL_TYPE_DETAIL', 'EB_DETAIL', 'MEGASIZE', 'MEA', 'BIOPSY_ID', 'GW']
			valid_opts = [ 'LAYER', 'MEA','CELLSPAN','GeW'] #  'GW']

		else:
			name = line[0] 
			for i in range(1,len(line)): 
				if header[i] in valid_opts: 
					if header[i] == 'CELLSPAN': 
						my_key['CELL-GROUP'][name] = line[i] 
					else:   
						my_key[header[i]][name] = line[i] 
	return my_key 

 




def load_cnts(stream,key_file,divisor=1.0,LOG=False,OBS = 0.22,pv=0.01):
	cnt_key = dd(lambda: {}) 
	cnt_idxs = dd(lambda: {}) 
	
	
	key = dd(lambda: dd(list))

	INIT = True 
	k, x_list = 0,[] 
	
	for line in stream:
		line = line.split() 
		if INIT: 
			feat = line[1]
			if line[3] == 'BIN': BIN = True
			else: 		     BIN = False 	
			INIT = False 

		if line[1] == feat:
			gene = line[0] 
			pval = float(line[2]) 
			sample_names= line[4].split(',') 
			raw, resids = [float(x) for x in line[6].split(',')],[float(x) for x in line[7].split(',')]
			obs=len([x for x in raw if x > 0])/float(len(raw))

			if obs < OBS: continue 
			if pval > pv: continue 
			if BIN: 
				attributes = line[5].split(',')


				x_list.append([pval,gene,sample_names,attributes,raw,resids])
			
			if line[3] == 'CONT': 
				attributes = [float(x) for x in line[5].split(',')]
				x_list.append([pval,gene,sample_names,attributes,raw,resids])

			if len(x_list) > 20: break
		else:
			break
	return feat,BIN,sorted(x_list)





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

    	pheno, BIN, x_list  = load_cnts(stream,options)


	if options.key: 
		key = read_key(options.key) 
	else:   key = None 


	if len(x_list) > 10: x_list = x_list[0:10] 
	if len(x_list) > 5: x_list = x_list[0:5]
	
	my_plot = plot(len(x_list),6,pheno,options) 
	for pval,gene,sample_names,attributes,raw,resids in x_list:
		if BIN: 
			my_plot.add_bin_gene(gene,pval,sample_names,attributes,raw,resids,key) 
		else:
			my_plot.add_cont_gene(gene,pval,sample_names,attributes,raw,resids,key) 






	plt.suptitle(pheno,fontweight='bold',fontsize=30) 
	plt.subplots_adjust(left=0.0, bottom=0.06, right=0.95, top=0.90,wspace=-0.25,hspace=1.0)
	plt.show() 
