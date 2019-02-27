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
mycolors = ['orange','darkblue','limegreen','purple','pink','yellow','cyan','magenta','Brown','DarkGray']
colors = blues + greens + reds+ yellows + grays + axillary
colors = mycolors+short_colors
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
        def __init__(self,xLen=1,yLen=1, options=None):
		#self.options = options
		#self.color_key = color_key 
	        sns.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'cornflowerblue'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
                self.fig.set_facecolor('skyblue') 
                self.fig.patch.set_facecolor('skyblue')
                matplotlib.rcParams['savefig.facecolor'] = 'skyblue'
                matplotlib.rcParams['ytick.labelsize'] = 7.5
                matplotlib.rcParams['xtick.labelsize'] = 7.5
                
                seaborn.set(rc={'axes.facecolor':'lightgrey', 'figure.facecolor':'cornflowerblue'})

                self.fig.patch.set_facecolor('skyblue')


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

		ax.set_title(pheno+'  '+str(round(R_i[0],2)))
		ax.set_xlabel(str(S_i[1]))
		ax2.set_title(str(round(R_r[0],3))+'   '+str(pval)) 
		ax.set_ylabel('Raw Gene Expression') 
		ax2.set_ylabel('Resid Gene Expression') 
		ax2.set_xlabel(str(S_r[1]))
		self.xLoc += 1








	def add_bin_gene(self,pheno,pval,samples,attributes,Yi,Yr,bMock=1):
		
		opts = list(set(attributes)) 
		groups = {opts[j]: [i for i in range(len(attributes)) if attributes[i] == opts[j]] for j in range(len(opts))}

		groups = {g: k for (g,k) in groups.items() if len(k)>15 and g not in ['boo']}


		gY = {opt: [Yi[j] for j in K] for opt,K in groups.items()} 
		gR = {opt: [Yr[j] for j in K] for opt,K in groups.items()} 


		ax = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,2), rowspan = bMock, colspan = 1)
		ax2 = plt.subplot2grid((self.xLen ,self.yLen), (self.xLoc,3), rowspan = bMock, colspan = 1)
		self.add_data(gY,ax) 
		self.add_data(gR,ax2) 

		ax.set_title(pheno) 
		ax2.set_title(str(pval)) 
		self.xLoc += bMock



	def add_data(self,gkey,ax,gene=None):

		opts = gkey.keys() 
		vals = gkey.values() 

		colors = ['DarkSalmon','DarkSlateGray','DeepPink','DarkTurquoise','ForestGreen','Crimson','Plum','Chocolate','FireBrick']
		color_key = {opts[j]: colors[j] for j in range(len(opts))} 
		

		vals = [gkey[opt] if opt in gkey else [] for opt in opts] 
		colors = [color_key[opt] if opt in color_key else 'k' for opt in opts] 
		
		pos, width = [self.xOffset+0.20*x for x in range(len(opts))], 0.2

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
    			plt.scatter(xJ, yJ, c=clevel, alpha=0.7,s=4,zorder=9)


                for patch, clevel in zip(self.bp['boxes'], colors):
                        patch.set_facecolor(clevel) #cm.prism(clevel))
                        patch.set_alpha(0.2)

		
		ax.set_xticklabels([opt.split('~')[-1] for opt in opts],fontsize=10,rotation=-30) 
		ax.set_xlim([pos[0]-0.2,pos[-1]+0.2])

















def resolve(gene,isos,cnts):



	T,E = [],[]  
	
	cnt_key = dd(lambda: dd(int)) 
	uniq_key  = dd(lambda: dd(int)) 


	if len(gene.split('>')) > 1: 
		ix = ",".join([ix.split('=')[-1] for ix in isos])
		T=[['INTERGENE_'+gene+'_'+ix]]
		E=[['EN_'+ix]]
	
	else:
		for i in isos:
			iT,iE = [],[] 
			chars = {char: len(i.split(char)) for char in ['/','|','@',':','-',',']} 
			if chars['/'] > 1:   iT,iE = parse([ix.split('=')[-1] for ix in i.split('/')])
			elif chars['|'] > 1: iT,iE = parse(i.split('|'))
			else:
				if    chars['@'] == 2:   t,e = i.split('@') 
				elif  chars[':'] == 2: t,e = i.split(':') 
				elif  chars[','] + 1 == chars['-']:  t,e = '>'.join(('TN_'+i).split(',')),'EN_'+i 
				else:
					print 'wtf'
					sys.exit() 
				iT,iE = [t],[e] 
			T.append(sorted(iT))
			E.append(sorted(iE))


	T_sets = [sorted(list(set(t))) for t in T]
	T_tuples = sorted(list(set([tuple(t) for t in T_sets])))
	T_universe = sorted(list(set([a for b in T_tuples for a in b])))

	if len(T_tuples) == 1:
		c_sums = [sum([cnts[i][j] for i in range(len(cnts))]) for j in range(len(cnts[0]))]
		return T_universe,[",".join(['1' for tu in T_universe]) if cs > 0 else ",".join(['0' for tu in T_universe]) for cs in c_sums] 

	else:
		for j,TS in enumerate(T_sets):

			trs,exons,C = T[j],E[j], cnts[j] 
			for i,c in enumerate(C):
				cnt_key[",".join(TS)][i] += c 
			
		c_keys = cnt_key.keys() 
		maximizer = Maximizer(T_universe,c_keys) 
		for i in range(len(cnts[0])): 
		 	ci = [cnt_key[k][i] for k in c_keys]
			maximizer.add(i,ci)	
		
 		return maximizer.universe,maximizer.data



class ContributionCounter:
        #def __init__(self,cnt_key,color_key, options=None):
        def __init__(self,samples):

		self.samples = samples 
		self.sampleRange = range(len(self.samples))
#		self.w = open('ISOS_AND_FEATURE_ESTIMATES.cnts','w') 
		self.w = sys.stdout
		self.w.write('%-50s %100s\n' % ('---',"\t".join(samples)))

	def calculate(self,gene,key):

		
		self.gene,self.table, transcripts, self.intergene = gene, [dd(float) for s in self.samples], [] , None
		for k,vals in key.items():
			if k in ['TOTAL','EXONIC','FORMS']: 
				for i,f in enumerate(vals): self.table[i][k] = f  
			elif k == 'TRANSCRIPTS': transcripts = vals 
			elif k.split('_')[0] == 'INTERGENE':
				self.intergene=k 
				for i,f in enumerate(vals): self.table[i]['INTERGENE'] = f 
			elif k == 'EXPC':
				continue 
			else:
				print k,'huh'
				sys.exit() 			

	
		if 'EXPC' in key.keys() and 'TOTAL' in key.keys(): 
			for i in range(len(self.samples)):
				self.table[i]['EXONIC'] = key['EXPC'][i] * key['TOTAL'][i] 
			

	

		for i in range(len(self.samples)): 
			
			total,exonic,isos = self.table[i]['TOTAL'], self.table[i]['EXONIC'],self.table[i]['FORMS'] 


			if total - exonic > 0: self.table[i]['INTRONIC'] = total-exonic 
			if exonic > 0 and len(transcripts)>0:		 
				for t,iso in zip(transcripts,isos): 
					self.table[i][t] = exonic*iso 
		
	
		return self	
	
	def write(self):
		keys = [k for k in sorted(self.table[0].keys()) if k != 'FORMS'] 
		transcripts = {} 
		for k in keys: 
			vals = [int(round(self.table[i][k])) for i in self.sampleRange]
			if sum(vals) == 0: continue 
			valstr = '\t'.join([str(v) for v in vals])
			if k in ['INTRONIC','EXONIC','TOTAL']:	self.w.write('%-50s %20s\n' % (self.gene+'@'+k,valstr))
			elif k == 'INTERGENE':			self.w.write('%-50s %20s\n' % (self.intergene,valstr))
			else:					transcripts[self.gene+'@'+k] = valstr
			
		if len(transcripts)>1:
			for t,vs in transcripts.items():
				self.w.write('%-50s %20s\n' % (t,vs))

				  



		




 




def estimate_isos(stream,fout='iso_percs'):

	cnt_idxs = dd(lambda: {}) 
	

	key = dd(lambda: dd(list))

	INIT = True 
	b_list, c_list = [], [] 
	geneID = None 
	gene_key = {} 
#	w = open(fout+'.txt','w') 
#	w = sys.stdout
	for line in stream:

		line = line.split() 
		if line[0] == '---': 
			samples = line[2::] 
			contribs = ContributionCounter(samples) 
			continue 
		
		gene,info,cnt_data = line[0],line[1],line[2::] 

		#feature,cnts = line[0],[float(x) for x in line[1::]] 
		if gene != geneID: 
			if geneID != None:
				contribs.calculate(geneID,gene_key).write() 		
 
			
			gene_key = dd(bool) 
			geneID = gene 

		info,forms = info.split('=')[0],info.split('=')[-1] 
		
		if info == 'FORMS':
			if forms.split('_')[0].upper() != 'INTERGENE': 
				gene_key['TRANSCRIPTS'] = forms.split(',') 
				gene_key[info] = [[float(x) for x in X.split(',')] for X in cnt_data]
			else:	
				gene_key[forms] = [float(x) for x in cnt_data]
				
		else:
			if info == 'GENE': info = 'EXONIC'
			gene_key[info] = [float(x) for x in cnt_data] 	 




		







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


	##############################################################################
	######################   STEP 1: LOAD FILE  ################################## 
	##############################################################################

	# colors = ['black','darkblue','limegreen','olive','purple','pink','Chartreuse','brown','cyan','black','Brown','DarkGray']
	(options, args) = parser.parse_args()

	if len(args) == 1:
		stream = open(args[0]) 
	else: 
		stream = sys.stdin


	estimate_isos(stream) 

	sys.exit() 


	plt.suptitle(gene.split(';')[-1])
	plt.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.95,wspace=0.3,hspace=0.9)
	plt.show() 
