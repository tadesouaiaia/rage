#!/usr/bin/env python
from __future__ import division
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as mlt 
import statsmodels.sandbox.stats.multicomp as mpt 
import seaborn as sns
import seaborn
from math import log
import math
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from random import shuffle
from sklearn.cluster import KMeans      
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import MinMaxScaler
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle as Rect
from matplotlib.lines import Line2D
import os
import sys
from collections import defaultdict as dd
from collections import Counter as cc
from scipy import stats
from math import log 
import math
from scipy.stats import chisquare
import numpy as np 
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors




#import statsmodels.sandbox.stats.multicomp.fdrcorrection0 as fdr

#statsmodels.sandbox.stats.multicomp.multipletests(p,alpha=0.05,method='fdr_bh')
#statsmodels.sandbox.stats.multicomp.fdrcorrection0(p,alpha=0.05)


#def parse_out_file(line): 

#			--- RS CV obs len parent maxG maxMeans maxObs maxChi maxP | params

#			['ENSG00000000971;chr1;CFH', '0.007', '3.537', '276', '2455', 'FULLGROUP', 'FULLGROUP~ES', '0.34', '0.19', '2.6e-05', '1.3e-03', '|', 'FULLGROUP~AB', 'False', '2.51e-03', '2.51e-03', '|', 'FULLGROUP~EB', 'False', '3.00e-05', '3.00e-05']

X_PIXELS = 1392
Y_PIXELS = 1040



COLORS_1 = [ 'indigo', 'gold', 'hotpink', 'firebrick', 'indianred', 'sage', 'yellow', 'mistyrose', 'darkolivegreen', 'olive', 'darkseagreen', 'pink', 'tomato', 'lightcoral', 'orangered', 'navajowhite', 'lime', 'palegreen', 'darkslategrey', 'greenyellow', 'burlywood', 'seashell', 'mediumspringgreen', 'fuchsia', 'papayawhip', 'blanchedalmond', 'chartreuse', 'dimgray', 'black', 'peachpuff', 'springgreen', 'aquamarine', 'white', 'orange', 'lightsalmon', 'darkslategray', 'brown', 'ivory', 'dodgerblue', 'peru', 'darkgrey', 'lawngreen', 'chocolate', 'crimson', 'forestgreen', 'slateblue', 'lightseagreen', 'cyan', 'mintcream', 'silver', 'antiquewhite']

COLORS_2 = [ 'mediumorchid', 'skyblue', 'gray', 'darkturquoise', 'goldenrod', 'darkgreen', 'floralwhite', 'darkviolet', 'darkgray', 'moccasin', 'saddlebrown', 'grey', 'darkslateblue', 'lightskyblue', 'lightpink', 'mediumvioletred', 'slategrey', 'red', 'deeppink', 'limegreen', 'darkmagenta', 'palegoldenrod', 'plum', 'turquoise', 'lightgrey', 'lightgoldenrodyellow', 'darkgoldenrod', 'lavender', 'maroon', 'yellowgreen', 'sandybrown', 'thistle', 'violet', 'navy', 'magenta', 'dimgrey', 'tan', 'rosybrown', 'olivedrab', 'blue', 'lightblue', 'ghostwhite', 'honeydew', 'cornflowerblue', 'linen', 'darkblue', 'powderblue', 'seagreen', 'darkkhaki']

COLORS_3 = [ 'snow', 'sienna', 'mediumblue', 'royalblue', 'lightcyan', 'green', 'mediumpurple', 'midnightblue', 'cornsilk', 'paleturquoise', 'bisque', 'slategray', 'darkcyan', 'khaki', 'wheat', 'teal', 'darkorchid', 'deepskyblue', 'salmon', 'darkred', 'steelblue', 'palevioletred', 'lightslategray', 'aliceblue', 'lightslategrey', 'lightgreen', 'orchid', 'gainsboro', 'mediumseagreen', 'lightgray', 'mediumturquoise', 'darksage', 'lemonchiffon', 'cadetblue', 'lightyellow', 'lavenderblush', 'coral', 'purple', 'aqua', 'lightsage', 'whitesmoke', 'mediumslateblue', 'darkorange', 'mediumaquamarine', 'darksalmon', 'beige', 'blueviolet', 'azure', 'lightsteelblue', 'oldlace']

AVOID_1 = ['yellow', 'lime', 'lightblue', 'purple', 'k', 'cyan', 'grey', 'olive', 'red']

AVOID_COLORS = AVOID_1
COLORS = [c for c in COLORS_1+COLORS_2+COLORS_3 if c not in AVOID_COLORS] 







def parse_data(x):
	MICRONS_PER_PIXEL = 6.2
	MICRONS_PER_PIXEL = 30
	PIXEL_PER_MICRON = 6.2
	x = x.split('|') 
	loc = [float(a) for a in x[0].split(',')]
	radX,radY = [float(a) for a in x[1].split(',')]
	
	radX = (radX * X_PIXELS) / PIXEL_PER_MICRON
	radY = (radY * Y_PIXELS) / PIXEL_PER_MICRON
	rad = [radX,radY] 

	size = rad[0] * rad[1] * math.pi 
	return loc,rad,size 





class SizePlot:
        def __init__(self,options,xLen=3,yLen=4,key={}):

                self.options = options
                self.VERBOSE = True
                
        
                sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
                self.fig.set_facecolor('white') 
                self.fig.patch.set_facecolor('white')
                matplotlib.rcParams['savefig.facecolor'] = 'white'
                matplotlib.rcParams['ytick.labelsize'] = 7.5
                
                seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})

                #self.fig.patch.set_facecolor('lightgrey')
                self.fig.patch.set_facecolor('white')
                self.xLen, self.yLen = xLen,yLen 
                self.xLoc, self.yLoc = 0,0 
		self.color_key = {'NA': 'k', 'UNK': 'grey','+': 'cyan','-': 'olive','16': 'lightblue','17': 'lime','18': 'red','14':'yellow','19':'purple'} 

		self.color_key['SVZ'] = 'crimson'
		self.color_key['SP'] = 'firebrick'
		self.color_key['CP'] = 'gold'
		self.color_key['MZ'] = 'olive'

		
		self.color_offset = 0


	def finish(self,out='BASIC',title='figure',xLen=2,yLen=3,kind='RAW'):
		out_name = out+'.png'
		if kind == 'REL': 	
			plt.suptitle(title+' (RELATIVE SIZE)',fontsize=25,fontweight='bold') 
		else: plt.suptitle(title,fontsize=25,fontweight='bold') 
   		plt.subplots_adjust(left=0.05, bottom=0.075, right=0.95, top=0.92,wspace=0.10,hspace=0.45)
		
		plt.show() 
		plt.savefig(out_name,dpi=200) 
		plt.clf() 

                self.xLen, self.yLen = xLen,yLen 
		self.xLoc,self.yLoc = 0, 0 


	def get_color(self,x): 

		if len(x.split(',')) == 2: 
			a,b = x.split(',') 
			if a not in self.color_key: 
				self.color_key[a] = COLORS[self.color_offset]
				self.color_offset += 1
				
			if b not in self.color_key: 
				self.color_key[b] = COLORS[self.color_offset]
				self.color_offset += 1

			return self.color_key[a],self.color_key[b] 
	
		else:
			try: return self.color_key[x] 

			except KeyError: 
				self.color_key[x] = COLORS[self.color_offset]
				self.color_offset +=1
			return self.color_key[x] 








        def add_hist(self,h_data,TITLE=None,DNAME=None,MINSIZE=30):

		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		if type(h_data) == list:
			y = sorted(h_data) 
			if len(self.ax.hist(y,bins='auto')[0]) < 2:
				self.ax.clear() 
				self.ax.hist(np.hstack(y),bins=min(4,len(cc(y))))

		else:

			x,y = h_data.keys(), sorted(h_data.values())

			X_key = dd(lambda: dd(int))


			R,pv = 'NA','NA'


			yOffset = 0 
			xOffset = 0 
			xLoc    = 0 
			xticks = [] 
			if DNAME == 'REL': 
				for X,Y in h_data.items(): 
					if len(Y) < MINSIZE: continue 

					for y in Y: 
						yAdd = int(y) 
						if yAdd < -5: yAdd = -5 
						if yAdd > 5:  yAdd = 5 
						X_key[X][yAdd] += 1 

				

				Y_key = dd(lambda: dd(int))
				for xx in X_key: 
					xsum = sum(X_key[xx].values()) 
					for yy in X_key[xx]: Y_key[xx][yy] = X_key[xx][yy] / float(xsum) 		

		
				try:
					my_means = []  
					my_x = [int(xx) for xx in X_key.keys()]
					for xx in X_key.keys():
						x_len = sum(X_key[xx].values()) 
						x_total = sum([int(a)*b for a,b in X_key[xx].items()]) 
						my_means.append(x_total/x_len) 
					R,pv = stats.pearsonr(my_x, my_means)
				except ValueError: 
					my_ex = [] 
				
				
				X_key = Y_key

				xJump = len(X_key)
				hRange = [-3,-2,-1,0,1,2,3,4,5] 
				items,labels = [],[]
				color_idx = 0 
				for X in X_key:
					xLoc = 0
					labels.append(X) 
					clr = self.get_color(X)  
					if type(clr) == tuple:  
						a,b = clr 
						items.append(Rect((0,0),1,1,fc=a,ec=b,hatch='*'))
						for s in hRange:
							self.ax.bar(xLoc+xOffset,X_key[X][s],width=1,bottom=0,color=a,ec=b,hatch='*')

							if xOffset == 0: xticks.append((xLoc+(xJump/2),s))			

							xLoc += xJump

						
 
						xOffset +=1 
					else:
						items.append(Rect((0,0),1,1,fc=clr))
						for s in hRange:
							self.ax.bar(xLoc+xOffset,X_key[X][s],width=1,bottom=0,color=clr)
							if xOffset == 0: xticks.append((xLoc+(xJump/2),s))			
							xLoc += xJump 
						xOffset +=1 



			else:
				for X,Y in h_data.items(): 
					if len(Y) < MINSIZE: continue 
					for y in Y:
						yAdd = y%100 
						if yAdd >= 50: yAdd = 50
						else: 	       yAdd = 0 
						yRnd =(int(y/100) *100)+yAdd
						if yRnd > 1000: yRnd = 1000 
						X_key[X][yRnd] += 1 


				Y_key = dd(lambda: dd(int))
				for xx in X_key: 
					xsum = sum(X_key[xx].values()) 
					for yy in X_key[xx]: Y_key[xx][yy] = X_key[xx][yy] / float(xsum) 		

				try:
					my_means = []  
					my_x = [int(xx) for xx in X_key.keys()]
					for xx in X_key.keys():
						x_len = sum(X_key[xx].values()) 
						x_total = sum([int(a)*b for a,b in X_key[xx].items()]) 
						my_means.append(x_total/x_len) 
					R,pv = stats.pearsonr(my_x, my_means)
				except ValueError: 
					my_ex = [] 
				
				
				X_key = Y_key


				xJump = len(X_key)
				hRange = range(0,1000,50) 
				items,labels = [],[]
				color_idx = 0 
				for X in X_key:
					xLoc = 0
					labels.append(X) 
					clr = self.get_color(X)  
					if type(clr) == tuple:  
						a,b = clr 
						items.append(Rect((0,0),1,1,fc=a,ec=b,hatch='*'))
						for s in hRange:
							self.ax.bar(xLoc+xOffset,X_key[X][s],width=1,bottom=0,color=a,ec=b,hatch='*')
							if xOffset == 0: xticks.append((xLoc+(xJump/2),s))			
							xLoc += xJump 
						xOffset +=1 
					else:
						items.append(Rect((0,0),1,1,fc=clr))
						for s in hRange:
							self.ax.bar(xLoc+xOffset,X_key[X][s],width=1,bottom=0,color=clr)
							if xOffset == 0: xticks.append((xLoc+(xJump/2),s))			
							xLoc += xJump 
						xOffset +=1 



			xt,xl= [xx[0] for xx in xticks], [str(xx[1]) for xx in xticks]
			self.ax.set_xticks(xt)  
			self.ax.set_xticklabels(xl)  
			

			if R != 'NA': 
				rSS = str(round(R,5))
				pSS = str(round(pv,8))
				self.ax.set_title('CORR: '+rSS+' pv= '+pSS)
			if TITLE:
				nc,fs,cs,hp,bb1,bb2 = 1,12 ,0.5,0.5,0.9,1.1
				if len(items) > 5: 
					if len(items) < 10: 		    
						nc,fs,cs,hp = 2,10,0.4,0.4
					elif len(items) < 15: 		    
						nc,fs,cs,hp,bb1,bb2 = 3,10,0.35,0.35,0.85,1.1
					elif len(items) < 25:
						nc,fs,cs,hp,bb1,bb2 = 4,9,0.25,0.25,0.85,1.1
					elif len(items) < 35:
						nc,fs,cs,hp,bb1,bb2 = 5,8,0.20,0.20,0.8,1.1
					else:
						nc,fs,cs,hp,bb1,bb2 = 6,7,0.15,0.15,0.75,1.1
				leg = self.ax.legend(items,labels,title=TITLE,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')
			
		self.yLoc += 1 	
		if self.yLoc == self.yLen:
			self.xLoc +=1 
			self.yLoc = 0 








        def add_boxes(self,h_data,TITLE=None,DNAME=None,MINSIZE=30,BW=True):

		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		self.xOffset = 0 
		#clr = self.get_color(X)  
		iKey = dd(int) 
		iKey['16'], iKey['17'],iKey['18'],iKey['GT_18'] = 1,2,3,4
		iKey['IZ'], iKey['SP'],iKey['CP'],iKey['CR'],iKey['MZ'] = 1,2,3,4,5
		iKey['SINGLE'], iKey['BRIEF'],iKey['REPETITIVE'] = 1,2,3
		iKey['*SVZ*'], iKey['IZ'],iKey['*IZ*'] = 0.5,1,2


		data = sorted([[iKey[d[0]],d[0],d[1]] for d in h_data.items()])

		data = [opt for opt in data if data[1] not in ['NA','UNK','UNKNOWN','UNCLEAR','NON']]

		opts,vals = [],[]
			
 
		for d in data: 

			d[2] = [float(dx) for dx in d[2] if dx != 'NA']
			
			if d[1] in ['NA','UNK','UNKNOWN','UNAVAIL','UNCLEAR','DOUBLETS']: continue 
			elif d[1] == 'GT_18':  opts.append('>18') 
			elif d[1] == 'SUB_16': opts.append('<16') 
			elif len(d[1].split('_')) > 1: continue 
			else:
				opts.append(d[1]) 

			vals.append(d[2])  


		colors = [self.get_color(X) for X in opts]
		colors = ['pink','orange','purple','orange']	
		if BW: colors = ['k' for X in opts] 
		if opts == ['SVZ', '*SVZ*', 'IZ', '*IZ*']: opts = ['$SVZ$','$SVZ_{novel}$','$IZ$','$IZ_{novel}$']
 	


		items,labels = [],[]
		if type(colors[0]) == tuple:
			colors,alt_colors = [c[0] for c in colors], [c[1] for c in colors]
		else:
			alt_colors = [c for c in colors]
		R,pv = 'NA','NA'


                pos, width = [self.xOffset+0.20*x for x in range(len(opts))], 0.16
                self.bp = self.ax.boxplot(vals,positions=pos, widths=width, patch_artist=True, showmeans=True,whis=1)
                clevels = np.linspace(0., 3., len(opts)+1)


                xticks = []     
                means = sorted([(np.mean(v),opt) for (v,opt) in zip(vals,opts) if opt[0:2] == 'si'])
                medians = sorted([(np.percentile(v,60),opt) for (v,opt) in zip(vals,opts) if opt[0:2] == 'si'])
                sorts = [(sorted(v,reverse=True),opt) for (v,opt) in zip(vals,opts) if opt[0:2] == 'si']

                for i,opt in enumerate(opts):
                        clr,ac = colors[i], alt_colors[i] 

 
                        self.bp['boxes'][i].set_edgecolor(clr) 
                        self.bp['boxes'][i].set_linewidth(1) 
                        plt.setp(self.bp['medians'][i], color=clr, linewidth=3)
                        plt.setp(self.bp['means'][i], marker='h',markersize=9,markerfacecolor=clr) 
                        plt.setp(self.bp['caps'][(i*2)+1], color=clr,linewidth=1) 
                        plt.setp(self.bp['caps'][(i*2)], color=clr,linewidth=1) 
                        plt.setp(self.bp['whiskers'][i*2], color=clr,linewidth=1) 
                        plt.setp(self.bp['whiskers'][1+(i*2)], color=clr,linewidth=1)   
                        plt.setp(self.bp['fliers'][i], markerfacecolor=clr, markeredgecolor = clr, marker='s',markersize=2.0)           

			if clr == ac:	items.append(Rect((0,0),1,1,fc=clr))
			else:   	items.append(Rect((0,0),1,1,fc=a,ec=b,hatch='*'))
			labels.append(opt) 



                clevels = np.linspace(0., 1., len(vals))
                xJitters = [np.random.normal(pos[i],0.01,len(vals[i])) for i in range(len(vals))]
                yJitters = [[vals[i][j]*np.random.normal(1,0.005,len(vals[i]))[j] for j in range(len(vals[i]))] for i in range(len(vals))]

                for xJ, yJ, clevel in zip(xJitters, yJitters, colors):
                        plt.scatter(xJ, yJ, c=clevel, alpha=0.7,s=4,zorder=9)


                for patch, clevel in zip(self.bp['boxes'], colors):
                        patch.set_facecolor(clevel) #cm.prism(clevel))
                        patch.set_alpha(0.5)

		

		yMin,yMax = self.ax.get_ylim() 
		yStep = (yMax - yMin)/30.0

                
       #         self.ax.set_xticklabels([opt.split('~')[-1] for opt in opts],fontsize=10,fontweight='bold',rotation=-30) 
		for j,xpos in enumerate(self.ax.get_xticks()):
			#self.ax.text(xpos,yMin-yStep,opts[j].split('~')[-1],fontweight='bold',rotation=60,verticalalignment='top',horizontalalignment='center') 
			self.ax.text(xpos,yMin+yStep,opts[j].split('~')[-1],fontsize=25,fontweight='bold',verticalalignment='bottom',horizontalalignment='center') 


                self.ax.set_xlim([pos[0]-0.2,pos[-1]+0.2])
#		self.ax.plot((pos[0],pos[-1]),(yMin-yStep*2,yMin-2*yStep))

		self.ax.axis('off') 

		title_key = {'LOC': 'LAYER','FS': 'FIRING STYLE'} 

		if TITLE:

			if TITLE in title_key: TITLE = title_key[TITLE]
			nc,fs,cs,hp,bb1,bb2 = 1,12 ,0.5,0.5,0.2,1.1
			if BW: self.ax.set_title(TITLE) 
			else:
				if len(items) > 4: 
					if len(items) < 10: 		    
						nc,fs,cs,hp = 2,10,0.4,0.4
					elif len(items) < 15: 		    
						nc,fs,cs,hp,bb1,bb2 = 3,10,0.35,0.35,0.85,1.1
					elif len(items) < 25:
						nc,fs,cs,hp,bb1,bb2 = 4,9,0.25,0.25,0.85,1.1
					elif len(items) < 35:
						nc,fs,cs,hp,bb1,bb2 = 5,8,0.20,0.20,0.8,1.1
					else:
						nc,fs,cs,hp,bb1,bb2 = 6,7,0.15,0.15,0.75,1.1
				leg = self.ax.legend(items,labels,title=TITLE,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')
			



















		self.yLoc += 1 	
		if self.yLoc == self.yLen:
			self.xLoc +=1 
			self.yLoc = 0 






















def run_script(data_file,options):
	k=0
	key = dd(lambda: {}) 
	efiz = dd(lambda: {}) 
	cell_key = dd(lambda: {}) 
	nick_key = {'MEDIUM': 'MED', 'UNKNOWN': 'UNK','SIZE-UNCLEAR': 'UNCLEAR','NONBIG': 'NON','LARGE': 'BIG','IZ_SVZ': 'SVZ_IZ'} 
	for line in open(data_file): 
		line = line.split() 
		if line[0] == '---': 
 			headers = line
		else:
			cell = line[0] 
			for i in range(1,len(headers)):

				if headers[i] in ['GW','OBS','TOT','area','dia1','dia2','relFC','relFM']:

	
					try:
						v = float(line[i])
					except ValueError: 	
						v = line[i] 
				else:
					v = line[i] 
				if v == 'CR': v = 'MZ'
				key[headers[i]][cell] =  v
				cell_key[cell][headers[i]] = v 





	d_key = dd(lambda: dd(list)) 
	r_key = dd(lambda: dd(list)) 
	m_key = dd(lambda: dd(list)) 

	ed_key = dd(lambda: dd(list)) 
	er_key = dd(lambda: dd(list)) 
	em_key = dd(lambda: dd(list)) 


	split_keys = [('SURELOC',p) for p in ['MEGASIZE','ORIGSIZE']]
	split_keys += [('LOC',p) for p in ['MEGASIZE','ORIGSIZE']]
	split_keys += [('vGW',p) for p in ['MEGASIZE','ORIGSIZE','LOC','SURELOC']]



	keys = ['GW', 'BIO', 'FS', 'LOC', 'relFC', 'TOT', 'BGW', 'refFM', 'MEGATYPE', 'dia1', 'dia2', 'OBS', 'size']
	keys = ['GW', 'BIO', 'FS', 'LOC','TOT', 'BGW', 'MEGATYPE', 'OBS' ]
	cont_keys = ['GW','TOT', 'BGW', 'OBS' ]
	bin_keys = ['FS', 'LOC', 'BGW', 'MEGATYPE' ]
	
	within_keys = ['IZ_SVZ_SIZE','SP_CP_GW'] 


	cont_groups         = dd(list) 
	cont_groups_rel         = dd(list) 
	cont_groups_no_mega = dd(list) 
	cont_groups_no_mega_rel = dd(list) 
	cont_groups_cp_sp = dd(list) 
	cont_groups_cp_sp_rel = dd(list) 
	cont_groups_iz_svz = dd(list) 
	cont_groups_iz_svz_rel = dd(list) 

	


	bin_groups_special = dd(lambda: dd(list)) 
	bin_groups = dd(lambda: dd(list)) 
	bin_groups_rel = dd(lambda: dd(list)) 
	bin_groups_no_mega = dd(lambda: dd(list)) 
	bin_groups_no_mega_rel = dd(lambda: dd(list)) 
	bin_groups_cp_sp = dd(lambda: dd(list)) 
	bin_groups_cp_sp_rel = dd(lambda: dd(list)) 
	bin_groups_iz_svz = dd(lambda: dd(list)) 
	bin_groups_iz_svz_rel = dd(lambda: dd(list)) 
	for c in cell_key:

		ck = cell_key[c]	


		sz,d1,d2,relFC,relFM = ck['area'],ck['dia1'],ck['dia2'],ck['relFC'],ck['relFM']

		for k in bin_keys: 
			bin_groups[k][ck[k]].append(sz) 
			bin_groups_rel[k][ck[k]].append(relFC) 


		for k in cont_keys: 
			cont_groups[k].append([ck[k],sz])
			cont_groups_rel[k].append([ck[k],relFC])


		cBGW, cLoc, cMEGATYPE = ck['BGW'],ck['LOC'], ck['MEGATYPE']
		if cMEGATYPE != 'MEGA': 

			for k in bin_keys: 
				bin_groups_no_mega[k][ck[k]].append(sz) 
				bin_groups_no_mega_rel[k][ck[k]].append(relFC) 

			for k in cont_keys: 
				cont_groups_no_mega[k].append([ck[k],sz])
				cont_groups_no_mega_rel[k].append([ck[k],relFC])



		if cLoc in ['SP','CP','SP_CP']: 
			for k in bin_keys: 
				bin_groups_cp_sp[k][ck[k]].append(sz) 
				bin_groups_cp_sp_rel[k][ck[k]].append(relFC) 

			for k in cont_keys: 
				cont_groups_cp_sp[k].append([ck[k],sz])
				cont_groups_cp_sp_rel[k].append([ck[k],relFC])


		if cLoc in ['IZ','SVZ','IZ_SVZ','SVZ_IZ']:

			 
			for k in bin_keys: 
				bin_groups_iz_svz[k][ck[k]].append(sz) 
				bin_groups_iz_svz_rel[k][ck[k]].append(relFC) 

			for k in cont_keys: 
				cont_groups_iz_svz[k].append([ck[k],sz])
				cont_groups_iz_svz_rel[k].append([ck[k],relFC])

		if cLoc in ['IZ','SVZ'] and cMEGATYPE in ['MINI','MEGA']: 
			if cMEGATYPE == 'MEGA': 
				k = '*'+cLoc+'*'
			else:
				k = cLoc 
			bin_groups_special['GSIZE'][k].append(sz) 


		continue  
		if  len(cell_key[c]['near']) > 0:
			N = [x[-1] for x in cell_key[c]['near'] if x[-1] > 0.0005]
			if len(N)  > 0: 

				
				Nmed, Nmean = np.median(N), np.mean(N) 
				if cz > Nmed: FM = cz / Nmed
				else: 	      FM = -1 * (Nmed/cz) 
				if cz > Nmean: FC = cz / Nmean
				else: 	      FC = -1 * (Nmean/cz) 
				

				for k in key.keys(): 
					if key[k][c] != 'NA': 
						r_key[k][key[k][c]].append(FC) 
						m_key[k][key[k][c]].append(FM) 
						if ms != 'NA': em_key[k][key[k][c]].append(ms)
						if ms != 'NA': er_key[k][key[k][c]].append(ms)

				for x,y in split_keys:
					if key[x][c] != 'NA' and key[y][c] != 'NA': 
						r_key[(x,y)][(key[x][c],key[y][c])].append(FC) 
						m_key[(x,y)][(key[x][c],key[y][c])].append(FM) 

						if ms != 'NA':  er_key[(x,y)][(key[x][c],key[y][c])].append(ms) 
						if ms != 'NA':  em_key[(x,y)][(key[x][c],key[y][c])].append(ms) 
			cell_key[c]['relsize'] = (FC,FM) 


	

	
	sp = SizePlot(options,xLen=1,yLen=2) 
#	sp.add_hist(sizes['ALL'])
#	sp.add_hist(sizes['REL'])
#	sp.finish('ALL_SIZES',title='ALL (RAW/REL) CELL SIZES',xLen=3,yLen=2) 





	basic_types = ['LOC', 'SURELOC','MEGASIZE','ORIGSIZE', 'vGW','vBIO']  
	split_types = ['LOC_&_ORIGSIZE', 'SURELOC_&_ORIGSIZE', 'vGW_&_ORIGSIZE', 'LOC_&_MEGASIZE', 'SURELOC_&_MEGASIZE', 'vGW_&_MEGASIZE']
	cr_types = ['MZ_ID','MZ_GW','CR_GW'] 
	gw_types = ['CP_GW-','SP_GW-','SVZ_IZ_GW-']
	gw_pos_types = ['GW+','SVZ_GW+','IZ_GW+']



	sp = SizePlot(options,xLen=1,yLen=1) 
	for b in bin_groups_special.keys(): sp.add_boxes(bin_groups_special[b],b,'RAW')
	sp.finish('BASIC_SIZE_TYPES',title='SIZES BY TYPE',xLen=1,yLen=1 ) 


	sys.exit() 	



	sp = SizePlot(options,xLen=2,yLen=2) 
	for b in bin_keys: sp.add_boxes(bin_groups[b],b,'RAW')
	sp.finish('BASIC_SIZE_TYPES',title='SIZES BY TYPE',xLen=2,yLen=2 ) 





	for b in bin_groups_no_mega: sp.add_boxes(bin_groups_no_mega[b],b,'RAW')
	sp.finish('BASIC_SIZES_WITHOUT_WOLSELEY_CELLS',title='SIZES BY TYPE WITHOUT WOLSELEY CELLS',xLen=2,yLen=2 ) 


	for b in bin_groups_cp_sp: sp.add_boxes(bin_groups_cp_sp[b],b,'RAW')
	sp.finish('BASIC_WITHIN_CP_SP',title='SIZES WITHIN CP/SP',xLen=2,yLen=2 ) 


	for b in bin_groups_iz_svz: sp.add_boxes(bin_groups_iz_svz[b],b,'RAW')
	sp.finish('BASIC_WITHIN_IZ_SVZ',title='SIZES WITHIN IZ/SVZ',xLen=2,yLen=2 ) 


	for b in bin_keys: sp.add_boxes(bin_groups_rel[b],b,'REL')
	sp.finish('BASIC_SIZE_TYPES',title='SIZES BY TYPE',xLen=2,yLen=2,kind='REL' ) 

	for b in bin_groups_no_mega: sp.add_boxes(bin_groups_no_mega_rel[b],b,'REL')
	sp.finish('BASIC_SIZES_WITHOUT_WOLSELEY_CELLS',title='SIZES BY TYPE WITHOUT WOLSELEY CELLS',xLen=2,yLen=2,kind='REL' ) 


	for b in bin_groups_cp_sp: sp.add_boxes(bin_groups_cp_sp_rel[b],b,'REL')
	sp.finish('BASIC_WITHIN_CP_SP',title='SIZES WITHIN CP/SP',xLen=2,yLen=2,kind='REL' ) 


	for b in bin_groups_iz_svz: sp.add_boxes(bin_groups_iz_svz_rel[b],b,'REL')
	sp.finish('BASIC_WITHIN_IZ_SVZ',title='SIZES WITHIN IZ/SVZ',xLen=2,yLen=2,kind='REL' ) 








	sys.exit() 






	sys.exit() 
	for b in split_types:	sp.add_boxes(size_types[b],b,'RAW')
	sp.finish('MULTI_SIZE_TYPES',title='SIZES BY MULTITYPE',xLen=3,yLen=1) 



	for b in cr_types: sp.add_boxes(size_types[b],b,'RAW',MINSIZE=5)
	sp.finish('MZ_SPECIFIC_SIZE_TYPES',title='SPECIFC SIZE BY MZ',xLen=3,yLen=1 ) 
	for b in gw_types: sp.add_boxes(size_types[b],b,'RAW',MINSIZE=2)
	sp.finish('SPECIFIC_SMALL_SIZE_LOCS',title='SPECIFC SIZE BY SMALL LOCS',xLen=3,yLen=1 ) 
	for b in gw_pos_types: sp.add_boxes(size_types[b],b,'RAW',MINSIZE=2)
	sp.finish('SPECIFIC_MEGA_SIZE_LOCS',title='SPECIFC SIZE BY MEGA LOCS',xLen=3,yLen=2 ) 

	sys.exit() 
	for b in basic_types: 	sp.add_hist(rel_types[b],b,'REL')
	sp.finish('REL_BASIC_SIZE_TYPES',title='REL SIZES BY TYPE',xLen=2,yLen=3 ) 
	for b in split_types:	sp.add_hist(rel_types[b],b,'REL')
	sp.finish('REL_MULTI_SIZE_TYPES',title='REL SIZES BY MULTITYPE',xLen=3,yLen=1) 

	for b in cr_types: sp.add_hist(rel_types[b],b,'REL',MINSIZE=5)
	sp.finish('REL_MZ_SPECIFIC_SIZE_TYPES',title='REL SPECIFC SIZE BY MZ',xLen=3,yLen=1 ) 
	for b in gw_types: sp.add_hist(rel_types[b],b,'REL',MINSIZE=2)
	sp.finish('REL_SPECIFIC_SMALL_SIZE_LOCS',title='REL SPECIFC SIZE BY SMALL LOCS',xLen=3,yLen=1 ) 
	for b in gw_pos_types: sp.add_hist(rel_types[b],b,'REL',MINSIZE=2)
	sp.finish('REL_SPECIFIC_MEGA_SIZE_LOCS',title='REL SPECIFC SIZE BY MEGA LOCS',xLen=3,yLen=2 ) 

	sys.exit() 






	w=sys.stdout


	ws = '%-30s %12s %20s %20s %12s %30s %30s\n'
	wx = '%-30s %12s %20s %20s %12s %30s %30s %20s %20s %20s\n'
	for k in d_key.keys():
#		continue 
		for x in d_key[k].keys(): 
			xd = d_key[k][x]
			


			efz = ed_key[k][x]
			eMed,eMean = np.median(efz), np.mean(efz) 
			if np.isnan(eMed):  eMed = 0 
			if np.isnan(eMean): eMean = 0 
			

			if type(k) == tuple: 
				kj = '&'.join(k) 
				xj = ','.join(x) 
				w.write(wx % (kj,'RAW-SIZE', xj ,len(xd),'MEAN/MED',np.mean(xd),np.median(xd),EF,eMed,eMean))
			else:
				w.write(wx % (k, 'RAW-SIZE',x ,len(xd),'MEAN/MED',np.mean(xd),np.median(xd),EF,eMed,eMean))



	for k in r_key.keys():
		for x in r_key[k].keys(): 
			xd = r_key[k][x]
			zd = m_key[k][x]
			efz = er_key[k][x]
			eMed,eMean = np.median(efz), np.mean(efz)
			if np.isnan(eMed):  eMed = 0 
			if np.isnan(eMean): eMean = 0 
			if type(k) == tuple: 
				kj = '&'.join(k) 
				xj = ','.join(x) 
				w.write(wx % ( kj,'MEAN-REL-SIZE', xj,len(xd),'MEAN/MED',np.mean(xd),np.median(xd),EF,eMed,eMean))
				w.write(wx % ( kj, 'MED-REL-SIZE',  xj,  len(zd),'MEAN/MED',np.mean(zd),np.median(zd),EF,eMed,eMean ))
			else:
				w.write(wx  % (k,'MEAN-REL-SIZE',x ,len(xd),'MEAN/MED',np.mean(xd),np.median(xd),EF,eMed,eMean))
				w.write(wx  % (k,'MED-REL-SIZE',x ,len(zd),'MEAN/MED',np.mean(zd),np.median(zd),EF,eMed,eMean ))















if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args[0],options)	














