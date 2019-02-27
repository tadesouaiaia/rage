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
                
		self.xOffset = 0 
        
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
		self.color_key['ADULT'] = 'red'
		self.color_key['ADULT'] = 'red'
		self.color_key['ES'] = 'cyan'
		self.color_key['CR'] = 'green'
		self.color_key['MINI'] = 'gray'

		self.iKey = dd(int) 

		self.iKey['MEGA'],self.iKey['ES'],self.iKey['CR'],self.iKey['ES'],self.iKey['ADULT'] = 1,2,3,4,5
		
		self.color_offset = 0


	def finish(self,out='BASIC',title='figure',xLen=2,yLen=3,kind='RAW'):
		out_name = out+'.png'
		plt.suptitle(title.upper(),fontsize=40,fontweight='bold') 
   		plt.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.90,wspace=0.10,hspace=0.05)
		
		plt.savefig(out_name,dpi=200) 
		plt.show() 
		plt.clf() 

                self.xLen, self.yLen = xLen,yLen 
		self.xLoc,self.yLoc = 0, 0 
		self.xOffset = 0 


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








        def add_boxes(self,h_data,TITLE=None,DNAME=None,MINSIZE=30,BW=True,MULTI=False):

		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

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
		

		if not MULTI: 
			self.yLoc += 1 	
			if self.yLoc == self.yLen:
				self.xLoc +=1 
				self.yLoc = 0 












	def add_cartoon(self):


		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		self.ax.set_xlim(0,100)
		self.ax.set_ylim(0,100)

		self.ax.plot((0,40),(40,40),'-',linewidth=5,color='k')
		self.ax.plot((40,60),(40,40),'--',linewidth=5,color='k')
		self.ax.plot((60,100),(40,40),'-',linewidth=5,color='k')
		self.ax.text(50,50,'PRE-SYNAPSE',horizontalalignment='center',verticalalignment='bottom',fontweight='bold',fontsize=30)
		self.ax.text(50,30,'POST-SYNAPSE',horizontalalignment='center',verticalalignment='top',fontweight='bold',fontsize=30)

		self.ax.axis('off') 
		self.xLoc+=1



        def add_multi_boxes(self,h_data,TITLE=None,HEADER=None,MINSIZE=30,BW=True,LEG=False,STEP=30):



		self.items, self.xOffset, self.yTop ={}, 0 , 0 
		box_names =  h_data.keys()

		self.axes = [plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)]
		box_names = [b[1] for b in sorted([(b.split(';')[-1],b) for b in box_names]) if b[0] not in ['GABBR1']][0:10]
		self.scale_key = dd(list) 
			
		for i in range(0,len(box_names),STEP):
			
			for x,b in enumerate(box_names[i:i+STEP]): 
				self.add_multi_box(h_data[b],b,SCALE=5)

			if i+STEP >= len(box_names):
				self.add_multi_box({opt: [sum([vals[x][y] for x in range(len(vals))]) for y in range(len(vals[0]))] for opt,vals in self.scale_key.items()},'SCALED_TOTAL')
				
			self.axes[-1].set_xlim(0-0.2,self.xOffset)
			self.axes[-1].set_ylim(0,self.yTop) 
			self.axes[-1].set_xticks([])
			self.axes[-1].axis('off') 
			self.xLoc += 1 
			self.xOffset = 0 
			try: self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))
			except: break 



		nc,fs,cs,hp,bb1,bb2 = len(self.items),20 ,1,1,0.5,-.20
		if LEG:
			leg_data = sorted([[self.iKey[x[0]],x[0],x[1]] for x in self.items.items()])
			labels,items = [ld[1] for ld in leg_data],[ld[2] for ld in leg_data]
			
			#data = [d for d in sorted([[self.iKey[d[0]],d[0],d[1]] for d in h_data.items()]) if d[1] not in ['MA','UN']]

			leg = self.axes[-1].legend(items,labels,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')
		self.axes[-1].axis('off') 
#		self.xLoc +=1












	def add_multi_box(self,h_data,TITLE,BW=False,SCALE=0):

		data = [d for d in sorted([[self.iKey[d[0]],d[0],d[1]] for d in h_data.items()]) if d[1] not in ['MA','UN']]
		opts,vals,items,labels,R,pv = [d[1] for d in data],[[v for v in  d[2]] for d in data],[],[],'NA','NA'
		colors, sigs,jitters = [self.get_color(X) for X in opts], [False for X in opts],[] 
		if type(colors[0]) == tuple:	colors,alt_colors = [c[0] for c in colors], [c[1] for c in colors]
		else:				alt_colors = [c for c in colors]
		
                self.pos, self.width = [self.xOffset+0.20*x for x in range(len(opts))], 0.16
		all_vals = [a for b in vals for a in b] 
		valMin,valMean,valMax = min(all_vals), np.mean(all_vals), max(all_vals) 

		for i in range(len(vals)): 
			if np.mean(vals[i]) > valMean: 
				other_vals = [a for b in [vals[j] for j in range(len(vals)) if j != i] for a in b] 
				if stats.ttest_ind(vals[i],other_vals)[-1]  < 0.05: sigs[i] = True
                	jitters.append([[self.pos[i] + np.random.normal(0,0.02) for j in range(len(vals[i]))],[v*np.random.normal(1,0.005) for v in vals[i]]])

			if SCALE > 0: 
				scaler = MinMaxScaler(feature_range=(0,SCALE))
				self.scale_key[opts[i]].append([vs for vs in scaler.fit_transform(np.array(vals[i],dtype=float).reshape(-1,1)).reshape(1,-1)[0]])
			


		self.fill_boxes(opts,vals,colors,alt_colors,JITTERS=jitters,SIG_STARS=sigs)

		if valMax > self.yTop: self.yTop = valMax
		self.xOffset = self.pos[-1]+2*self.width                 
		if TITLE:	self.axes[-1].text(np.mean(self.pos),0-1,TITLE.split(';')[-1],fontsize=16,fontweight='bold',horizontalalignment='center',verticalalignment='top') 
		
	
		return 



	def fill_boxes(self,opts,vals,colors,alt_colors,JITTERS=None,SIG_STARS=None):


                self.bp = self.axes[-1].boxplot(vals,positions=self.pos, widths=self.width, patch_artist=True, showmeans=True,whis=0.7)
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
			if clr == ac:	self.items[opt] = Rect((0,0),1,1,fc=clr)
			else:   	self.items[opt] = Rect((0,0),1,1,fc=a,ec=b,hatch='*')


			

                for patch, clevel in zip(self.bp['boxes'], colors):
                        patch.set_facecolor(clevel); patch.set_alpha(0.5) 

		if JITTERS != None:
			for i,(xJ,yJ) in enumerate(JITTERS):
                        	plt.scatter(xJ, yJ, c=colors[i], alpha=0.7,s=4,zorder=9)

		if SIG_STARS != None:
			for i,SIG in enumerate(SIG_STARS):
				if SIG: plt.scatter(self.pos[i],max(vals[i])*1.2,s=100,marker='*',color='gold')



def quick_order(gene,G,key,skey):


	if 'NA' in key: key.pop('NA') 

	obs  = {Y: len([x for x in X if x > 0])/float(len(X)) for Y,X  in key.items()}
	avg  = {Y: np.mean(X) for Y,X  in key.items()}
	
	data = sorted([[k,avg[k],obs[k],key[k]] for k in key.keys()],key=lambda P: P[1],reverse=True)
	all_data = []



	for d in data:
		all_data.extend([(d[-1][j],d[0]) for j in range(len(d[-1]))])
	all_data.sort(reverse=True) 
	d_obs = len([x for x in all_data if x>0])/float(len(x))  	
	dc = cc([x[1] for x in all_data])		
	e100, e250 = {},{} 

	c100 =  cc([x[1] for x in all_data[0:100]])
	c250 =  cc([x[1] for x in all_data[0:250]])
	c400 =  cc([x[1] for x in all_data[0:400]])

	exp100,obs100 = [],[]
	exp250,obs250 = [],[]
	exp400,obs400 = [],[]

	mults = [] 
	for k,x in dc.items():
		kr = x/len(all_data) 
		e100[k] = int(kr*100) 
		e250[k] = int(kr*250) 
		o100,o250,o400 = 0,0 ,0
		if k in c100.keys():o100 = c100[k]  
		if k in c250.keys():o250 = c250[k]  
		if k in c400.keys():o400 = c400[k]  
		
		exp100.append(kr*100)
		exp250.append(kr*250) 
		exp400.append(kr*400) 
		obs100.append(o100) 
		obs250.append(o250) 
		obs400.append(o400)

		m100 = o100/(kr*100) 
		m250 = o250/(kr*250) 
		m400 = o400/(kr*400) 

		mults.append((k,m100,m250,m400,m100+m250+m400))


	p100 = stats.chisquare(obs100,f_exp=exp100)[-1]
	p250 = stats.chisquare(obs250,f_exp=exp250)[-1]
	p400 = stats.chisquare(obs400,f_exp=exp400)[-1]
	mults.sort(reverse=True,key=lambda X: X[-1])
	
	SD =  sorted([a for b in skey.values() for a in b],reverse=True)
	mHi,mLo = mults[0][0],mults[-1][0] 

	STOP = 100 
	loOuts = {} 
	hiOuts = {} 

	if mHi == 'UN': return 

	for i in range(len(SD)): 
		v,g,s = SD[i] 
		if i < STOP and v > 0 and g in [mLo,'UN']:
			hiOuts[s] = STOP - i 
	if 'NA' in key: key.pop('NA') 
			print gene,s,g,G,mHi,mLo,"TOO-HIGH",v,STOP-i

		elif g == mHi:			
			if i+STOP >= len(SD) or v == 0: 
				hiOuts[s] = 1		
				print gene, s,g,G,mHi,mLo,"TOO-LOW",v,1

	




def quick_dex(gene,G,key):

	if 'UN' in key: key.pop('UN') 
	if 'NA' in key: key.pop('NA') 
	obs  = {Y: len([x for x in X if x > 0])/float(len(X)) for Y,X  in key.items()}
	avg  = {Y: np.mean(X) for Y,X  in key.items()}

	data = sorted([[k,avg[k],obs[k],key[k]] for k in key.keys()],key=lambda P: P[1],reverse=True)
	data_obs = sorted([[k,avg[k],obs[k],key[k]] for k in key.keys()],key=lambda P: P[2],reverse=True)
	maxO,maxA = data_obs[0][0],data[0][0]

	max_obs = data_obs[0][2] 
	max_avg = data[0][1] 
	all_data = [a for b in [data[m][-1] for m in range(len(data))] for a in b] 

	all_obs = len([x for x in all_data if x>0])/float(len(all_data))



	global_key = {} 


	if G == 'BIOS': 
		obs_list= [d[0].split('_')[-1] for d in data_obs]
		avg_list= [d[0].split('_')[-1] for d in data]
		MEGA_RATE = len([m for m in obs_list if m == 'MEGA'])/len(obs_list)  

		o_obs = len([m for m in obs_list[0:18] if m == 'MEGA'])
		a_obs = len([m for m in avg_list[0:18] if m == 'MEGA'])
		o_miss = 18 - o_obs 
		a_miss = 18 - a_obs  			
		fexp=[MEGA_RATE*18,(1-MEGA_RATE)*18]
		CHI_OBS = stats.chisquare([o_obs,o_miss],f_exp=[MEGA_RATE*18,(1-MEGA_RATE)*18])[1]
		CHI_AVG = stats.chisquare([a_obs,a_miss],f_exp=[MEGA_RATE*18,(1-MEGA_RATE)*18])[1]

		oCnt = o_obs - fexp[0] 
		aCnt =  a_obs - fexp[0] 
		fail_obs = dd(int) 
		fail_avg = dd(int) 
		print gene,'CHI','MEGA_BIO','CNTS',round(oCnt,2),round(aCnt,2),'pvs',CHI_OBS,CHI_AVG
		if CHI_OBS <0.05 and CHI_AVG < 0.05: 

			if o_obs > fexp[0] and a_obs > fexp[0]: 

				obs_fail = [(18-z,d[0]) for z,d in enumerate(data_obs[0:18]) if d[0].split('_')[-1] != 'MEGA']
				avg_fail = [(18-z,d[0]) for z,d in enumerate(data[0:18]) if d[0].split('_')[-1] != 'MEGA']

				fail_obs = {x[1]: x[0] for x in obs_fail} 
				fail_avg = {x[1]: x[0] for x in avg_fail} 
				return True,'HI',fail_obs, fail_avg				


			elif o_obs < fexp[0] and a_obs < fexp[0]: 

				obs_fail = [(z,d[0]) for z,d in enumerate(data_obs[18::]) if d[0].split('_')[-1] != 'MEGA']
				avg_fail = [(z,d[0]) for z,d in enumerate(data[18::]) if d[0].split('_')[-1] != 'MEGA']
				fail_obs = {x[1]: x[0] for x in obs_fail} 
				fail_avg = {x[1]: x[0] for x in avg_fail} 
				return True,'LO',fail_obs, fail_avg				


			 

		return False,'NO','NA','NA' 
	print gene,G,'SUMMARY','OBS/CV',round(all_obs,3),round(stats.variation(all_data),4),'max_obs/avg',round(max_obs,3),round(max_avg,3),'maxIDS',maxO,maxA
	for i in range(len(data)): 
		
		iN,iA,iO,iC = data[i]
		others = [a for b in [data[m][-1] for m in range(len(data)) if m != i] for a in b] 
		GFC = iA /  (np.mean(others)+0.01)
		GP = stats.ttest_ind(others,iC)[1] 
		global_key[iN] = GFC,GP 			

		print gene,G,'GLOBAL',iN,'OBS/AVG',iO,iA,"FC/PV",GFC,GP 


		for j in range(i+1,len(data)):
			jN,jA,jO,jC = data[j] 
			PV = stats.ttest_ind(iC,jC)[1] 
			FC = iA/(jA+0.0001) 	
			if iA > jA: 
				print gene,G,"LOCAL",iN,jN,'FC/PV',round(FC,3),PV 
			else: 
				print gene,G,"LOCAL",jN,iN,'FC/PV',round(jA/(iA+0.001),2),PV

	return False,'NO','NA','NA' 



def run_script(data_file,options):
	k=0
	TOPICS = ['CT','DCT','SCT','SDCT'] 

	key = dd(lambda: {}) 

	for line in open(options.key):
		line = line.split() 
		if line[0] == '---': 
 			headers = line  
		else:
			s = line[0] 
			for i in range(1,len(line)): 
				key[headers[i]][s] = line[i]  
				
	anno_key = dd(list) 
	genes,cnts,annos = [],[],[]
	loc_key = dd(lambda: dd(list))
	nkey = dd(lambda: dd(lambda: {}))
	
	topic_key = {} 
	topic_key = dd(lambda: dd(lambda: dd(lambda: {})))

	chi_hi = dd(list) 
	chi_lo = dd(list) 

	for line in open(data_file): 
		line = line.split() 
		if line[0] == '---': 
 			samples = line[1::] 
			#for k in key.keys():
			for k in TOPICS:
				for i,s in enumerate(samples): 
					loc_key[k][key[k][s]].append(i) 	


		else:
			
			nkey = dd(lambda: dd(lambda: {}))
			
			gene,cnts = line[0],[log(float(x)+1,2) for x in line[1::]]
			for k in TOPICS:
				group_key = {} 
				sample_key = {} 
				for g in loc_key[k].keys(): 
					group_key[g] = [cnts[x] for x in loc_key[k][g]] 
					sample_key[g] = [(cnts[x],g,samples[x]) for x in loc_key[k][g]]
					
				quick_order(gene,k,group_key,sample_key) 
				continue 	
				chiBool,loc,obs,avg = quick_dex(gene,k,group_key) 

				if chiBool: 
					if loc == 'HI':
						for oo,xx in obs.items(): chi_hi[oo].append(xx) 	
					else:
						for oo,xx in obs.items(): chi_lo[oo].append(xx) 	

	return 
	for x in chi_hi.keys()+chi_lo.keys(): 

		xHi,xLo = chi_hi[x],chi_lo[x] 
		print x,'BIOFAILS','hi',sum(xHi),len(xHi),np.mean(xHi),'lo',sum(xLo),len(xLo),np.mean(xLo)





	sys.exit() 
				















if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args[0],options)	













