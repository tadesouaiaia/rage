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





class IsoPlot:
        def __init__(self,options,prefix,xLen=2,yLen=1,key={}):

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

		self.prefix = prefix
		self.color_offset = 0

	def setLens(self,xLen,yLen):
                self.xLen, self.yLen = xLen,yLen 
		self.xLoc,self.yLoc = 0, 0 
		self.xOffset = 0 


	def finish(self,gene,pt,out='BASIC',title='figure',xLen=2,yLen=3,kind='RAW'):
	
		k_key = {'CM': 'FETAL CELLTYPE', 'CT': 'CELL TYPE', 'FS': 'FIRING STYLE'} 
	
		out_name = 'isodex_'+self.prefix+'_'+gene.split(';')[-1]+'_'+pt+'.png'

		title = gene.split(';')[-1]+' '+k_key[pt] 


		plt.suptitle(title,fontsize=25,fontweight='bold') 
   		plt.subplots_adjust(left=0.05, bottom=0.075, right=0.95, top=0.90,wspace=0.15,hspace=0.50)
		
#		plt.show() 
		plt.savefig(out_name,dpi=200) 
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















	def add_legend(self):
		items =  self.items.values() 
		labels = self.items.keys() 
		nc,fs,cs,hp,bb1,bb2 = 1,14 ,0.5,0.5,0.9,0.8
		if len(items) > 4: 
			if len(items) < 10: 		    
				nc,fs,cs,hp = 2,12,0.4,0.4
			elif len(items) < 15: 		    
				nc,fs,cs,hp,bb1,bb2 = 3,10,0.35,0.35,0.85,1.1
			else:
				nc,fs,cs,hp,bb1,bb2 = 6,7,0.15,0.15,0.75,1.1
		leg = self.ax.legend(items,labels,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')




	def add_transcript_boxes(self,gene,exonic_data,valid_key,TITLE=None,BW=False,LEG=0):

	

		self.items = {} 

		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		self.add_iso_boxes(exonic_data,'EXONIC')
		self.ax.set_xlim(0-0.2,self.xOffset)
		self.ax.axis('off') 
		self.add_legend()



		self.xLoc,self.yLoc = 1, 0 
		self.xOffset = 0 

		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		k_idx,offsets = 0, []  
		for k,k_data in valid_key.items():
			self.add_iso_boxes(k_data,k) 
			k_idx += 1 
			if k_idx % 5 == 0 and k_idx != len(valid_key): 
				self.ax.set_xlim(0-0.2,self.xOffset)
				self.ax.axis('off') 
				self.xLoc += 1
				self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
				self.xOffset = 0 

		self.ax.set_xlim(0-0.2,self.xOffset)
		self.ax.axis('off') 
		self.xLoc +=1












	def add_iso_boxes(self,box_data,TITLE=None,BW=False):




		
		opts,counts,vals = [], [],[] 
		for obs,avg,name,cnts in box_data:


			opts.append(name) 
			counts.append(cnts)
			vals.append([log(c+1,2) for c in cnts])


		colors = [self.get_color(X) for X in opts]
#		colors = ['pink','orange','purple','orange']	
		if BW: colors = ['k' for X in opts] 
		#if opts == ['SVZ', '*SVZ*', 'IZ', '*IZ*']: opts = ['$SVZ$','$SVZ_{novel}$','$IZ$','$IZ_{novel}$']



		items,labels = [],[]
		if type(colors[0]) == tuple:
			colors,alt_colors = [c[0] for c in colors], [c[1] for c in colors]
		else:
			alt_colors = [c for c in colors]
		R,pv = 'NA','NA'


                pos, width = [self.xOffset+0.20*x for x in range(len(opts))], 0.16
                self.bp = self.ax.boxplot(vals,positions=pos, widths=width, patch_artist=True, showmeans=True,whis=1)
                clevels = np.linspace(0., 3., len(opts)+1)



		self.xOffset = pos[-1]+2*width 
		xMid = np.mean(pos) 
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

			if clr == ac:	self.items[opt] = Rect((0,0),1,1,fc=clr)

			else:   	self.items[opt] = Rect((0,0),1,1,fc=a,ec=b,hatch='*')



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

                
#		for j,xpos in enumerate(self.ax.get_xticks()):
#			self.ax.text(xpos,yMin-yStep,opts[j].split('~')[-1],fontweight='bold',rotation=60,verticalalignment='center',horizontalalignment='center') 
			#self.ax.text(xpos,yMin-yStep,opts[j].split('~')[-1],fontsize=25,fontweight='bold',verticalalignment='bottom',horizontalalignment='center') 
#		self.ax.set_xticks([]) 




		if TITLE:
			self.ax.text(xMid,yMax,TITLE.split(';')[-1],fontsize=20,fontweight='bold') 
		

		
	
		return 









def quickstat(list1,list2,maxPV):

	obsI,avgI,nI,cI = list1
	obsJ,avgJ,nJ,cJ = list2 

	double = sorted([(nI,cI,avgI,obsI),(nJ,cJ,avgJ,obsJ)],key=lambda X:X[0])
	
	nX,cX,aX,oX = double[0] 
	nY,cY,aY,oY = double[1] 

	pv = stats.ttest_ind(cX,cY)[1]

	if oX+oY < 0.5: return 'NA','NA','NA'

	if pv < maxPV: 
		if aY > aX: 
			return (nX,nY),1,pv 
		else:
			return (nX,nY),0,pv 
	return 'NA','NA','NA'






def plot_iso_data(gene,tran_key,options,data_file,modT=0.15,minT=0.05,maxPV = 0.01, maxTrans=10):

	prefix = data_file.split('/')[-1].split('.')[0]

	TEST = False 
	isoplot = IsoPlot(options,prefix) 
	INVALID_OPTS = ['SG','NA','OB']
	tK = [k for k in tran_key if k != 'EXONIC']
	for k in ['CM','CT','FS']: 
		opts = tran_key['EXONIC'][k].keys() 
		valid_opts = []
		top_key   = {} 
		valid_key = {} 
		pv_lists = []
		tup_key = dd(lambda: [[],[]])
		for t in ['EXONIC']+tK:
			k_data = [] 
			for opt in opts: 
				if t == 'EXONIC' or opt in valid_opts: 
					cnts = tran_key[t][k][opt]
					obs = len([c for c in cnts if c>0])/len(cnts) 
					avg = np.mean(cnts) 
					k_data.append([obs,avg,opt,cnts]) 
			
			
			obs_k = sorted(k_data,reverse=True)
			avg_k = sorted(k_data,reverse=True,key=lambda X: X[1])			
			obsV,obsN = obs_k[0][0],obs_k[0][2]
			avgV,avgN = avg_k[0][1],avg_k[0][2]

			
			
			if t == 'EXONIC': 
				if obsV < modT or obs_k[1][0] < minT: return 
				else:   
					valid_opts = [kd[2] for kd in k_data if kd[0] > minT and kd[2] not in INVALID_OPTS]
					exonic_data = sorted([kd for kd in k_data if kd[2] in valid_opts],key=lambda X: X[2])
					exMax = [(obsN,avgN),(obsV,avgV)] 

			else:
				#if obsV < modT or avgN not in valid_opts: continue 
				val_avg = [kd for kd in avg_k if kd[2] in valid_opts] 				
				obs_avg = [kd for kd in obs_k if kd[2] in valid_opts] 				

				
				if TEST:
					if val_avg[0][2] == obs_avg[0][2]:
						obsC,avgC,nC,cC = val_avg[0] 

						if nC !=  exMax[0][0] and nC != exMax[0][1]: 
							try: 
								xO  = [vx for vx in obs_avg if vx[2] == exMax[0][0]][0]
								xA  = [vx for vx in val_avg if vx[2] == exMax[0][1]][0]

								obsD = obsC - xO[0] 
								valD = avgC / (0.001+xA[1]) 
											
								pv1 = stats.ttest_ind(cC,xO[-1])[1]
								pv2 = stats.ttest_ind(cC,xA[-1])[1]

					 
								print gene,prefix,'DIFFS',t,nC,'|',obsC,avgC,obsD,valD,pv1,pv2
							except IndexError:
								badness = True



				pvs = [1] 
				for i in range(len(val_avg)):
					obsI,avgI,nI,cI = val_avg[i] 
					for j in range(i+1,len(val_avg)):
						obsJ,avgJ,nJ,cJ = avg_k[j] 

						tup,idx,pv = quickstat(val_avg[i],val_avg[j],maxPV) 
						if pv != 'NA': 
							tup_key[tup][idx].append(t) 
							pvs.append(pv) 

				valid_key[t] = sorted([kd for kd in k_data if kd[2] in valid_opts],key=lambda X: X[2])
				pv_lists.append((min(pvs),t))

		if len(valid_opts) > 1 and len(valid_key) > 1: 

			if TEST: 
				d_trans = []
				for A,B in [(A,B) for (A,B) in tup_key.values() if len(A) > 0 and len(B)>0]: d_trans.extend(A+B) 
				if len(d_trans)>0:
					print gene, prefix, 'DTRANS', len(d_trans) 




#			vk =  valid_key.keys()  
#			pass_trans,m = [], 0 
#			pv_trans = [tn for tn in sorted(pv_lists)]
#			if len(vk) > maxTrans: 
#				d_trans = [] 
#				for A,B in [(A,B) for (A,B) in tup_key.values() if len(A) > 0 and len(B)>0]: d_trans.extend(A+B) 
#				pass_trans = list(set(d_trans)) 
#				
#				while m < len(pv_trans) and len(pass_trans) < maxTrans: 
#					pass_trans.append(pv_trans[m][1]) 
#					m+=1 
#
#				valid_key = {kn: valid_key[kn] for kn in valid_key.keys() if kn in pass_trans}
				 

				
			
			xlen = 1+int(0.99+((len(valid_key))/5))
			if not TEST: 			
				isoplot.setLens(xLen=xlen,yLen=1)
				isoplot.add_transcript_boxes(gene,exonic_data,valid_key) 
				isoplot.finish(gene,k)


		else:
			if TEST:
				print gene,prefix,'not covered'

			continue 








def run_script(data_file,options):
	geneID = None
	group_key = dd(lambda: dd(list))
	index_key = dd(lambda: dd(list))
	sample_key = dd(lambda: {}) 
	for line in open(options.key): 
		line = line.split() 
		if line[0] == '---': 
 			headers = line
		else:
			s = line[0] 
			for i in range(1,len(headers)):
				sample_key[line[0]][headers[i]] = line[i] 


	for line in open(data_file):

		line = line.split() 
		if line[0] == '---': 
			samples = line[1::] 
			for h in headers[1::]:
				h_key = dd(list) 
				i_key = dd(list) 
				for i in range(len(samples)): 
					if samples[i] in sample_key:
						hs = sample_key[samples[i]][h]
						h_key[hs].append(samples[i]) 
						i_key[hs].append(i) 
				group_key[h] = h_key
				index_key[h] = i_key 
			
		else:
			
                        gene,f = line[0].split('@')
                        if f == 'TOTAL': continue
		#	if gene.split(';')[-1] not in ["IP6K2","GTF3A","MACF1","PDE4DIP","SPTAN1"]: continue 
 
#			if geneID == None and gene != 'ENSG00000197694;chr9;SPTAN1': continue 

                        if gene != geneID: 

                                if geneID != None:
                                        plot_iso_data(geneID,tran_key,options,data_file) 

                                geneID = gene   
				tran_key = {} 
			
			cnts,cnt_key = [float(x) for x in line[1::]],{} 

			for h in index_key:
				cnt_key[h] = {k: [cnts[j] for j in index_key[h][k]] for k in index_key[h].keys()}


				

			tran_key[f] = cnt_key 





			
	sys.exit() 
	#['LOC', 'rad', 'MEGASIZE', 'loc', 'SURELOC', 'vBIO', 'ORIGSIZE', 'near', 'relsize', 'rel', 'vGW', 'data', 'size']  
















if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args[0],options)	














