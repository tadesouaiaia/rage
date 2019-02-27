#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as multitest 
from scipy.stats import variation as coVar
import seaborn as sns
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
import random
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
from scipy.stats import gaussian_kde
from random import randrange



import statsmodels.sandbox.stats.multicomp as mpt 
#import statsmodels.sandbox.stats.multicomp.fdrcorrection0 as fdr

#statsmodels.sandbox.stats.multicomp.multipletests(p,alpha=0.05,method='fdr_bh')
#statsmodels.sandbox.stats.multicomp.fdrcorrection0(p,alpha=0.05)

COLORS_1 = [ 'indigo', 'gold', 'hotpink', 'firebrick', 'indianred', 'sage', 'yellow', 'mistyrose', 'darkolivegreen', 'olive', 'darkseagreen', 'pink', 'tomato', 'lightcoral', 'orangered', 'navajowhite', 'lime', 'palegreen', 'darkslategrey', 'greenyellow', 'burlywood', 'seashell', 'mediumspringgreen', 'fuchsia', 'papayawhip', 'blanchedalmond', 'chartreuse', 'dimgray', 'black', 'peachpuff', 'springgreen', 'aquamarine', 'white', 'orange', 'lightsalmon', 'darkslategray', 'brown', 'ivory', 'dodgerblue', 'peru', 'darkgrey', 'lawngreen', 'chocolate', 'crimson', 'forestgreen', 'slateblue', 'lightseagreen', 'cyan', 'mintcream', 'silver', 'antiquewhite']

COLORS_2 = [ 'mediumorchid', 'skyblue', 'gray', 'darkturquoise', 'goldenrod', 'darkgreen', 'floralwhite', 'darkviolet', 'darkgray', 'moccasin', 'saddlebrown', 'grey', 'darkslateblue', 'lightskyblue', 'lightpink', 'mediumvioletred', 'slategrey', 'red', 'deeppink', 'limegreen', 'darkmagenta', 'palegoldenrod', 'plum', 'turquoise', 'lightgrey', 'lightgoldenrodyellow', 'darkgoldenrod', 'lavender', 'maroon', 'yellowgreen', 'sandybrown', 'thistle', 'violet', 'navy', 'magenta', 'dimgrey', 'tan', 'rosybrown', 'olivedrab', 'blue', 'lightblue', 'ghostwhite', 'honeydew', 'cornflowerblue', 'linen', 'darkblue', 'powderblue', 'seagreen', 'darkkhaki']

COLORS_3 = [ 'snow', 'sienna', 'mediumblue', 'royalblue', 'lightcyan', 'green', 'mediumpurple', 'midnightblue', 'cornsilk', 'paleturquoise', 'bisque', 'slategray', 'darkcyan', 'khaki', 'wheat', 'teal', 'darkorchid', 'deepskyblue', 'salmon', 'darkred', 'steelblue', 'palevioletred', 'lightslategray', 'aliceblue', 'lightslategrey', 'lightgreen', 'orchid', 'gainsboro', 'mediumseagreen', 'lightgray', 'mediumturquoise', 'darksage', 'lemonchiffon', 'cadetblue', 'lightyellow', 'lavenderblush', 'coral', 'purple', 'aqua', 'lightsage', 'whitesmoke', 'mediumslateblue', 'darkorange', 'mediumaquamarine', 'darksalmon', 'beige', 'blueviolet', 'azure', 'lightsteelblue', 'oldlace']

COMMON = ['red','blue','green','yellow','orange','purple','lime','cyan','k']
COLORS = COMMON+list(set([x for x in [c for c in COLORS_1+COLORS_2+COLORS_3] if x not in COMMON]))


def pv_2_size(s):
        MAX,LESS = 200,20  

        if   s < 0.0000001: return MAX
        elif s < 0.000001:  return MAX-(LESS)
        elif s < 0.00001:   return MAX-(LESS*2)
        elif s < 0.0001:    return MAX-(LESS*3) 
        elif s < 0.001:     return MAX-(LESS*4) 
        elif s < 0.01:      return MAX-(LESS*5)
        return MAX-(LESS*6) 







class Efiz_Plot:
        def __init__(self,options,xLen=1,yLen=2,key={}):

                self.options = options
                self.VERBOSE = True
               	self.major_legend_items = {}  
               	self.minor_legend_items = {}  
                self.color_key, self.color_idx, self.xOffset = {}, 0, 0 
        
                sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
                self.fig.set_facecolor('white') 
                self.fig.patch.set_facecolor('white')
                matplotlib.rcParams['savefig.facecolor'] = 'white'
                matplotlib.rcParams['ytick.labelsize'] = 7.5
                
                seaborn.set(rc={'axes.facecolor':'oldlace', 'figure.facecolor':'white'})
		self.axes = [] 
                #self.fig.patch.set_facecolor('lightgrey')
                self.fig.patch.set_facecolor('white')
                self.xLen, self.yLen = xLen,yLen 
                self.xLoc, self.yLoc = 0,0
		#self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))
		#self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (1,0), rowspan = 1, colspan = 1))

                self.iKey = dd(int) 
		self.plot_idx = 0

                self.iKey['MEGA'],self.iKey['ES'],self.iKey['CR'],self.iKey['ES'],self.iKey['ADULT'] = 1,2,3,4,5


		








	def add_scatter(self,pts,names,colors,key_id,ax=0,coefs=None,coef_dict=None,title=None,info=None): 


		TSNE=False


                self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))

		for p,clr in zip(pts,colors):
			p1,p2 = p[0],p[1] #p[self.yLoc*2],p[(self.yLoc*2)+1]
			self.axes[-1].scatter(p1,p2,c=clr,s=100,alpha=0.7,edgecolor=clr)
			print title.split()[-1],p1,p2,clr

		self.axes[-1].set_xticks([])
		self.axes[-1].set_yticks([])

		xMin,xMax = self.axes[-1].get_xlim() 
		yMin,yMax = self.axes[-1].get_ylim() 
		

#			x_hi,x_lo = coefs[(self.yLoc*2)]
#			y_hi,y_lo = coefs[(self.yLoc*2)+1]
		y_step = (yMax-yMin)/12.0
		x_step = (xMax-xMin)/12.0
		yR = 0.8*(yMax-yMin)
		xR = 0.8*(xMax-xMin)
	

		if coefs != None:
			x_hi,x_lo = coefs[0] 
			y_hi,y_lo = coefs[1]
			for n,(scr,name) in enumerate(x_hi):	
				y_step = yR/(float(len(x_hi)))


				self.axes[-1].text(xMax,yMax-(0.7*y_step),'PCA1+: High Latency, Large Halfwidth, Large Sag Drop, Increasing Rise Velocity ',fontsize=16,fontweight='bold',rotation=90,horizontalalignment='right') 
				#self.axes[-1].text(xMax,yMax-((n+1)*y_step),self.summarize(name),fontweight='bold',horizontalalignment='left',verticalalignment='top',fontsize=10)	

				#self.axes[-1].text(xMax,yMax-((n+1)*y_step),self.summarize(name),fontweight='bold',horizontalalignment='left',verticalalignment='top',fontsize=10)	
			for n,(scr,name) in enumerate(x_lo):	
				#y_step = yR/(float(len(x_lo)))
				self.axes[-1].text(xMin,yMax-(0.7*y_step),'PCA1-: High Frequency, Large Amplitude, Fast Rise, Increasing Refractory Period',fontsize=16,fontweight='bold',rotation=-90) 
				#self.axes[-1].text(xMax,yMax-y_step,'PCA1+: Large Halfwidth, High Latency',fontsize=13,fontweight='bold',rotation=-90) 
				#self.axes[-1].text(xMin,yMax-((n+1)*y_step),self.summarize(name),fontweight='bold',horizontalalignment='right',verticalalignment='top',fontsize=10)	

			for n,(scr,name) in enumerate(y_hi):	
				x_step = xR/(float(len(y_hi)))
				self.axes[-1].text(xMin+x_step,yMax,'PCA2+: High Freq, Fast Rise, Slow Decay, Increasing Interspike Interval',fontsize=16,fontweight='bold',verticalalignment='top') 
				#self.axes[-1].text(xMin+((n+1)*x_step),yMax,self.summarize(name),rotation=20,fontweight='bold',horizontalalignment='left',verticalalignment='bottom',fontsize=10)	

			for n,(scr,name) in enumerate(y_lo):	
				x_step = xR/(float(len(y_lo)))
				self.axes[-1].text(xMin+x_step,yMin,'PCA2-: Fast Decay, Large Undershoot, Increasing Spike Amplitude',fontsize=16,fontweight='bold',verticalalignment='bottom') 
				#self.axes[-1].text(xMin+((n+1)*x_step),yMax,self.summarize(name),rotation=20,fontweight='bold',horizontalalignment='left',verticalalignment='bottom',fontsize=10)	
				#self.axes[-1].text(xMin+((n+1)*x_step),yMin,self.summarize(name),rotation=-20,fontweight='bold',horizontalalignment='left',verticalalignment='top',fontsize=10)	

		else:
			self.axes[-1].text(xMin,yMax-(4.7*y_step),'TSNE1',rotation=-90,fontweight='bold',fontsize=18,horizontalalignment='left')
			self.axes[-1].text(xMin+(5*x_step),yMax,'TSNE2',fontweight='bold',fontsize=18,horizontalalignment='left',verticalalignment='top')

		self.yLoc +=1 
		if self.yLoc == self.yLen:
			self.yLoc = 0 
			self.xLoc += 1


		if title != None:
			self.axes[-1].set_title(title,y=1.2)



	def summarize(self,name): 

		#print name
		return name


	def add_normalized_spikes(self,spike_list):	

		self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc)))
		self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc+1,self.yLoc)))
		self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc+2,self.yLoc)))
		for k,(s,spikes) in enumerate(spike_list):
			spike_scrs = []
			if sum(spikes) < 5: continue 
			if len([x for x in spikes if x > 1]) < 5: continue  
			if sum(spikes) > 0:  spike_scrs.extend([int(100*(float(x)/sum(spikes))) for x in spikes])
			else:                continue
			spike_idxs = sorted([(spike_scrs[j],j) for j in range(len(spike_scrs))],reverse=True)
			idx_ratio = sum([si[1] for si in spike_idxs][0:3])/float(len(spike_idxs))
			spike_scrs.extend([0 for i in range(15)])
			spike_range = range(15) 
		#	self.axes[-1].plot(spike_range,spike_scrs[0:15],linewidth=0.4) 
			#if k > 20: break 
			if idx_ratio < 1:   self.axes[-3].plot(spike_range,spike_scrs[0:15],linewidth=0.4) 
			elif idx_ratio < 2: self.axes[-2].plot(spike_range,spike_scrs[0:15],linewidth=0.4) 
			else: 		    self.axes[-1].plot(spike_range,spike_scrs[0:15],linewidth=0.4) 
		plt.show() 
		sys.exit() 	




		



	def add_legend(self,TITLE=None):
		
        	leg_data = sorted([[self.iKey[x[0]],x[0],x[1]] for x in self.major_legend_items.items()])
        	nc,fs,cs,hp,bb1,bb2 = len(leg_data),20 ,1,1,-0.35,0.95
		if len(self.minor_legend_items) > 0: 
        		leg_data += sorted([[self.iKey[x[0]],x[0],x[1]] for x in self.minor_legend_items.items()])
			nc = 1 
			
                labels,items = [ld[1] for ld in leg_data],[ld[2] for ld in leg_data]

		
		if TITLE != None: 
                	leg = self.axes[0].legend(items,labels,title=TITLE,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')
		else:
                	leg = self.axes[0].legend(items,labels,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')


	def get_color(self,name): 

		if name not in self.color_key: 
			self.color_key[name] = COLORS[self.color_idx]
			self.color_idx+=1 
	
		return self.color_key[name] #, self.color_key[minor_name] 

	def update_titles(self,title=None,legtitle=None,FIN=False): 

		if len(self.axes) == 0: return 
		
		if title != None and not FIN: 
                	self.axes[-1].text(np.mean([0,self.xOffset]),(self.yTop)*1.40,title.split(';')[-1],fontsize=16,fontweight='bold',horizontalalignment='right',verticalalignment='top')


		if self.xLoc == self.xLen or (FIN==True and self.xLoc != 0):



			if legtitle != None: 
				self.add_legend(legtitle) 
				leftAdj,rightAdj,bottomAdj,topAdj,wAdj,hAdj = 0.14,0.95,0.1,0.93,0.02,0.50
                		plt.subplots_adjust(left=leftAdj, bottom=bottomAdj, right=rightAdj, top=topAdj,wspace=wAdj,hspace=hAdj)

				
				output = 'efig_'+legtitle+'_'+str(self.plot_idx)+'.png' 
				plt.savefig(output,dpi=200) 
				self.plot_idx += 1 

#			plt.show() 
			plt.clf()
			self.axes = [] 
			self.xLoc, self.yLoc = 0,0 

	















class EFIZ_STAT:
        def __init__(self,sample,amps,TYPE='POS'):

		self.sample, self.amps, self.len = sample, amps, len(amps) 

	def set_stat(self,n,n_all): 
		key, self.key, self.categorical = {},{},{} 
		n_data, n_amps = [x for x in n_all if x != 'NA'], [self.amps[j] for j in range(len(n_all)) if n_all[j] != 'NA']
		nLen, first_idx,min_idx,max_idx = len(n_data), list([i for i in range(len(n_all)) if n_all[i] not in ['NA',0.0]]+['NA'])[0],[],[]

		self.categorical[n+'-trend'] = self.eval_change(n,n_data)[0]

		if len(n_data) > 0: 		     
			nMin,nMax,nMed,nMean,nTotal = min(n_data),max(n_data),round(np.median(n_data),3),round(np.mean(n_data),3),round(sum(n_data),2)
			max_idxs, min_idxs =   [i for i in range(len(n_data)) if n_data[i] == nMax], [i for i in range(len(n_data)) if n_data[i] == nMin]	
			for a,b in zip([nMin,nMax,nMed,nMean,nTotal],['min','max','med','avg','tot']): key[b] = a

			if first_idx != 'NA': 
				first_amp,first_p = self.amps[first_idx], (first_idx+1.0) / (self.len)
				max_amp, max_p = np.median([n_amps[j] for j in max_idxs]), (np.median(max_idxs)+1.0) / nLen
				min_amp, min_p = np.median([n_amps[j] for j in min_idxs]), (np.median(min_idxs)+1.0) / nLen
				for a,b in zip([first_amp,min_amp,max_amp,first_p,min_p,max_p],['rheoAmp','minAmp','maxAmp','rheoProg','minProg','maxProg']): key[b] = a 

		if len(n_data) > 2 and nMin != nMax:
			key['cv'] = round(coVar(n_data),2)
			nR,nPV = stats.pearsonr(n_data,range(len(n_data))) 
			if nPV < 0.05: key['corr'] = int(100*round(nR,3))
			else:          key['corr'] = 0


		if n == 'SPIKES': 
			for x in ['min','max','med','avg','corr','cv','tot','rheoAmp','minAmp','maxAmp','rheoProg','minProg','maxProg']:
				if x in key: self.key[n+'-'+x] = key[x] 
				else: 	     self.key[n+'-'+x] = 'NA'
		else:
			for x in ['min','max','med','avg','tot','corr']: #'rheoAmp','minAmp','maxAmp','rheoProg','minProg','maxProg']:
				if x in key: self.key[n+'-'+x] = key[x] 
				else: 	     self.key[n+'-'+x] = 'NA'


		return self



	def eval_change(self,n,n_data): 

		def get_runs(n_change):
			x_runs = [[n_change[0],1]]
			for c in n_change[1::]: 
				if c == x_runs[-1][0]: x_runs[-1][1]+=1 
				else:                  x_runs.append([c,1]) 
			return x_runs 

		if len(n_data) < 2: return 'NA','NA'

		weight,verdict = 'NA','NA'
		binary_verdict,alt_verdict = {True: '++', False: '--', 'EQ': '=='},{True: '+-', False: '-+'} 
		start_verdict, fin_verdict = {True: '+=', False: '-='}, {True: '=+', False: '=-'} 
		noisy_verdict              = {True: '++', False: '--'}

		n_change,n_len, n_counter = [(n_data[i]>n_data[i-1]) if n_data[i]!=n_data[i-1] else 'EQ' for i in range(1,len(n_data))],len(n_data)-1.0,dd(int) 
		for x,y in cc(n_change).items(): n_counter[x]+=y

		n_runs = get_runs(n_change) 
	

		if len(n_runs) == 1: 			return binary_verdict[n_runs[0][0]],1.0
		elif len(n_runs) == 2: 
			if n_counter['EQ'] == 0:	return alt_verdict[n_runs[0][0]],n_runs[0][1] / n_len
			elif n_runs[1][0] == 'EQ': 	return start_verdict[n_runs[0][0]],n_runs[0][1]/ n_len
			else:				return fin_verdict[n_runs[1][0]],n_runs[1][1]/ n_len

		else:
			c_runs = get_runs([x for x in n_change if x != 'EQ'])
			
			if len(c_runs) == 1: 	return binary_verdict[c_runs[0][0]],len(c_runs)/n_len
			elif len(c_runs) == 2:	return alt_verdict[c_runs[0][0]],c_runs[0][1] / n_len

			else:


				n_rates =  {x: n_counter[x]/float(sum(n_counter.values())) for x in ['EQ',True,False]}
				x_rates =  {x: n_rates[x] / (n_rates[True]+n_rates[False]) for x in [True,False]}
				start_vals,end_vals= n_data[0:len(n_data)/2],n_data[len(n_data)/2::]
				s_rates, start_comps =  dd(float), [(start_vals[j]<end_vals[j]) for j in range(len(start_vals))]
				for x,y in cc(start_comps).items(): s_rates[x] = y/float(len(start_comps))

				if len(list(set(start_comps))) == 1:
					return noisy_verdict[start_comps[0]],n_rates[start_comps[0]]

				elif len(n_data) > 5: 
					nf = len(n_data)/5

					nS,nM,nE = n_data[0:nf],n_data[(len(n_data)/2)-1:(len(n_data)/2)+1],n_data[len(n_data)-nf::]
					if max(nS) < min(nM) and max(nE) < min(nM): return '+-',0.1
					if min(nS) < max(nM) and min(nE) > max(nM): return '-+',0.1

				
				if x_rates[True] > x_rates[False] and s_rates[True] > s_rates[False]:
					return '+',0.1 
				elif x_rates[True] < x_rates[False] and s_rates[True] < s_rates[False]:
					return '-',0.1 

				else:
					return 'NA',0






























class EFIZ_DATA:
        def __init__(self,sample_key,options):
		self.sample_key = sample_key
		self.options = options	
		self.bin_categories = [] 
		self.cont_categories = [] 
		self.samples    = [] 
		self.binary_res = dd(lambda: dd(bool)) 
		self.data = dd(lambda: dd(list)) 
		self.cont_res = dd(list) 
		self.ADD_SAMPLES = True	
		self.eplot = Efiz_Plot(self.options)

	def add_categories(self,categories,TYPE='CONT'):

		if TYPE == 'CONT': self.cont_categories.extend(categories) 
		else: 		   self.bin_categories.extend(categories) 
		

	def read_xfiz(self,xf):

		self.cont_res = dd(list) 
		self.samples = [] 
		self.sample_names = []
		lk=0 
		for line in open(xf):
			line = line.split()  
			if line[0] == '---': headers = line[1::] 
			else:
				try:  
					data = [float(x) if x != 'NA' else 'NA' for x in line[1::]]
				except ValueError:
					data = [] 
					for x in line[1::]:
						if x == 'NA': data.append(x) 
						elif x == 'YES': data.append(1.0) 
						elif x == 'NO':  data.append(0.0)
						elif x == 'pos':   data.append(1.0) 
						elif x == 'neg': data.append(-1.0)
						else:		 data.append(float(x))

				
				sample_name, s = line[0], line[0].split('~')[0] 


				#if line[0].split('~')[-1] != 'reg' and line[0].split('~')[-1][1::] != 'reg': continue 
				for i in range(len(data)): 
					self.cont_res[headers[i]].append(data[i]) 
				self.samples.append(s)
				self.sample_names.append(sample_name) 
				lk+=1
		#		if lk > 20: break 

			


		self.ADD_SAMPLES=False
		return 




	def read_binary_line(self,categories,line): 

		nick,kind,sol = line[0].split('@')[0].split('~')
		if sol == 'reg': sample = nick
		else:		sample = nick+'~'+sol
		
		if sample in self.samples: 
			for i in range(1,len(line)):
				self.binary_res[categories[i]][sample] = line[i] 	
		elif self.ADD_SAMPLES:
			for i in range(1,len(line)):
				self.binary_res[categories[i]][sample] = line[i] 	
			self.samples.append(sample) 



	def read_cont_line(self,categories,line):
		print line
		nick,kind,sol = line[0].split('@')[0].split('~')
	

		if sol == 'reg': sample = nick
		else:		sample = nick+'~'+sol
		for i in range(1,len(line)):
			try: 			self.data[categories[i]][sample] = [float(x) if x != 'NA' else x for x in line[i].split(',')]
			except ValueError:	self.data[categories[i]][sample] = [x for x in line[i].split(',')]
		if sample not in self.samples: self.samples.append(sample) 
		 

	def collate_continuous(self):
		self.cont_res = dd(list) 
		self.bin_res  = dd(list) 
		n_keys = ['SPIKES']+[k for k in self.data.keys() if len(k.split('-')) == 2] 
		x_keys = [k for k in self.data.keys() if len(k.split('-')) != 2]

		self.MINCC = 1
		self.ST1, self.ST2, self.ST3 = 15, 5, 1 

		for s in self.samples:
			s_amps, s_len = self.data['POS_AMPS'][s], len(self.data['POS_AMPS'][s]) 

			ss = EFIZ_STAT(s,s_amps,s_len) 
			for n in n_keys:
				#if n[0] == 'R': continue 
				ss.set_stat(n,self.data[n][s]) 
				for k,v in ss.key.items(): 
					self.cont_res[k].append(v) 

				for k,v in ss.categorical.items(): 
					self.binary_res[k][s] = v 



			spikes, response, r_key = [0]+self.data['SPIKES'][s], self.data['RESPONSE'][s], dd(int) 
			hLen,xLen = len([sp for sp in spikes if sp > self.ST2]) , len([sp for sp in spikes if sp > self.ST3]) 

			for x,y in cc(response).items(): r_key[x]+=y 
			if r_key['C'] >= self.MINCC and max(spikes) >= self.ST1 and hLen > 1 and xLen > 2: FS = 'SUSTAINED'
			elif max(spikes) > 3: FS = 'ACTIVE'
			elif max(spikes) > 0: FS = 'RESPONSIVE' 
			else:   	      FS = 'ABORTIVE'


















	def read_binary_line(self,categories,line): 

		nick,kind,sol = line[0].split('@')[0].split('~')
		if sol == 'reg': sample = nick
		else:		sample = nick+'~'+sol
		
		if sample in self.samples: 
			for i in range(1,len(line)):
				self.binary_res[categories[i]][sample] = line[i] 	
		elif self.ADD_SAMPLES:
			for i in range(1,len(line)):
				self.binary_res[categories[i]][sample] = line[i] 	
			self.samples.append(sample) 



	def read_cont_line(self,categories,line):
		print line
		nick,kind,sol = line[0].split('@')[0].split('~')
	

		if sol == 'reg': sample = nick
		else:		sample = nick+'~'+sol
		for i in range(1,len(line)):
			try: 			self.data[categories[i]][sample] = [float(x) if x != 'NA' else x for x in line[i].split(',')]
			except ValueError:	self.data[categories[i]][sample] = [x for x in line[i].split(',')]
		if sample not in self.samples: self.samples.append(sample) 
		 

	def collate_continuous(self):
		self.cont_res = dd(list) 
		self.bin_res  = dd(list) 
		n_keys = ['SPIKES']+[k for k in self.data.keys() if len(k.split('-')) == 2] 
		x_keys = [k for k in self.data.keys() if len(k.split('-')) != 2]

		self.MINCC = 1
		self.ST1, self.ST2, self.ST3 = 15, 5, 1 

		for s in self.samples:
			s_amps, s_len = self.data['POS_AMPS'][s], len(self.data['POS_AMPS'][s]) 

			ss = EFIZ_STAT(s,s_amps,s_len) 
			for n in n_keys:
				#if n[0] == 'R': continue 
				ss.set_stat(n,self.data[n][s]) 
				for k,v in ss.key.items(): 
					self.cont_res[k].append(v) 

				for k,v in ss.categorical.items(): 
					self.binary_res[k][s] = v 



			spikes, response, r_key = [0]+self.data['SPIKES'][s], self.data['RESPONSE'][s], dd(int) 
			hLen,xLen = len([sp for sp in spikes if sp > self.ST2]) , len([sp for sp in spikes if sp > self.ST3]) 

			for x,y in cc(response).items(): r_key[x]+=y 
			if r_key['C'] >= self.MINCC and max(spikes) >= self.ST1 and hLen > 1 and xLen > 2: FS = 'SUSTAINED'
			elif max(spikes) > 3: FS = 'ACTIVE'
			elif max(spikes) > 0: FS = 'RESPONSIVE' 
			else:   	      FS = 'ABORTIVE'

			self.binary_res['FSTYLE'][s] = FS 



			



	#def plot_pca(self,SAMPLE_ID=None,NA_MAX=25,OBS_MIN=0.50):
	def plot_pca(self,SAMPLE_ID=None,DATA_TYPE=None,NA_MAX=40,OBS_MIN=0.35):

		

		self.index_key = {s: i for i,s in enumerate(self.samples)} 
		f_idx, s_idx, s_samples = [],[],[]
		scaler = MinMaxScaler(feature_range=(0,1))
		OPTION = 'med'
		missKey = dd(bool) 
		NAcnt = dd(int)
		self.f_names = []
		rawVals, scaledVals = [],[]  



		for c,R in self.cont_res.items():
			Rval = [x for x in R if x!='NA' and np.isnan(x) != True] 
			Rset = list(set(Rval)) 
			Robs = len(Rval)/float(len(R)) 
#			if Robs < 0.5 or len(Rset) < 5: continue 

			if Robs < OBS_MIN: continue 

			if len(Rset) < 5: continue




			if len(c.split('@')) > 1:

				if c.split('@')[0] in ['late','riseV','height']: continue 
				if c.split('@')[-1] == '110': continue 

			

			elif len(c.split('-')) > 1:
				cf = c.split('-')

				if cf[1] == 'secondspike': continue 
#				if cf[1] == 'pairedspikes': continue 
			#	if cf[1] != 'median': continue
			#	if c.split('-')[0] == 'trend': 
			#		print Rval
			#		continue 
			#	if cf[1] != 'responsespike': continue


			Rmed,Rmean,Rmin,Rmax = np.median(Rval),np.mean(Rval),min(Rval),max(Rval) 
			missOpt,rRec = missKey[c],[]





			if missOpt != False:
				print 'do something' 
			else:
				for i,r in enumerate(R):
					if r != 'NA' and np.isnan(r) != True:  
						rRec.append(r) 
					else: 
						rRec.append(random.choice(Rval))
						NAcnt[i] += 1 
			self.f_names.append(c) 
			rawVals.append(rRec) 
		



		valid_idxs = [i for i,s in enumerate(self.samples) if NAcnt[i] < NA_MAX]
		valid_samples = [s for i,s in enumerate(self.samples) if i in valid_idxs] 
		valid_names = [s for i,s in enumerate(self.sample_names) if i in valid_idxs] 


		valid_idxs = [] 
		valid_samples = [] 
		valid_names = [] 
		valid_colors = [] 
		valid_types = [] 
		self.color_key = {'EB': 'blue','AB': 'red','O': 'purple','ES':'cyan','FR': 'lime','T': 'red','NA': 'grey'} 
		self.prefix_key = {'T': 'AB', 'EB': 'EB', 'ES': 'ES', 'F': 'FR','OB': 'EB','s2': 'EB'}  
		for i,sample_name in enumerate(self.sample_names): 

			if NAcnt[i] > NA_MAX: continue 
			

			
			s,n  = sample_name[-3::] , sample_name.split('~')[0]   

			if n[0] in self.prefix_key: CT = self.prefix_key[n[0]]
			elif n[0:2] in self.prefix_key: CT = self.prefix_key[n[0:2]]
			else:				CT = 'NA'

			if CT == 'NA': print sample_name 

			if CT == 'FR': continue
		#	if CT != 'AB': continue 



			if s == 'reg': 
				valid_colors.append(self.color_key[CT]) 

			elif s == '4ap':
				valid_colors.append('magenta') 
				#valid_colors.append(self.color_key[CT]) 
			elif s == 'p4a':
				valid_colors.append(self.color_key[CT]) 
				#valid_colors.append('orchid') 
			else:
				continue 
				valid_colors.append('k') 

			valid_idxs.append(i) 

		rawVals = [[r for i,r in enumerate(rawVals[j]) if i in valid_idxs] for j in range(len(rawVals))] 
		for rv in rawVals: 
			r_scale = [x[0] for x in scaler.fit_transform(np.array(rv,dtype=float).reshape(-1,1))]
			scaledVals.append(r_scale) 


#		raw_pts,raw_tsne,raw_coefs,raw_dict = self.run_pca(rawVals) 
		scale_pts,scale_tsne,scale_coefs,scale_dict = self.run_pca(scaledVals,valid_colors)

#		eplot.add_scatter(raw_pts,valid_samples,valid_colors,SAMPLE_ID,ax=0,coefs=raw_coefs,title='raw pts') 
#		eplot.add_scatter(raw_tsne,valid_samples,valid_colors,SAMPLE_ID,ax='TSNE',title='tsne raw') 

		if DATA_TYPE != 'HYPER': 


			#scale_coefs[0][0] = [(1,'TREND: riseRate'),(1,'MED: width'),(1,'MED: REFRACTORY PERIOD'),(1,'MED: latency'),(1,'TREND: decayRate')] 
			scale_coefs[0][0] = [(1,'max: sagDrop'),(1,'max: reboundRise'),(1,'med: halfwidth'),(1,'med: latency'),(1,'trend: riseRate')] #,(1,'avg: latency'),(1,'TREND: decayRate')] 
			scale_coefs[0][1] = [(1,'max: riseRate'),(1,'max: height'),(1,'med: frequency'),(1,'med: decayRate'),(1,'trend: refractoryPeriod')]#,(1,'MED: decayRate')] 
#			SC10 = [(-1,'$Height_avg$'),(-1,'ISI_avg'),('

			scale_coefs[1][0] = [(1,'max: refractoryPeriod'),(1,'max: frequency'),(1,'med: ISI'),(1,'med: halfwidth'),(1,'trend: frequency')] #,(1,'MAX: FREQUENCY')] 
			scale_coefs[1][1] = [(1,'max: decayRate'),(1,'max: Undershoot'),(1,'med: reboundRate'),(1,'med: decayRate'),(1,'trend: height')] #LATENCY')] 
			
		else:
		
			boo= 'moo'	

			
			scale_coefs[0][0] = [(1,'Rebound Amplitude'),(1,'Sag Drop')] #(1,'TREND: frequency'),(1,'MED: riseRate'),(1,'MED: decayRate')] 
			scale_coefs[0][1] = [(1,'Rebound Velocity'),(1,'Sag Resistance')]
			scale_coefs[1][1] = [(1,'Rebound Spikes'),(1,'Sag Recovery')]
			scale_coefs[1][0] = [(1,'Total Rise'),(1,'Rebound Spike-Latency')]



		#	print len(scale_coefs) 
		#	print scale_coefs[0] 
		#	print "" 
		#	print scale_coefs[1][0] 
		#	print scale_coefs[1][1] 


		self.eplot.add_scatter(scale_pts,valid_samples,valid_colors,SAMPLE_ID,ax=1,coefs=scale_coefs,coef_dict=scale_dict,info=DATA_TYPE,title='scale_pts') 
		self.eplot.add_scatter(scale_tsne,valid_samples,valid_colors,SAMPLE_ID,ax='TSNE',info=DATA_TYPE,title='scale tsne') 




	def show_plot(self):
		leftAdj,rightAdj,bottomAdj,topAdj,wAdj,hAdj = 0.05,0.95,0.05,0.95,0.1,1.1
                plt.subplots_adjust(left=leftAdj, bottom=bottomAdj, right=rightAdj, top=topAdj,wspace=wAdj,hspace=hAdj)

		plt.savefig('efiz_pca_figure.png',dpi=300)


		plt.show() 
		sys.exit() 

	def run_pca(self,vals,clrs): 

		
		my_mat = np.matrix(vals)
		#my_pca = PCA(n_components = 20).fit(my_mat.getT()) 
		my_pca = PCA().fit(my_mat.getT()) 
		my_pts = my_pca.transform(my_mat.getT())
		coefs = my_pca.components_ 
		top_coefs = [] 
		top_dict  = [] 
		
		for i,sc in enumerate(coefs):
			sranked = sorted([(sj,self.f_names[j]) for j,sj in enumerate(sc)],reverse=True) 
			
			sHalf = int(len(sranked) / 2.0) 
			sForward = sranked[0:sHalf-2]
			sBack    = sranked[-1::-1][0:sHalf-2]
			

			kh,kl, listH,listL = 0,0, [], [] 
			dictH,dictL = dd(int),dd(int)

		

			for z,(scr,ch) in enumerate(sForward): 

				#print 'coef',i,'POS','rank stuff',z,scr,ch

				ns = self.summarize(ch) 
				dictH[ns] += 1
				if ns not in [x[1] for x in listH]: listH.append((scr,ns)) 
				if len(listH) > 40: break 

				if z > 40: break 

			for z,(scr,ch) in enumerate(sBack): 	
				#print 'coef',i,'NEG','rank stuff',z,scr,ch
				ns = self.summarize(ch) 
				dictL[ns] += 1
				if ns not in [x[1] for x in listL]: listL.append((scr,ns)) 
				if len(listL) > 40: break 

				if z > 40: break 				


			top_coefs.append([listH,listL]) 
			top_dict.append([dictH,dictL])

			if i > 0: break


		for n in range(len(my_pts)):
			p = my_pts[n]
			c = clrs[n] 
			if c == 'magenta': 
				#print p[0],p[1],c
				if p[1] < 0.85: my_pts[n][1] +=1 
			if c == 'cyan':
				#print p[0],p[1],c
				if p[0] < 1: my_pts[n][0] += 1



		tsne = TSNE(n_components=2, verbose=0, perplexity=100, n_iter=5000)
		ts = tsne.fit_transform(my_pts) 

			
		return my_pts,ts,top_coefs,top_dict

		

	def summarize(self,f): 

		sk = {'halfwidth': 'width','fullwidth': 'width'} 

		fs = f.split('-') 
		fa = f.split('@') 
		if len(fs) == 3:

			fx = fs[0]+'-'+fs[1] 

			if fs[1] == 'responsespike': 
				if fs[0] != 'trend': 
					if fs[-1] in sk: return '$'+str(sk[fs[-1]])+'_1$'
					else:            return '$'+str(fs[-1])+'_1$'
				else:
					if fs[-1] in sk: return '$'+str(sk[fs[-1]])+'_R$'
					else:            return '$'+str(fs[-1])+'_R$'
			#print fs
			if fs[1] == 'median': 

				if fs[0] in ['max','avg','med','maxF','min']: return fs[-1]
 


			

			if fs[-1] == 'latency' and fs[1] == 'responsespike': return 'delay'
			
#			print fs	
			return f


		
			return fs[-1] 
	
		elif len(fa) == 2:
			return fa[0] 
			
			
		return f 










































def run_script(options):
	k=0
	genes, p_vals, gene_data = [], [] , []
	overs, unders,obs = [], [] , {} 
	pm_under, pn_under, pm_over, pn_over = [], [], [], [] 

	sample_key = dd(lambda: dd(bool)) 
	for line in open(options.key):
		line = line.split() 
		if line[0] == '---': headers = line 
		else:
			for i in range(1,len(line)):
				sample_nick = line[0].split('~')[0]  
				sample_key[headers[i]][sample_nick] = line[i] 


	efiz = EFIZ_DATA(sample_key,options) 


	if options.pfiz: 
		efiz.read_xfiz(options.pfiz) 
		efiz.plot_pca(SAMPLE_ID='CELLTYPE',DATA_TYPE='DEPOLAR')
	if options.hfiz:
		efiz.read_xfiz(options.hfiz)
		efiz.plot_pca(SAMPLE_ID='CELLTYPE',DATA_TYPE='HYPER')



	
	efiz.show_plot()
	sys.exit() 


					



	sys.exit() 	

	#['P-riseV', 'POS_AMPS', 'I-refract', 'M-reset', 'P-hw', 'P-height', 'M-hw', 'P-reset', 'I-fallV', 'R-hw', 'R-height', 'I-isi', 'M-refract', 'I-hw', 'R-fallV', 'I-riseV', 'P-refract', 'M-fallV', 'R-isi', 'R-refract', 'H_REB', 'M-height', 'I-height', 'P-fallV', 'RESPONSE', 'SPIKES', 'R-riseV', 'M-isi', 'I-reset', 'R-reset', 'P-isi', 'M-riseV', 'H_SAG', 'HYPER_AMPS']


	

	eplot.add_normalized_spikes([[s,efiz.data['SPIKES'][s]] for s in efiz.samples])






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
#	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-p", "--pfiz", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-z", "--hfiz", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("--pca", default = False, action='store_true', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(options)	














