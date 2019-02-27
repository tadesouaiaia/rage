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
        def __init__(self,options,xLen=6,yLen=4,key={}):

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
                
                seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
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













	def add_scatter(self,pts,key,names,ax=0): 

		for p,n in zip(pts,names):

			cl,cf = key['CL'][n], key['CF'][n] 
			if n[0:2] == 'EB': 
				clr = 'green'
				if cl == 'mega': clr = 'purple' 
				if cf == 'MF':   clr = 'lime'  
			elif n[0:2] == 'ES': clr = 'cyan' 
			elif n[0:2] == 'OB': clr = 'black'
			elif n[0] == 'T': clr = 'red' 
			else:
				clr = 'orange'
			if ax == 0: p1,p2 = p[0],p[1] 
			elif ax == 1: p1,p2 = p[2],p[3] 

			self.axes[ax].scatter(p1,p2,c=clr) 


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

	


	def add_grouped_boxes(self,group_id,group_key): 


		my_key = dd(lambda: dd(list)) 
		for k in group_key.keys():
			major,minor =  k.split('&') 
			if minor in ['NA',False,'False'] or major in ['NA',False,'False']: continue 
			my_key[major][minor] = group_key[k]	







		
		major_keys, minor_keys = my_key.keys(), sorted(list(set([a for b in [my_key[k].keys() for k in my_key.keys()] for a in b])))

		



		self.scale_key, self.xOffset, self.yTop = dd(list) , 0 , 0  
                self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))
		
		for m in major_keys: 
			box_data = [] 
			for minor in minor_keys: 
				box_data.append(my_key[m][minor]) 
			self.make_minor_boxes(m,minor_keys,box_data)
		
		self.axes[-1].set_xlim(0-self.width,self.xOffset)  
                #self.axes[-1].text(np.mean([0,self.xOffset]),(self.yTop)-0.1,gene.split(';')[-1],fontsize=16,fontweight='bold',horizontalalignment='left',verticalalignment='top')
                self.axes[-1].axis('off') 

		self.yLoc += 1 
		if self.yLoc == self.yLen: 
			self.yLoc = 0 
			self.xLoc +=1





	def make_minor_boxes(self,major_id,opts,vals): 

		#self.ax.bar(xLoc+xOffset,X_key[X][s],width=1,bottom=0,color=a,ec=b,hatch='*')
		#items.append(Rect((0,0),1,1,fc=a,ec=b,hatch='*'))

		major_clrs,minor_clrs,global_sigs,local_sigs,jitters =[self.get_color(major_id) for X in opts],[self.get_color(X) for X in opts], [1 for X in opts],[1 for X in opts],[] 
		self.pos, self.width = [self.xOffset+0.20*x for x in range(len(opts))], 0.16

		all_vals = [a for b in vals for a in b]


                try: 
			valMin,valMean,valMax, valMeans = min(all_vals), np.mean(all_vals), max(all_vals)*1.25, [np.mean(v) for v in vals]
		except ValueError:
                	valMin,valMean,valMax, valMeans = 0,0,0,[0 for v in vals]
				


		for i in range(len(vals)):

                	other_vals = [a for b in [vals[k] for k in range(len(vals)) if k != i] for a in b] 
                	pv = stats.ttest_ind(vals[i],other_vals)[-1] 
                        if valMeans[i] > valMean and pv < global_sigs[i]: global_sigs[i] = pv
        
                        jitters.append([[self.pos[i] + np.random.normal(0,0.02) for m in range(len(vals[i]))],[v*np.random.normal(1,0.005) for v in vals[i]]])
                        for j in range(i+1,len(vals)): 
                                pv = stats.ttest_ind(vals[i],vals[j])[-1] 
                                if valMeans[i] > valMeans[j] and pv < local_sigs[i]: local_sigs[i] = pv 
                                elif valMeans[i] < valMeans[j] and pv < local_sigs[j]: local_sigs[j] = pv 



                self.major_legend_items[major_id] = Rect((0,0),1,1,fc=major_clrs[0]) #,ec=clr,linewidth=2,hatch='+')
		self.fill_boxes(opts,vals,minor_clrs,major_clrs,JITTERS=jitters,SIG_STARS=(global_sigs,local_sigs))
		if valMax > self.yTop: self.yTop = valMax
                self.xOffset = self.pos[-1]+2*self.width 

                self.axes[-1].text(np.mean(self.pos),0-(self.yTop)*0.1,major_id.split(';')[-1],fontsize=16,fontweight='bold',horizontalalignment='center',verticalalignment='top')


 

        def fill_boxes(self,opts,vals,colors,alt_colors,JITTERS=None,SIG_STARS=None):

                G_THRES,L_THRES = 0.05, 0.05
                self.bp = self.axes[-1].boxplot(vals,positions=self.pos, widths=self.width, patch_artist=True, showmeans=True,whis=0.7)
                for i,opt in enumerate(opts):


                        clr,ac = colors[i], alt_colors[i] 
                        self.bp['boxes'][i].set_edgecolor(clr) 
                        self.bp['boxes'][i].set_linewidth(1) 
                        plt.setp(self.bp['medians'][i], color=clr, linewidth=3)
                        plt.setp(self.bp['means'][i], marker='h',markersize=9,markerfacecolor=clr,markeredgecolor=ac) 
                        plt.setp(self.bp['caps'][(i*2)+1], color=clr,linewidth=1) 
                        plt.setp(self.bp['caps'][(i*2)], color=clr,linewidth=1) 
                        plt.setp(self.bp['whiskers'][i*2], color=clr,linewidth=1) 
                        plt.setp(self.bp['whiskers'][1+(i*2)], color=clr,linewidth=1)   
                        plt.setp(self.bp['fliers'][i], markerfacecolor=clr, markeredgecolor = ac, markeredgewidth=2,marker='s',markersize=2.0)           



                        
                for box, opt, clr, ac in zip(self.bp['boxes'], opts, colors, alt_colors):
			if clr == ac: 
                   		box.set_facecolor(clr)
                        	self.major_legend_items[opt] = Rect((0,0),1,1,fc=clr)
			else: 
				
				box.set(ec=clr,facecolor=ac,linewidth=4,hatch='+') 				
                        	self.minor_legend_items[opt] = Rect((0,0),1,1,fc='white',ec=clr,linewidth=2,hatch='+')


			box.set_alpha(0.5) 



                if JITTERS != None:
                        for i,(xJ,yJ) in enumerate(JITTERS):
                                plt.scatter(xJ, yJ, c=colors[i], alpha=0.7,s=4,zorder=9)
                if SIG_STARS != None:
                        for i,(g,l) in enumerate(zip(SIG_STARS[0],SIG_STARS[1])): 
                                if l < L_THRES: plt.scatter(self.pos[i],max(vals[i])*1.1,s=pv_2_size(l),marker='*',color='silver',edgecolor='black',linewidth=1)
                                if g < G_THRES: plt.scatter(self.pos[i],(max(vals[i])*1.1)+2,s=1.5*pv_2_size(l),marker='*',color='gold',edgecolor='black',linewidth=1)




















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
		

	def add_categories(self,categories,TYPE='CONT'):

		if TYPE == 'CONT': self.cont_categories.extend(categories) 
		else: 		   self.bin_categories.extend(categories) 
		

	def read_xfiz(self,xf):

		self.cont_res = dd(list) 
		self.samples = [] 
		for line in open(xf):
			line = line.split()  
			if line[0] == '---': headers = line[1::] 
			else: 
				data = [float(x) if x != 'NA' else 'NA' for x in line[1::]]

				s = line[0].split('~')[0] 
				if line[0].split('~')[-1] != 'reg': continue 
				for i in range(len(data)): 
					self.cont_res[headers[i]].append(data[i]) 
				self.samples.append(s) 

		return 

	def read_binary_line(self,categories,line): 

		nick,kind,sol = line[0].split('@')[0].split('~')
		if sol == 'reg': sample = nick
		else:		sample = nick+'~'+sol
		for i in range(1,len(line)):
			self.binary_res[categories[i]][sample] = line[i] 	
		if sample not in self.samples: self.samples.append(sample) 



	def read_cont_line(self,categories,line):
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



			
	def bin_dex(self,b):

		b_list = sorted([[round(np.mean(bv),2),round(float(len([vv for vv in bv if vv>0]))/len(bv),3),bv,bk] for bk,bv in self.b_key.items()],reverse=True)
		b_list = [bx for bx in b_list if len(bx[2])>10]
		if max([bl[1] for bl in b_list]) < 0.15: return 

		for x in range(len(b_list)): 
			xM,xO,xC,xN = b_list[x] 
			if len(xC) < 10 or xN == 'NA': return
			diff_cnts = [w for u in [b_key[vv] for vv in b_key.keys() if vv != xN] for w in u] 
			diff_mean, diff_obs = round(np.mean(diff_cnts),3), round(len([pp for pp in diff_cnts if pp > 0]) / float(len(diff_cnts)),3)
			pv  = stats.ttest_ind(xC,diff_cnts)[-1]
			if xM > diff_mean: FC = xM / (diff_mean+0.001) 						
			else: 		   FC = -1*(diff_mean / (xM+0.001)) 
			print self.gene,b,'GLOBAL',xN,'ELSE','|',len(xC),len(diff_cnts),'|',xO,diff_obs,'|',xM,diff_mean,'|',round(FC,5),round(pv,8) 
			for y in range(x+1,len(b_list)):
				
				yM,yO,yC,yN = b_list[y] 
				FC = xM/ (yM+0.001)
				print self.gene,b,'PAIRS',xN,yN,'|',len(xC),len(yC),'|',xO,yO,'|',xM,yM,'|',round(FC,5),round(pv,8)






	def group_bin_dex(self,b):



		valid_keys = [k for k in self.b_key.keys() if k.split('&')[0] not in ['NA'] and k.split('&')[1] not in ['NA']]
		#if float(len([x for x in ALL_VALS if x > 0])) / len(ALL_VALS) < 30: return False		
		#else: return True



		m_key = dd(lambda: {}) 
		for k,V in self.b_key.items():
			k1,k2 = k.split('&') 
			if k1 not in ['False','NA',False] and k2 not in ['False','NA',False]: 
				m_key[k1][k2] = V 

		
		global_tests,global_mids,local_tests,local_mids = [],[],[] ,[]
		
		
		for group in m_key.keys(): 

			b_list = sorted([[round(np.mean(bv),2),round(float(len([vv for vv in bv if vv>0]))/len(bv),3),bv,bk] for bk,bv in m_key[group].items()],reverse=True)			
			if max([bl[1] for bl in b_list]) < 0.15: continue

			for x in range(len(b_list)): 
				xM,xO,xC,xN = b_list[x]
 
				if len(xC) < 10 or xN == 'NA': continue
				diff_cnts = [w for u in [m_key[group][vv] for vv in m_key[group].keys() if vv != xN] for w in u] 
				if len(diff_cnts) < 10: continue 


				diff_mean, diff_obs = round(np.mean(diff_cnts),3), round(len([pp for pp in diff_cnts if pp > 0]) / float(len(diff_cnts)),3)
				pv  = stats.ttest_ind(xC,diff_cnts)[-1]
				if xM > diff_mean: FC,FS = xM / (diff_mean+0.001),'+'			
				else: 		   FC,FS = -1*(diff_mean / (xM+0.001)),'-' 
				#print self.gene,b,'GLOBAL',xN,'ELSE','|',len(xC),len(diff_cnts),'|',xO,diff_obs,'|',xM,diff_mean,'|',round(FC,5),round(pv,8) 
		
				if pv < 0.025: 
					global_tests.append((xN,xO,len(xC),group,FS,FC,pv))				
				if pv < 0.05 and FC>0:	
					global_mids.append(group+'-'+xN)
				for y in range(x+1,len(b_list)):
				
					yM,yO,yC,yN = b_list[y] 
					FC = xM/ (yM+0.001)
					pv  = stats.ttest_ind(xC,yC)[-1]
					#print self.gene,b,'PAIRS',xN,yN,'|',len(xC),len(yC),'|',xO,yO,'|',xM,yM,'|',round(FC,5),round(pv,8)

					if pv < 0.01: 
						local_tests.append((xN,xO,len(xC),yN,group,FC,pv))	
							

		if len(global_tests) > 0: 
			mpv = min([xx[-1] for xx in global_tests])
			mlen = max([xx[2] for xx in global_tests])
			mpo = max([xx[1] for xx in global_tests])
			return True,'GLOBAL',len(global_tests),",".join(sorted(global_mids)),mpo,mlen,mpv
		elif len(local_tests)>0:
			mpv = min([xx[-1] for xx in local_tests])
			mlen = max([xx[2] for xx in local_tests])
			mpo = max([xx[1] for xx in  local_tests])
			return True,'LOCAL',len(local_tests),",".join(sorted(global_mids)),mpo,mlen,mpv
			
		return False, False, False , False ,False,False,False



#
#		if len(global_tests) < 1 and len(local_tests) < 2: return False
#		return True		




	def add_cnts(self,cnt_file,SAMPLE_ID=None):

		
		self.cnt_plot = Efiz_Plot(self.options) 

		self.index_key = {s: i for i,s in enumerate(self.samples)} 
		f_idx, s_idx, s_samples = [],[],[]
		for line in open(cnt_file):
			line = line.split() 
			if line[0] == '---': 
				for i,s in enumerate(line[1::]): 
					if s in self.samples: 
						f_idx.append(i) 
						s_idx.append(self.index_key[s]) 
						s_samples.append(s) 
				continue 
					

			self.gene, self.cnts, self.b_key = line[0], [log(float(x)+1,2) for i,x in enumerate(line[1::]) if i in f_idx], dd(list) 
		
			self.obs = len([x for x in self.cnts if x > 0]) /float(len(self.cnts))

		#	if self.obs < 0.4: continue
	
			for b in self.binary_res.keys():
				if b[0] == 'X': continue 

				self.b_key = dd(list)  
				for c,s in zip(self.cnts,s_samples):
					if SAMPLE_ID!=None:	self.b_key[str(self.sample_key[SAMPLE_ID][s])+'&'+self.binary_res[b][s]].append(c) 
					else: 			self.b_key[self.binary_res[b][s]].append(c) 

				if SAMPLE_ID == None: self.bin_dex(b) 
				else:
					self.cnt_plot.add_grouped_boxes(b,self.b_key) 
			
			self.cnt_plot.add_title(self.gene) 
			self.cnt_plot.add_legend(TITLE=b) 
			plt.show() 
			sys.exit() 


 
			for c,R in self.cont_res.items():
				s_vals = [R[j] for j in s_idx] 	
				v_cnts,v_vals = [self.cnts[k] for k in range(len(s_vals)) if s_vals[k] != 'NA'],[s_vals[k] for k in range(len(s_vals)) if s_vals[k] != 'NA'] 
				PR,pv = stats.pearsonr(v_cnts,v_vals) 
				PS,sv = stats.spearmanr(v_cnts,v_vals) 
				if pv < 0.1 or sv < 0.1:  
					print self.gene,c,round(PR,4),round(PS,4),'|',len(v_cnts),'|',round(pv,8),round(sv,8)







	def plot_cnts(self,cnt_file,SAMPLE_ID=None):

		

		self.index_key = {s: i for i,s in enumerate(self.samples)} 
		f_idx, s_idx, s_samples = [],[],[]
		lk=0 

		self.genes, self.counts, self.obs = [], [], [] 
		for line in open(cnt_file):
			line = line.split() 
			if line[0] == '---': 
				for i,s in enumerate(line[1::]): 
					if s in self.samples: 
						f_idx.append(i) 
						s_idx.append(self.index_key[s]) 
						s_samples.append(s) 
				continue 
					

			self.gene, self.cnts, self.b_key = line[0], [log(float(x)+1,2) for i,x in enumerate(line[1::]) if i in f_idx], dd(list) 
			obs = len([x for x in self.cnts if x > 0]) /float(len(self.cnts))


			#if self.gene.split('@')[-1] in ['TOTAL','INTRONIC','EXONIC']: continue 
			#if len(self.gene.split('.')) > 1: continue 
			#if len(self.gene.split('orf')) > 1: continue 
			#if len(self.gene.split('-'))>1:     continue 
			if obs < 0.15: continue


			self.genes.append(line[0]) 
			self.counts.append(self.cnts) 
			self.obs.append(obs) 
	
			lk+=1
#			if lk>100: break 



	
		for b in self.binary_res.keys():


			if b[0] == 'X': continue 
			self.cnt_plot = Efiz_Plot(self.options) 
		
			for i in range(len(self.genes)): 
				
				#if self.genes[i] != 'ENSG00000151320;chr14;AKAP6': continue
				self.gene, self.cnts, obs = self.genes[i], self.counts[i], self.obs[i] 
				gene_nick = self.gene.split(';')[-1][0:3]  
				self.b_key = dd(list)  
				for c,s in zip(self.cnts,s_samples):
					if SAMPLE_ID!=None:	self.b_key[str(self.sample_key[SAMPLE_ID][s])+'&'+self.binary_res[b][s]].append(c) 
					else: 			self.b_key[self.binary_res[b][s]].append(c) 


				
				gbool, gtype, glen, gstr, gobs,xlen, gpv = self.group_bin_dex(b) 
				if gene_nick not in ['KCN','CAC','SLC','HCN']:
					if not gbool or gobs < 0.33: continue 
					if gpv > 0.005 and glen < 2: continue 
					if xlen < 10: 		     continue 
				elif not gbool: continue 
				

				print self.gene,self.cnt_plot.plot_idx, b , gbool, gtype, glen, gstr,gobs,xlen, gpv 

				self.cnt_plot.add_grouped_boxes(b,self.b_key) 			
				self.cnt_plot.update_titles(self.gene,b) 




			self.cnt_plot.update_titles(self.gene,b,FIN=True) 

		for c,R in self.cont_res.items():


			self.cnt_plot = Efiz_Plot(self.options) 
			s_vals = [(R[j],s_samples[i],self.sample_key[SAMPLE_ID][s_samples[i]]) for i,j in enumerate(s_idx)] 	


			for i in range(len(self.genes)): 
				
				self.gene, self.cnts, obs = self.genes[i], self.counts[i], self.obs[i] 
				gene_nick = self.gene.split(';')[-1][0:3]  
				c_key, self.b_key = dd(list),  dd(list)  
				for j,v in enumerate(s_vals): 
					if v[0] != 'NA' and v[2] not in [False,'NA','False']:
						c_key[v[2]].append([v[0],self.cnts[j]])
					
				for k,V in c_key.items():
					vLen = float(len(V)) 
					if vLen<10: continue 					
					else: 
						vL = int(vLen/3)
						
					V.sort() 
					for xLo,yLo in V[0:vL]:
						self.b_key[k+'&LOW'].append(yLo)
					for xMi,yMi in V[vL:2*vL]: 
						if xMi == xLo: self.b_key[k+'&LOW'].append(yMi)
						else:	       self.b_key[k+'&MID'].append(yMi)
					for xHi,yHi in V[2*vL::]:
						if xHi == xMi: self.b_key[k+'&MID'].append(yMi)
						else:	       self.b_key[k+'&HI'].append(yHi)
						
				gbool,gtype,glen,gstr,gobs,xlen,gpv = self.group_bin_dex('b') 
				if gene_nick not in ['KCN','CAC','SLC','HCN']:
					if not gbool or gobs < 0.33: continue 
					if gpv > 0.0015 and glen < 2: continue 
					if xlen < 10: 		     continue 
				elif not gbool: continue 


				print self.gene,self.cnt_plot.plot_idx, c , gbool, gtype, glen, gstr, gobs,xlen, gpv 
				self.cnt_plot.add_grouped_boxes(c,self.b_key) 			
				self.cnt_plot.update_titles(self.gene,c) 

			self.cnt_plot.update_titles(self.gene,c,FIN=True) 





















































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



	if options.xfiz: 
		efiz.read_xfiz(options.xfiz) 

	elif options.efiz: 	 
		for line in open(options.efiz): 
			line = line.split() 
			if line[0] == '---': categories = line  #efiz.add_categories(line)
			else:		     efiz.read_cont_line(categories,line) 
		efiz.add_categories(categories[1::],TYPE='CONT') 	
		efiz.collate_continuous() 

	elif options.bfiz: 
		for line in open(options.bfiz): 
			line = line.split() 
			if line[0] == '---': categories = line # efiz.add_categories(line,TYPE='BINARY') 
			else: 		     efiz.read_binary_line(categories,line)  
		efiz.add_categories(categories[1::],TYPE='BINARY') 





	if options.cnts: 
		efiz.plot_cnts(options.cnts,SAMPLE_ID='CELLTYPE') 
#		efiz.add_cnts(options.cnts,SAMPLE_ID='CELLTYPE') 
		#efiz.add_cnts(options.cnts) 

					



	sys.exit() 	

	#['P-riseV', 'POS_AMPS', 'I-refract', 'M-reset', 'P-hw', 'P-height', 'M-hw', 'P-reset', 'I-fallV', 'R-hw', 'R-height', 'I-isi', 'M-refract', 'I-hw', 'R-fallV', 'I-riseV', 'P-refract', 'M-fallV', 'R-isi', 'R-refract', 'H_REB', 'M-height', 'I-height', 'P-fallV', 'RESPONSE', 'SPIKES', 'R-riseV', 'M-isi', 'I-reset', 'R-reset', 'P-isi', 'M-riseV', 'H_SAG', 'HYPER_AMPS']


	

	eplot.add_normalized_spikes([[s,efiz.data['SPIKES'][s]] for s in efiz.samples])






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-c", "--cnts", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-e", "--efiz", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-b", "--bfiz", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-x", "--xfiz", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("--plot", default = False, action='store_true', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(options)	














