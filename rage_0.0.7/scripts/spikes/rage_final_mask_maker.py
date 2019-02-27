#!/usr/bin/env python

import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
import scipy.stats as stats
from scipy.stats import variation as coVar 

from random import random
import numpy as np
from random import randrange
from math import fabs
from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
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
from scipy.stats import chisquare
import numpy as np 
fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}




#matplotlib.rcParams['xtick.labelsize'] = 8.5
matplotlib.rcParams['ytick.labelsize'] = 6.5


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))




def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation

#    print get_colors(modified_z_score, plt.cm.jet)
    return modified_z_score > thresh, modified_z_score, get_colors(modified_z_score, plt.cm.seismic) 
#    return modified_z_score > thresh, modified_z_score, get_colors(modified_z_score, plt.cm.jet) 


















class mplot:
        def __init__(self,xLen,yLen,key={}):



		
		if 'r_key' in key: self.r_key = key['r_key']
		if 'p_key' in key: self.p_key = key['p_key']
	
		sns.set(rc={'axes.facecolor':'lightgrey', 'figure.facecolor':'white'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
		self.fig.set_facecolor('skyblue') 
    		self.fig.patch.set_facecolor('skyblue')
    		matplotlib.rcParams['savefig.facecolor'] = 'white'
		matplotlib.rcParams['ytick.labelsize'] = 7.5
		
		seaborn.set(rc={'axes.facecolor':'lightgrey', 'figure.facecolor':'white'})

                self.fig.patch.set_facecolor('white')
		self.xLen, self.yLen = xLen,yLen 
		self.xLoc, self.yLoc = 0,0 
		self.options = options

		if self.options and self.options.prefix != None: self.prefix = options.prefix
		else:						 self.prefix = 'model_plot'
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		self.pv_key = dd(int) 


	def add_pts(self,pts,names,coefs,varE,cords=[0,1]): 

		x,y = cords[0], cords[1]
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		for i,p in enumerate(pts): 
			

			size = 5 
			if len(names[i]) == 4: 
				name,current,order,spikes = names[i] 
			else:
				name,order,current,spikes,activity = names[i] 
			if name[0:2] == 'EB':   	color='red'
			elif name[0:2] == 'ES': 	color = 'cyan'
			elif name[0:1] =='T': 			    	color = 'green'
			elif name[0:2] == 'OB': 				color == 'purple'
			else: color= 'k' 

			text = name.split('~')[0]+'_'+str(order+1)+'/'+str(spikes)
			self.ax.text(p[x],p[y],text,color=color,fontsize=9,alpha=0.75) 
			
			#self.ax.text(p[x],p[y],text,color=color,fontsize=9,alpha=0.75) 
			size = 5 + int(current / 12.0)
			self.ax.scatter(p[x],p[y],color=color,s=size) 

		self.ax.set_xlabel('PCA'+str(x+1)+' = '+str(round(varE[x],3)))
		self.ax.set_ylabel('PCA'+str(y+1)+' = '+str(round(varE[y],3)))
			
		self.ax.set_xticks([]) 
		self.ax.set_yticks([]) 

 
		self.label_axis(coefs[x]) 
		self.label_axis(coefs[y],'Y') 



		self.xLoc +=1
		return 



	def label_axis(self,coefs,axis='X'): 
		x_pos =  sorted([(c[0],c[-1]) for c in coefs if c[1] > 0], reverse=True) 
		x_neg =  sorted([(c[0],c[-1]) for c in coefs if c[1] < 0], reverse=True) 
	
		gMax = 50.0
	
		x0,x1 = self.ax.get_xlim() 
		y0,y1 = self.ax.get_ylim() 
		
		yStep = (y1-y0)/gMax
		xStep = 2*((x1-x0)/gMax)
		fMax = 12
		fReal = 11
		fMin  = 6 
		maxItems = 100
		if axis.upper() != 'Y':

			

			gInit,xInit,yInit,k = gMax-3,x1,y1,0   
			for xS,xN in x_pos: 
				yInit-=yStep	
				fS = xS/float(x_pos[0][0])
				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit += (gMax-3)
					xInit += xStep
					yInit = y1-yStep 
					
 
				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin

				if fS > 0.90:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='red',verticalalignment='top',horizontalalignment='left') 
				elif fS > 0.80:	self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='orange',verticalalignment='top',horizontalalignment='left') 
				elif fS > 0.60:	self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='yellow',verticalalignment='top',horizontalalignment='left') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='green',verticalalignment='top',horizontalalignment='left') 
				else:		self.ax.text(xInit,yInit,xN,fontsize=fSize,alpha=fS,color='purple',verticalalignment='top',horizontalalignment='left') 


			gInit,xInit,yInit,k = gMax-3,x0,y1,0   
			for xS,xN in x_neg: 
				yInit-=yStep					
				fS = xS/float(x_neg[0][0])
				 

				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit+=(gMax-3)
					xInit -= (1.1*xStep) 
					yInit = y1-yStep 	
				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin
				
				
				if fS > 0.85:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='k',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.65:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='grey',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='purple',verticalalignment='top',horizontalalignment='right') 
				else:		self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='blue',verticalalignment='top',horizontalalignment='right') 

				


		else: 
			
			gInit,xInit,yInit,k = gMax*1.1,x0,y1,0 
			for xS,xN in x_pos: 
				xInit+=(yStep*1.3)	
				fS = xS/float(x_pos[0][0])
				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit += gMax
					xInit = x0+(2*yStep)
					yInit = y1+yStep 
					

				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin

				if fS > 0.85:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='k',verticalalignment='bottom',horizontalalignment='right') 
				elif fS > 0.65:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='grey',verticalalignment='bottom',horizontalalignment='right') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='purple',verticalalignment='bottom',horizontalalignment='right') 
				else:		self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='blue',verticalalignment='bottom',horizontalalignment='right') 

			gInit,xInit,yInit,k = gMax*1.1,x0,y0,0 
			for xS,xN in x_neg: 
				xInit+=(yStep*1.3) 	
				fS = xS/float(x_neg[0][0])
				fSize = int(fS*fMax) 
				k+=1 
				if k > maxItems: break 
				elif k > gInit:
					gInit += gMax
					xInit =  x0+(2*yStep) 
					yInit = y0-yStep 
				if   fSize > fReal: fSize = fReal 
				if   fSize < fMin: fSize = fMin

				if fS > 0.85:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='k',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.65:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='grey',verticalalignment='top',horizontalalignment='right') 
				elif fS > 0.40:	self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='purple',verticalalignment='top',horizontalalignment='right') 
				else:		self.ax.text(xInit,yInit,xN,fontweight='bold',fontsize=fSize,alpha=fS,color='blue',verticalalignment='top',horizontalalignment='right') 







































































def qc(name): 
	if name[0:2] == 'EB': return 'red'
	elif name[0:2] == 'ES': return 'cyan'
	elif name[0:2] == 'OB': return 'purple'
	elif name[0] == 'T': return 'green'
	else:			return 'k'

def run_pca(data,names): 

	#vals = [[log(x+1.0) for x in sv] for sv in s_vals] 
	fit_mat = np.matrix(data) 
	fit_pca = PCA(n_components = 5).fit(fit_mat)
	varE =  fit_pca.explained_variance_ratio_ 
	pts = fit_pca.transform(fit_mat)
	coefs = [sorted([(c*c,c,f) for (c,f) in zip(comps,range(len(comps)))],reverse=True) for comps in fit_pca.components_]



	plot = mplot(1,1) 
	plot.add_pts(pts,names,coefs,varE,cords=[0,1]) 	

	return pts,coefs












def scan_key(spike_file,options):

	masks = {}
	k=0
	for line in open(options.masks):
		line = line.split() 
		
		m_pts = [float(p) for p in line[1].split(',')]
		lo_pts = [float(p) for p in line[2].split(',')]
		hi_pts = [float(p) for p in line[3].split(',')]
#		if int(line[0]) in [12,8,9,2,6,1]: continue 
		if int(line[0]) in [11,7,6,13,5,8,10]: continue 
		
		masks[int(line[0])] = m_pts 


	s_data = dd(list) 
	for line in open(spike_file): 
		line = line.split() 
		if len(line) == 2: name = line[0] 
		else: 		   name = "@".join(line[0:4]) 

		cell = name.split('~')[0] 
		s_pts = [float(p) for p in line[-1].split(',')][0:500]
		s_data[cell].append((name,np.array(s_pts)))

	spike_assignments = dd(bool) 

	pair_scrs = dd(int) 
	match_grades = dd(lambda: dd(int))
	verdicts = {}  
	for cell in s_data.keys(): 
		mask_outs = dd(list) 
		mask_choices = [] 
		
		for name,pts in s_data[cell]: 
			mask_ranks = []  
			mask_rank = dd(list) 
			mask_val = dd(list) 
			mask_scrs = dd(list)  
			for m,MASK in masks.items(): 
				mask_scrs['spike'].append((np.linalg.norm(pts[200:400]-MASK[200:400]),m))
				mask_scrs['pre'].append((np.linalg.norm(pts[0:300]-MASK[0:300]),m))
				mask_scrs['post'].append((np.linalg.norm(pts[300:500]-MASK[300:500]),m))
				mask_scrs['full'].append((np.linalg.norm(pts-MASK),m))
				mask_ranks.append([mask_scrs['spike'][-1][0] , mask_scrs['pre'][-1][0], mask_scrs['post'][-1][0], mask_scrs['full'][-1][0],m]) 
				mask_val[m] = mask_ranks[-1][0:4]
			mr = sorted(mask_ranks) 
			for i,x in enumerate(sorted(mask_ranks)):			mask_rank[x[-1]].append(i) 
			for i,x in enumerate(sorted(mask_ranks,key=lambda X: X[1])):	mask_rank[x[-1]].append(i) 
			for i,x in enumerate(sorted(mask_ranks,key=lambda X: X[2])):	mask_rank[x[-1]].append(i) 
			for i,x in enumerate(sorted(mask_ranks,key=lambda X: X[3])):	mask_rank[x[-1]].append(i) 
				
			for m in masks: 
				mask_outs[m].append((np.mean(mask_rank[m]),np.mean(mask_val[m])))
			my_choice = []
			for m in masks.keys(): 

				if sum(mask_rank[m]) == 0 and sum(mask_val[m]) < 1000:  
					my_choice = [(0,m)] 
					break
				elif min(mask_rank[m]) > 2 and sum(mask_val[m]) > 3000: 
					continue 
				else: 
					my_choice.append((np.mean(mask_val[m]) * (1+np.mean(mask_rank[m])),np.mean(mask_val[m]),np.mean(mask_rank[m]),m))

			mask_choices.append(sorted(my_choice)) 


		rank_list, scr_list = [], [] 
		for m in mask_outs:
			rs,ss = np.median([x[0] for x in mask_outs[m]]),np.median([x[1] for x in mask_outs[m]]) 			
			rank_list.append((rs,m))
			scr_list.append((ss,m))
		rank_list.sort() 
		scr_list.sort() 
		VERDICT = None 
		tops, my_tops = [x[0][-1] for x in mask_choices], [] 
		if (rank_list[0][1] == scr_list[0][1]):
			
			CC = sorted(cc(tops).items(),key=lambda X: X[1],reverse=True) 
			if len(list(set(tops))) == 1: 
				VERDICT = tops  
				
			elif CC[0][1] / float(len(tops)) >= 0.90: 
				VERDICT = [CC[0][0] for t in tops]
			else:
				for i,c in enumerate(mask_choices):
					if len(c) == 1 and c[0][0] == 0:
						my_tops.append(c[0][-1]) 
					else:
						if c[0][-1] == CC[0][0]: my_tops.append(c[0][-1]) 
						elif (c[0][0]*3<c[1][0]) and (c[0][1]*2 < c[1][1]) and (c[0][2] * 2 < c[1][2]): my_tops.append(c[0][-1]) 
						elif (c[0][1] * 2 > c[1][1]) and (c[1][-1] == CC[0][0]):  my_tops.append(c[1][-1]) 
						elif CC[0][0] not in [c[0][-1],c[1][-1]] and (c[0][1] > 500 and c[0][2]) > 1: my_tops.append('NA') 				
						elif rank_list[0][1] == CC[0][0] and CC[0][0] in [c[0][-1], c[1][-1],c[2][-2]]: my_tops.append(CC[0][0]) 
						else: my_tops.append('NA') 
		else:
			for i,c in enumerate(mask_choices):
				if len(c) == 1: my_tops.append(c[0][-1]) 
				elif (c[0][0]*3<c[1][0]) and (c[0][1]*2 < c[1][1]) and (c[0][2] * 2 < c[1][2]): my_tops.append(c[0][-1]) 
				else: my_tops.append('NA') 

		if not VERDICT: VERDICT = my_tops 
		verdicts[cell] = VERDICT


	xLen,yLen = len(masks) / 3, 3 
	xLoc, yLoc = 0,0 
	axes = {} 
	id_key = {'NA': 'NA'} 	
	m_ids = ['A','B','C','D','E','F']
	for i,m in enumerate([1,4,3,12,9,2]): 
		MASK = masks[m] 
		if i == 0:
			MX = range(len(MASK))
		M_ID = m_ids[i] 
		ax = plt.subplot2grid((xLen,yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
		ax.plot(MX,MASK,color='k',linewidth=5,zorder=3,linestyle='--')
#		ax.set_title(str(m)+'      '+" ".join([str(s) for s in [match_grades[m][1],match_grades[m][2],match_grades[m][3]]])) 
		id_key[m] = M_ID 
		axes[m] = ax
		yLoc += 1
		if yLoc == yLen:
			yLoc =0 
			xLoc +=1 
#		xLoc +=1 
#		if xLoc == xLen:
#			xLoc = 0 
#			yLoc += 1
	collecter = dd(list) 
	MS = MX[25:-25]
	clr_counts = dd(lambda: dd(int)) 
	w = open( spike_file.split('.')[0]+'_spike_verdicts.out', 'w') 
	for CI,cell in enumerate(s_data.keys()): 
		verdict = verdicts[cell] 
		w.write('%s %s\n' % (cell," ".join([id_key[v] for v in verdict])))
		mode = sorted(cc(verdict).items(),key=lambda X: X[1],reverse=True)[0]
		clr = qc(cell)
		if clr == 'purple': continue  
		for j,(name,pts) in enumerate(s_data[cell]): 
			info, v = name.split('@'), verdict[j] 
			if v == 'NA': continue 
			dist = np.linalg.norm(pts[25:-75]-masks[v][25:-75])
			lt = id_key[v] 
			lw,alp = 1,0.3 
			if lt == 'E': lw = 1.75
			if lt == 'E': alp = 0.5
			p_max = max(pts[25:200])	
			if lt in ['C','B','D'] and p_max > 20: continue 
			if dist > 1500: continue 
			elif dist > 600:
				my_pts = [np.mean([pts[n],masks[v][n],masks[v][n],masks[v][n]]) for n in range(len(pts))]
				axes[v].plot(MS,my_pts[25:-25],color=clr,linewidth=lw,zorder=2,alpha=alp)
			elif dist > 500 or (dist > 300 and lt == 'E'): 
				my_pts = [np.mean([pts[n],masks[v][n],masks[v][n]]) for n in range(len(pts))]
				axes[v].plot(MS,my_pts[25:-25],color=clr,linewidth=lw,zorder=2,alpha=alp)
			elif dist > 250: 
				my_pts = [np.mean([pts[n],masks[v][n]]) for n in range(len(pts))]
				axes[v].plot(MS,my_pts[25:-25],color=clr,linewidth=lw,zorder=2,alpha=alp)
			else:	
				axes[v].plot(MS,pts[25:-25],color=clr,linewidth=1,zorder=2,alpha=alp)

			clr_counts[v][clr] += 1

#		if CI > 90: 	break 

	
	for i,m in enumerate([1,4,3,12,9,2]): 
	
		M_ID = id_key[m] 
		axes[m].set_ylim(-40,250)  
		axes[m].set_xlim(25,440)  

		xMin,xMax = axes[m].get_xlim() 

		dHeight = 110
		dX =      235
		wid =     35


		axes[m].plot([dX,400],[dHeight,dHeight],color='k',linewidth=3) 
		axes[m].text(200,dHeight+5,M_ID,fontsize=35,fontweight='bold') 

		eb,ab,es = clr_counts[m]['red'], clr_counts[m]['green'], clr_counts[m]['cyan']

		if M_ID in ['C','F']:
			eb,ab,es = 1+int(eb/10.0), 1+int(ab/10.0), 1+int(es/10.0) 
		elif M_ID != 'E':
			eb,ab,es = 1+int(eb/5.0), 1+int(ab/5.0), 1+int(es/5.0) 

		
		axes[m].bar(dX+20,ab,wid,dHeight,color='green') 
		axes[m].bar(dX+20+wid,eb,wid,dHeight,color='red') 
		axes[m].bar(dX+20+wid+wid,es,wid,dHeight,color='cyan') 
		axes[m].text(dX+20+wid,dHeight+eb,'yo') 

		axes[m].axis('off') 

	plt.suptitle('Spike Shape Analysis') 
	plt.subplots_adjust(left=0.05, bottom=0.05, right=0.92, top=0.95,wspace=0.1,hspace=0.1)
#	plt.savefig(spike_file.split('.')[0]+'_finalmaskfig.png')
	plt.show()
	sys.exit() 









				






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-m", "--masks", default = None, type='string', help="Output Filename Prefix")
	parser.add_option( "--show", default = False,action='store_true' , help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	scan_key(args[0],options)













