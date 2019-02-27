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
 	merge = [(12,14),(1,5),(0,6)] 
	good = [0,1,4,5,6,8,9,12,14]
	for line in open(options.masks):
		line = line.split() 	
		m_pts = [float(p) for p in line[1].split(',')]
		if int(line[0]) not in [0,6,3]:
			masks[int(line[0])] = np.array(m_pts)
	#	if int(line[0]) in [4,8,9,1,12]: 
	#		masks[int(line[0])] = np.array(m_pts)

	fin = []
	p_masks= [] 
	for a,b in merge:
		fin.append(a) 
		fin.append(b) 
#		p_masks.append([np.mean([masks[a][i],masks[b][i]]) for i in range(len(masks[a]))])

#	for x,y in masks.items(): 
#		if x not in fin: p_masks.append(y) 

#	masks = {i: p_masks[i] for i in range(len(p_masks))}




	
	s_names, s_spikes  = [], [] 
	for line in open(spike_file): 
		line = line.split() 
		if len(line) == 2: name = line[0] 
		else: 		   name = "@".join(line[0:4]) 

		s_pts = [float(p) for p in line[-1].split(',')]
		s_names.append(name) 
		s_spikes.append(np.array(s_pts)) 

	spike_assignments = dd(bool) 	
	for i in range(len(s_spikes)): 
		
		s_name, pts = s_names[i], s_spikes[i] 
		mask_scrs = dd(list)  
		for m,MASK in masks.items(): 
			mask_scrs['spike'].append((np.linalg.norm(pts[200:-200]-MASK[200:-200]),m))
			mask_scrs['close'].append((np.linalg.norm(pts[150:-150]-MASK[150:-150]),m))
			mask_scrs['full'].append((np.linalg.norm(pts[0:500]-MASK[0:500]),m))
					
		spike,close,full = sorted(mask_scrs['spike']),  sorted(mask_scrs['close']), sorted(mask_scrs['full'])

		diffs = [spike[1][0] - spike[0][0], close[1][0] - close[0][0], full[1][0] - full[0][0]] 
		scrs = [spike[0][0],close[0][0],full[0][0]] 


		


		if spike[0][1] == close[0][1] and close[0][1] == full[0][-1]:
			


			

			#if sum(scrs) < 400 and max(scrs) < 140 and sum(diffs) > 200 and min(diffs) > 75: 
			#if sum(scrs) < 1000: 
			if sum(scrs) < 500 and max(scrs) < 250 and sum(diffs) > 200 and min(diffs) > 50: 
				spike_assignments[i] = spike[0][1]  
				
			
		if spike[0][1] != close[0][1] and close[0][1] != full[0][-1] and min(scrs) > 250:
			spike_assignments[i] = 99  
		elif sum(scrs) > 800 and min(scrs) >300:
			spike_assignments[i] = 99  
			print 'whoa'
#			print spike[0][1], close[0][1] 

#			print spike[0][1], close[0][1] 

#			print spike[0][1], close[0][1] 

#			print spike[0][1], close[0][1] 
	
	
	assignments = list(set(spike_assignments.values())) 
	xLen,yLen = len(assignments) / 3, 3+1 
	xLoc, yLoc = 0,0 
	axes = {} 
	for a in assignments: 
		axes[a] = plt.subplot2grid((xLen,yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)

		if a != 99: 		axes[a].plot(range(len(masks[a])),masks[a],color='k',linewidth=5)
		axes[a].set_ylim(-100,105)
		axes[a].set_title(a) 
		xLoc +=1 
		if xLoc == xLen:
			xLoc = 0 
			yLoc += 1 
			
	w = open( spike_file.split('.')[0]+'_spike_match.out', 'w') 
	u = open( spike_file.split('.')[0]+'_spike_miss.out', 'w') 
	z = open( spike_file.split('.')[0]+'_spike_far.out', 'w') 

	match_key = dd(list) 
	match,missed = 0,0
	for i in range(len(s_spikes)): 

		n, pts = s_names[i], s_spikes[i] 
		if n[0] == 'O' or not spike_assignments[i]: 
			u.write('%s %s\n' % (n, ",".join([str(s) for s in pts])))
			missed += 1 
		elif spike_assignments[i] == 99: 
			z.write('%s %s\n' % (n, ",".join([str(s) for s in pts])))
			missed += 1 
		else:
			w.write('%s %s %s\n' % (n,spike_assignments[i], ",".join([str(s) for s in pts])))

			a = spike_assignments[i] 
			if n[0] == 'T': clr = 'green' 
			elif n[0:2] == 'EB': clr = 'red'
			elif n[0:2] == 'ES': clr = 'cyan'
			else:		     clr = 'purple'

			axes[a].plot(range(len(pts)), pts,color=clr,linewidth=0.5)
			match += 1
			match_key[a].append(i) 

	new_masks = {} 
	for m in match_key.keys():
		if m != 99: 
			print m, len(match_key[m]),match_key[m][0]  
			big_m = [[s_spikes[x][j] for x in match_key[m]] for j in range(len(s_spikes[match_key[m][0]]))][0:500]
			bm = [] 	
			lowA,lowB,hiA,hiB = [],[],[],[]
			for pj,pt in enumerate(big_m):
				if m == 8: 
					points = np.array(sorted(pt[0:-5]))  
				else:
					points = np.array(sorted(pt[1:-1]))  
				points = points[:,None] 
				median = np.median(points, axis=0)
				diff = np.sum((points - median)**2, axis=-1)
				diff = np.sqrt(diff)
				med_abs_deviation = np.median(diff)
				modified_z_score = 0.6745 * diff / med_abs_deviation

				if m == 8:	
					lowB.append(np.mean(pt) - 2*med_abs_deviation) 
					hiB.append(np.mean(pt) + 2*med_abs_deviation) 
					if pj > 470: 
						pt = np.mean(sorted(pt)[0:-22]) 
						lowA.append(np.mean(pt) - 1.6*med_abs_deviation) 
						hiA.append(np.mean(pt) + 0.6* med_abs_deviation) 
					elif pj > 455: 
						pt = np.mean(sorted(pt)[0:-20]) 
						lowA.append(np.mean(pt) - 1.5*med_abs_deviation) 
						hiA.append(np.mean(pt) + 0.6* med_abs_deviation) 
					elif pj > 430: 
						pt = np.mean(sorted(pt)[0:-15]) 
						lowA.append(np.mean(pt) - 1.4*med_abs_deviation) 
						hiA.append(np.mean(pt) + 0.7* med_abs_deviation) 
					elif pj > 380: 
						pt = np.mean(sorted(pt)[0:-12]) 
						lowA.append(np.mean(pt) - 1.3*med_abs_deviation) 
						hiA.append(np.mean(pt) + 0.8* med_abs_deviation) 
					elif pj < 320:  
						lowA.append(np.mean(pt) - 1.5*med_abs_deviation) 
						hiA.append(np.mean(pt) + 0.8* med_abs_deviation) 
					elif pj < 300:         
						lowA.append(np.mean(pt) - med_abs_deviation) 
						hiA.append(np.mean(pt) +  1.2* med_abs_deviation) 
					else:
						lowA.append(np.mean(pt) - 1.2*med_abs_deviation) 
						hiA.append(np.mean(pt) +   med_abs_deviation) 


				else:		
					lowA.append(np.mean(pt) - med_abs_deviation) 
					hiA.append(np.mean(pt) + med_abs_deviation) 
					lowB.append(np.mean(pt) - 2*med_abs_deviation) 
					hiB.append(np.mean(pt) + 2*med_abs_deviation) 
				bm.append(np.mean(pt)) 
			axes[m].plot(range(len(lowA)),lowA,color='orange',linewidth=2)
			axes[m].plot(range(len(lowA)),hiA,color='orange',linewidth=2)

			new_masks[m] = [bm, lowA, hiA]



	wm = open(spike_file.split('.')[0]+'_new_masks.MASKS','w') 

	for m in new_masks: 
		a,b,c = new_masks[m]

		wm.write('%s %s %s %s\n' % (m,",".join([str(s) for s in a]),",".join([str(s) for s in b]),",".join([str(s) for s in c])))


	print 'match/miss',match,missed

	plt.savefig(spike_file.split('.')[0]+'spike_mask_match.png') 

	plt.show() 









				






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-m", "--masks", default = None, type='string', help="Output Filename Prefix")
	parser.add_option( "--show", default = False,action='store_true' , help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	scan_key(args[0],options)













