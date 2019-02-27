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
	s_names  = [] 
	s_spikes = [] 
	for line in open(spike_file): 
		line = line.split() 
		if len(line) == 2: 
			name = line[0].split('@')[0] 
			current,order,total = [int(x) for x in line[0].split('@')[1::]]
		else:
			name,current,order,total = line[0], int(line[1]), int(line[2]) , int(line[3]) 
		s_pts = [float(p) for p in line[-1].split(',')]
		s_names.append([name,current,order,total])  
		#s_spikes.append(np.array(s_pts[150:-50])) 
		s_spikes.append(np.array(s_pts)) 
	s_dists = [] 
	s_key = dd(list)  
	pair_key = {} 



	pca_pts, pca_coeffs = run_pca([s[100:-100] for s in s_spikes],s_names) 
	plt.savefig(spike_file.split('.')[0]+'_spike_pt_pca_100.png') 
	if options.show: plt.show() 

	pca_key = {} 
	for i in range(len(s_spikes)): 
		for j in range(i+1,len(s_spikes)): 
			dist = np.linalg.norm(s_spikes[i][100:-100]- s_spikes[j][100:-100])  
			dist = np.linalg.norm(s_spikes[i][25:-25]- s_spikes[j][25:-25])  

			pi1,pi2= pca_pts[i][0],pca_pts[i][1]
			pj1,pj2 = pca_pts[j][0],pca_pts[j][0] 
			pca_dist = ((pi1-pj1)*(pi1-pj1))+((pi2-pj2)*(pi2-pj2))
				

			s_key[i].append([dist,j]) 
			s_key[j].append([dist,i]) 
			pair_key[(i,j)] = dist
			pair_key[(j,i)] = dist  	
			s_dists.append([dist,(i,j)]) 

			pca_key[(i,j)] = pca_dist
			pca_key[(j,i)] = pca_dist  	


	s_dists.sort()
	lowV,hiV = s_dists[int(len(s_dists)*0.30)][0], s_dists[int(len(s_dists)*0.80)][0]

	s_ranks = [] 
	for i in s_key.keys():
		s_key[i].sort() 
		s_ranks.append((len([x for x in s_key[i] if x[0] < lowV]),i))


	s_ranks.sort(reverse=True) 
	
	xLen,yLen = 4,4 
	axes = [plt.subplot2grid((xLen,yLen), (0,0), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (1,0), rowspan = 1, colspan = 1)]
	axes.extend([plt.subplot2grid((xLen,yLen), (2,0), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (0,1), rowspan = 1, colspan = 1)])
	axes.extend([plt.subplot2grid((xLen,yLen), (1,1), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (2,1), rowspan = 1, colspan = 1)])
	axes.extend([plt.subplot2grid((xLen,yLen), (0,2), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (1,2), rowspan = 1, colspan = 1)])
	axes.extend([plt.subplot2grid((xLen,yLen), (3,0), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (3,1), rowspan = 1, colspan = 1)])
	axes.extend([plt.subplot2grid((xLen,yLen), (3,2), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (3,3), rowspan = 1, colspan = 1)])
	axes.extend([plt.subplot2grid((xLen,yLen), (0,3), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (1,3), rowspan = 1, colspan = 1)])
	axes.extend([plt.subplot2grid((xLen,yLen), (2,3), rowspan = 1, colspan = 1), plt.subplot2grid((xLen,yLen), (2,2), rowspan = 1, colspan = 1)])

	AX = 0 
	X = range(len(s_spikes[0])) 	
	
	fin = dd(bool) 	

	bin_dict = dd(list) 

	pca_cutoff =  np.median(pca_key.values())
	pca_cutoff =  np.percentile(pca_key.values(),25)

	pca_q = np.percentile(pca_key.values(),10) 

	for j,(scr,n) in enumerate(s_ranks): 


		if   s_names[n][0][0:2] == 'EB': clr = 'red'
		elif s_names[n][0][0:2] == 'ES': clr = 'cyan'
		else:			      clr = 'green'
	
		if j == 0: 
			MATCH=True
			axes[AX].plot(X,s_spikes[n],color=clr,linewidth=0.5) 
			axes[AX].set_title('AXES '+str(AX))
			bin_dict[AX].append(n) 

			fin[n] = True 
			for d,k in s_key[n]: 
				if d > lowV: break 
				if pca_key[(n,k)] > pca_cutoff: break 
				if not fin[k]: 
					axes[AX].plot(X,s_spikes[k],color=qc(s_names[k][0]),linewidth=0.5) 
					bin_dict[AX].append(k) 
				fin[k] = True 

		else:
			ps =  pair_key[(n,s_ranks[j-1][1])]
			
			
			if MATCH and ps < lowV and pca_key[(n,k)] < pca_cutoff:
				axes[AX].plot(X,s_spikes[n],color=clr,linewidth=0.5) 
				bin_dict[AX].append(n) 
				fin[n] = True
				for d,k in s_key[n]: 
					#if d > lowV: break 
					if d > lowV or pca_key[(n,k)] > pca_cutoff: break 
					if not fin[k]: 
						axes[AX].plot(X,s_spikes[k],color=qc(s_names[k][0]),linewidth=0.5) 
						bin_dict[AX].append(k) 
						fin[k] = True 
				
			else:
				MATCH= False 
				if ps > hiV: 
					
					AX+=1 
					axes[AX].plot(X,s_spikes[n],color=clr,linewidth=0.5) 
					axes[AX].set_title('AXES '+str(AX))
					bin_dict[AX].append(n) 
					fin[n] = True 
					for d,k in s_key[n]: 
						if d > lowV or pca_key[(n,k)] > pca_cutoff: break 
						if not fin[k]: 
							axes[AX].plot(X,s_spikes[k],color=qc(s_names[k][0]),linewidth=0.5) 
							bin_dict[AX].append(k) 
							fin[k] = True 

			if AX + 2== len(axes):  break 

	for ax in axes: ax.set_ylim(-50,100)

	total_len = 0 
	print ""
	

	plt.savefig(spike_file.split('.')[0]+'_spike_types_100offset.png') 
	if options.show:	plt.show() 	
	MASKS = [] 
	for AX in bin_dict.keys():
		N  =  bin_dict[AX]
		total_len += len(N)
		if total_len < 1: continue 
		if len(N) < 1: continue 
		mean_pos = [] 
		scrs = []
		p_scrs, d_scrs = [],[]  
		for i in range(len(N)):  
			n = N[i]
			pca_scr= sorted([pca_key[(n,m)] for m in N  if m!=n ])[0:-2]
			d_scr = sorted([pair_key[(n,m)] for m in N if m!=n])[0:-2]
			p_scrs.append((pca_scr,n)) 
			d_scrs.append((d_scr,n)) 		

		p_scrs.sort()
		d_scrs.sort()
		p_pass = [x[1] for x in p_scrs[0:int(len(p_scrs)*0.75)]]
		d_pass = [x[1] for x in d_scrs[0:int(len(d_scrs)*0.75)]]
		b_pass = [a for a in p_pass if a in d_pass] 

		if len(b_pass) > 3:
			big_spike =  [sorted([s_spikes[b][u] for b in b_pass]) for u in range(len(s_spikes[b_pass[0]]))] 
			if len(b_pass) < 6:	big_mean = [np.mean(bs[1:-1]) for bs in big_spike] 
			elif len(b_pass) < 12:	big_mean = [np.mean(bs[2:-2]) for bs in big_spike] 
			elif len(b_pass) < 24:  big_mean = [np.mean(bs[4:-4]) for bs in big_spike] 		
			elif len(b_pass) < 30:  big_mean = [np.mean(bs[6:-6]) for bs in big_spike] 		
			else:	  		big_mean = [np.mean(bs[8:-8]) for bs in big_spike] 		
			MASKS.append((AX,big_mean)) 
		else:
			big_spike =  [sorted([s_spikes[b][u] for b in N]) for u in range(len(s_spikes[N[0]]))] 
			if len(b_pass) < 6:	big_mean = [np.mean(bs) for bs in big_spike] 
			MASKS.append((AX,big_mean))

	print len(MASKS) 
	xLen = int(len(MASKS)/2) 
	yLen = 3
	xLoc, yLoc = 0,0 
	
	w = open( spike_file.split('.')[0]+'_spike_types_masks.MASKS', 'w') 

	for AX,m in MASKS:

		ax  = plt.subplot2grid((xLen,yLen), (xLoc,yLoc), rowspan = 1, colspan = 1)
		
		ax.plot(range(len(m)),m,color='k',linewidth=2) 
		ax.set_title(AX) 
		yLoc +=1 
		if yLoc == 3:
			xLoc +=1 
			yLoc = 0 
		w.write('%s %s\n' % (AX, ",".join([str(s) for s in m])))


	plt.savefig(spike_file.split('.')[0]+'_spike_types_masks.png') 
	if options.show:	plt.show() 	



	sys.exit() 




		


				






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	parser.add_option( "--show", default = False,action='store_true' , help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	scan_key(args[0],options)













