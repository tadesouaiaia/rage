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

import random
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
#from ..Rage_Transforms import rage_KDE

sns.set(color_codes=True)
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
	def __init__(self,samples,options,xLen,yLen,key={}):


		self.samples, self.options = samples, options 

		
		if 'r_key' in key: self.r_key = key['r_key']
		if 'p_key' in key: self.p_key = key['p_key']
	

		self.fig = matplotlib.pyplot.gcf()
		self.fig.set_size_inches(19.5, 10.5)
		matplotlib.rcParams['ytick.labelsize'] = 7.5

		self.fig = matplotlib.pyplot.gcf()
		seaborn.set(rc={'axes.facecolor': options.plotColors['AXIS'],  'figure.facecolor':options.plotColors['FACE'], 'savefig.facecolor':options.plotColors['FACE']})
		self.fig.patch.set_facecolor(options.plotColors['FACE'])

		self.xLen, self.yLen = xLen,yLen*3
		self.xLoc, self.yLoc = 0,0 
		self.options = options

		if self.options and self.options.prefix != None: self.prefix = options.prefix
		else:						 self.prefix = 'model_plot'
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 3)
		self.pv_key = dd(int) 



	def update(self,key={}):
		if 'yadd' in key: yAdd = key['yadd']
		else:		  yAdd = 3

		if 'colspan' in key: colspan = key['colspan']
		else:		     colspan = 3 
		if 'clear_axes' in key: 
			self.ax.set_xticks([]) 
			self.ax.set_yticks([]) 
		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])
		if 'title' in key: 	
			self.ax.set_title("\_".join(key['title'].split("_")),fontweight='bold',y=0.99,x=0.0,horizontalalignment='left',fontsize=12)

		if 'xlim' in key: self.ax.set_xlim(key['xlim'][0],key['xlim'][1])


		self.yLoc += yAdd 
		if self.yLoc == self.yLen: 
			self.xLoc,self.yLoc = self.xLoc+1, 0 
		if self.xLoc >= self.xLen: return False
		
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = colspan)
		return True 







	def add_pca_pts(self,pca,key={}):
		if 'comps' in key:	comps = key['comps']
		else: 			comps = [0,1]
		if 'colspan' in key: self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = key['colspan'])

		items,labels,OBS,parents,children = [],[],dd(bool),[],[]
		
		my_data = [(i,p,self.samples[i]) for i,p in enumerate(pca['pts'])]
		
		for i,p,s in my_data:
			
			if self.samples[i].label != None:
				self.ax.plot(p[0],p[1],alpha=0.66,**self.samples[i].label.vals)

				if not OBS[self.samples[i].label.id]:
					sLab = self.samples[i].label.id
					sParents  = [sl.split('~')[0].split('=')[0] for sl in sLab.split('|')]
					sChildren = [sl.split('~')[-1].split('=')[-1] for sl in sLab.split('|')]
					sChildren = [",".join(sc.split(',')[0:2])+'...' if len(sc.split(','))>2 else sc for sc in sChildren ]
				
					parents.append(' & '.join(sParents)) 	
					labels.append('\n'.join(sChildren))
#					labels.append(S[i].label.id)
					items.append(self.samples[i].label.proxy)
					OBS[self.samples[i].label.id] = True

			else:
				self.ax.plot(p[0],p[1])

		self.ax.set_xlabel(pca['axes'][0]) 
		self.ax.set_ylabel(pca['axes'][1]) 

		self.ax.set_xticks([])
		self.ax.set_yticks([]) 	
		if len(items)>0 and self.yLoc ==2:
			ncol = len(labels)
			SL = len("".join(labels))
			yMin,yMax =self.ax.get_xlim() 
			xMin,xMax = self.ax.get_ylim() 
			xLen, yLen = xMax-xMin, yMax-yMin 
			title = list(set(parents))[0]
			if ncol < 12:
				if SL < 100: leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.15,  fontsize=15.5,ncol=ncol,bbox_to_anchor=(0.5,1.75),loc='upper center')
				else: 	     leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.175,  fontsize=10.5,ncol=ncol,bbox_to_anchor=(0.5,1.75),loc='upper center')
			elif ncol < 35:
				ncol = len(labels)/2
				if SL < 300: leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.14,  fontsize=10.5,ncol=ncol,bbox_to_anchor=(0.5,1.9),loc='upper center')
				else: 	     leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.14,  fontsize=9.0,ncol=ncol,bbox_to_anchor=(0.5,1.9),loc='upper center')
			elif ncol < 60:
				ncol = len(labels)/3
				if SL < 500: leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.125,  fontsize=9.0,ncol=ncol,bbox_to_anchor=(0.5,1.95),loc='upper center')
				else: 	     leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.125,  fontsize=8.5,ncol=ncol,bbox_to_anchor=(0.5,1.95),loc='upper center')
			else:
				ncol = len(labels)/4
				if SL < 500: leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.15,  fontsize=15.5,ncol=ncol,bbox_to_anchor=(0.5,1.75),loc='upper center')
				else: 	     leg = self.ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.15,  fontsize=10.5,ncol=ncol,bbox_to_anchor=(0.5,1.75),loc='upper center')

			if 'labels' in self.samples.notes.keys(): 
				title =  " & ".join([kind+'='+val for kind,val in self.samples.notes['labels'].items() if val != None]) 

			leg.set_title(title,prop={'size':12,'weight': 'bold'}) 


		return self		


		


	def change_limits(self,key):
		x0,x1 = self.ax.get_xlim()
		y0,y1 = self.ax.get_ylim()
		if 'x0' in key: x0 = key['x0'] 
		if 'x1' in key: x1 = key['x1'] 
		if 'y0' in key: y0 = key['y0'] 
		if 'y1' in key: y1 = key['y1'] 
		self.ax.set_xlim(x0,x1)
		self.ax.set_ylim(y0,y1)





	def add_labels1(self,title,xlabel,ylabel):
		self.ax.set_title(title) 
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

	def add_legend1(self,labels,colors):
		labs,items = [],[] 
		for a,b in zip(labels,colors):
			labs.append(a) 
			items.append(Rect((0,0),1,1,fc=b))

		plt.legend(items,labs,loc = 'upper left',ncol = len(labs)/3,bbox_to_anchor=(-0.1,1.16),fontsize=10)
		




	def add_pv_bars(self,pvs,options,sim_pvs=[]):


		r_width = 0.2 
		s_width = (1.0-r_width)/len(sim_pvs)
		myticks,myticklabels = [] , []
		for i,k in enumerate(self.p_key):

			myticks.append(i) 
			myticklabels.append(' < '+str(k))
			self.ax.bar(i,pvs[i],width=r_width,color=options.plotColors['RESULTS'])
			iLoc =i+r_width
			for b in sorted([bh[i] for bh in sim_pvs],reverse=True):
				self.ax.bar(iLoc,b,width=s_width,color=options.plotColors['SIMS']) 
				iLoc+=s_width
		self.ax.set_xticks(myticks) 	
		self.ax.set_xticklabels(myticklabels,horizontalalignment='left')		
		items = [Rect((0,0),1,1,fc=options.plotColors['RESULTS']),Rect((0,0),1,1,fc=options.plotColors['SIMS'])]
		self.ax.legend(items,['Results','Simulation'],loc = 'upper right',bbox_to_anchor=(1,1),fontsize=10)
		return self


	def add_rs_bars(self,rs,options,sim_rs=[]):

		r_width = 0.2 
		s_width = (1.0-r_width)/len(sim_rs)
		myticks,myticklabels = [] , []
		for i,k in enumerate(self.r_key):

			myticks.append(i) 
			myticklabels.append(' > '+str(k))
			self.ax.bar(i,rs[i],width=r_width,color=options.plotColors['RESULTS'])
			iLoc =i+r_width
			for b in sorted([bh[i] for bh in sim_rs],reverse=True):
				self.ax.bar(iLoc,b,width=s_width,color=options.plotColors['SIMS']) 
				iLoc+=s_width
		self.ax.set_xticks(myticks) 	
		self.ax.set_xticklabels(myticklabels,horizontalalignment='left')		
		items = [Rect((0,0),1,1,fc=options.plotColors['RESULTS']),Rect((0,0),1,1,fc=options.plotColors['SIMS'])]
		self.ax.legend(items,['Results','Simulation'],loc = 'upper right',bbox_to_anchor=(1,1),fontsize=10)
		return self



	def add_model_table(self,M,X,skrees=None,key={}):


		cell_data, row_names, row_keys = [], ['BIC','$R^2$','$R^2_{adj}$','Power-05','Power-001'], ['bic','rsq','rsa','pwr1','pwr2']

		for k in row_keys:

			cell_data.append([round(np.mean(M.out[k]),3), round(np.std(M.out[k]),3),round(np.percentile(M.out[k],90),3), round(np.percentile(M.out[k],99),3)])
	
		col_names = ['$ \overline{X}$','$S^2$','top 10%','top 1%']
		colors = plt.cm.BuPu(np.linspace(0, 0.5, len(row_names)))
		the_table = self.ax.table(cellText=cell_data,
			rowLabels=row_names,
			rowColours=colors,
			colLabels=col_names,
			colWidths = [0.2,0.2,0.2,0.2,0.2],
			loc='bottom left',
			bbox=(0.11,0.1,0.85,0.95))
#        	the_table.auto_set_font_size(False); the_table.set_fontsize(13.5); the_table.scale(1,2.0)

		x_scale,t2_fs,t2_sc,t2_cw = 1, 12, 1.5,  0.12
		if skrees != None: 
			init                      = skrees[0]
			covars, model, sims      = int(100*skrees[1]/init), int(100*skrees[2]/init), int(100*np.mean(skrees[3])/init)
			col_names,row_names = ['init','covariates','model','simulation'], ['Total Variance']
			cell_data = [[100,covars,model,sims]]
			# print skrees 
			colors = ['#d8daeb']

			table2 = self.ax.table(cellText=cell_data,rowLabels=row_names,rowColours=colors,colLabels=col_names,colWidths = [t2_cw,t2_cw,t2_cw,t2_cw],loc='bottom left',bbox=(0.13,-0.65,0.75,0.30))
			table2.auto_set_font_size(False)
			table2.set_fontsize(t2_fs)

			table2.scale(1,t2_sc) 


		self.ax.set_title('MODEL SUMMARY',fontweight='bold',fontsize=12,y=1.1) 
		self.ax.axis('off') 
		return self





	def add_predictor_table(self,M,Mp,Mc,X,options,key={}):
		
		nD,colors,row_names,cell_data = 0,[],[],[]
		col_names = ['VIF','$ \overline{pv}$','top 10%','top 1%']+['$ \overline{pv}$ (NM)','top 10%','top 1%']
		col_colors = ['white' for c in col_names]
		key_top,key_bottom= dd(list),dd(list) 		

		for n,P in M.pv_dict.items():
			if n == 'intercept': continue 
			m_stats = [M.vif[n]] + [round(np.mean(P),3),  '%2.1e' % np.percentile(P,10), '%2.1e' % np.percentile(P,1)]

			if X.PREDICTOR[n]:key_top[n] =  m_stats + [round(np.mean(Mp.pv_dict[n]),3),  '%2.1e' % np.percentile(Mp.pv_dict[n],10), '%2.1e' % np.percentile(Mp.pv_dict[n],1)]
			else: 		  key_top[n] =  m_stats + [round(np.mean(Mc.pv_dict[n]),3),  '%2.1e' % np.percentile(Mc.pv_dict[n],10), '%2.1e' % np.percentile(Mc.pv_dict[n],1)]

		row_names = key_top.keys()+key_bottom.keys() 
		cell_data = key_top.values()+key_bottom.values() 
		colors =  [options.plotColors['FOUND'] for x in key_top] + [options.plotColors['FALSEPOS'] for x in key_bottom]
		x_scale,t2_fs,t2_sc,t2_cw = 1, 12, 1.5,  0.12

		if len(row_names) < 4: 		colsteps,y_scale,fs = 0.11, 2 , 11
		elif len(row_names) < 6: 	colsteps,y_scale,fs = 0.10, 1.75 , 10
		elif len(row_names) < 8:	colsteps,y_scale,fs = 0.09, 1.5, 9
		elif len(row_names) < 12:
			colsteps,y_scale,fs = 0.075, 1.24, 9.5 
			t2_fs,t2_sc,t2_cw = 11, 1.3,  0.11
		elif len(row_names) < 15:
			colsteps,x_scale,y_scale,fs = 0.06, 1.5, 1.0, 9.0 
			t2_fs,t2_sc,t2_cw = 10, 1.2,  0.10
		else:
			colsteps,x_scale,y_scale,fs = 0.04, 2.85, 0.90, 6.5 
			t2_fs,t2_sc,t2_cw = 9, 1.1,  0.08

		maxRow = max([len(x) for x in row_names])
		if maxRow < 10: bX = 0.1
		elif maxRow < 15: bX = 0.2
		elif maxRow < 20: bX = 0.275
		else: 		  bX = 0.3


		cW = [colsteps for i in range(len(row_names)+10)]
		the_table = self.ax.table(cellText=cell_data,rowLabels=row_names,rowColours=colors,colLabels=col_names,colColours=col_colors,colWidths =cW,loc='center')
		the_table.auto_set_font_size(False)
		the_table.set_fontsize(fs)
		the_table.scale(x_scale,y_scale) 



		if 'sim_pvs' in key:

			cells = [M.pv_cnt,[round(np.mean([n[j] for n in key['sim_pvs']]),2) for j in range(len(M.pv_cnt))]]
			col_names, row_names = M.pv_key, ['Found','FalsePos']
			col_widths = [0.10, 0.10,0.10,0.10,0.10,0.10,0.10,0.10, 0.10, 0.10]
			colors =[options.plotColors['FOUND'], options.plotColors['FALSEPOS']]
			table2 = self.ax.table(cellText=cells,rowLabels=row_names,rowColours=colors,colLabels=col_names,colWidths = col_widths,loc='bottom left',bbox=(0.1,-0.65,0.825,0.3))
		self.ax.set_title('PREDICTOR CONTRIBUTIONS',fontweight='bold',fontsize=12,y=1.2) 
		self.ax.axis('off') 
		return self






	def save_mfig(self,mt,preds,covars,key={}):


		p_key, f_key, m_str, x = dd(int), {},  mt+':  y  \propto  ', 1
		for p in preds+covars:
			if p == 'intercept': m_str+=' B_o +'
			else:
				p_key[p.split('=')[0]]+=1
				f_key[p.split('=')[0]] = p 

		for i,(n,p) in enumerate(sorted([(b,a) for (a,b) in p_key.items()])):
			if n == 1: p = f_key[p] 
			if len(p.split('_')) > 1: p = '\_'.join(p.split('_')) 
			if i == 0:
				if n == 1:  m_str+= ' B_'+str(x)+'('+p+')'
				else:       m_str+= ' B_'+str(x)+'('+p+')_'+str(n)
			else:
				if n == 1:  m_str+= ' + B_'+str(x)+'('+p+')'
				else:       m_str+= ' + B_'+str(x)+'('+p+')_'+str(n)
			x+=1


		plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.90,wspace=0.25,hspace=1.05)
		plt.suptitle('$'+m_str+'$',fontsize=25, fontweight='bold')
		if 'axis_off' in key: plt.axis('off') 

		fig_name = self.options.prefix+"_"+mt+"_predictorModelAnalysis_"+"_".join(self.options.predictors)+"_covars_"+"_".join(covars)+".pdf" 
		self.fig.savefig(str(fig_name),dpi=200)
	
	
		if self.options.show: 
			plt.show() 




class eplot:
	def __init__(self,options,pv_discs,key={}):

#		if 'r_key' in key: self.r_key = key['r_key']
#		if 'p_key' in key: self.p_key = key['p_key']

		self.pv_discs = sorted(pv_discs)
		self.dK = []
		mycolor = 'cyan'	
		sns.set(rc={'axes.facecolor':'skyblue', 'figure.facecolor':'cornflowerblue'})
		self.fig = matplotlib.pyplot.gcf()
		self.fig.set_size_inches(19.5, 10.5)
		self.fig.set_facecolor(mycolor) 
		self.fig.patch.set_facecolor(mycolor)
		matplotlib.rcParams['savefig.facecolor'] = mycolor


		self.xLen, self.yLen = 11,1
		self.xLoc, self.yLoc = 0,0 
		self.options = options

		if self.options and self.options.prefix != None: self.prefix = options.prefix
		else:						 self.prefix = 'covar_plot'
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)




	def add_base_model(self,BM,Xp,BS):

		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)

		col_names = ['$ \overline{R^2}$','$P_{99}(R^2)$','$ \overline{pv}$','$P_{99}(pv)$','$pv<'+str(round(BM.pv_key[0],4))+'$','$pv<'+str(round(BM.pv_key[-1],8))+'$']
		row_names = ['$Model_{base}$','$Model_{sim}$']

		t1 = [round(BM.stats['rs'].mean,4), round(BM.stats['rs'].p1,3),round(BM.stats['pv'].mean,3),'%2.2e' % BM.stats['pv'].p1,BM.pv_cnt[0],BM.pv_cnt[-1]]
		
		t2 = [round(np.mean([s.stats['rs'].mean for s in BS]),3),round(np.mean([s.stats['rs'].p1 for s in BS]),3),round(np.mean([s.stats['pv'].mean for s in BS]),3),'%2.2e' % np.mean([s.stats['pv'].p1 for s in BS])]
		t2.extend([round(np.mean([s.pv_cnt[0] for s in BS]),3),round(np.mean([s.pv_cnt[-1] for s in BS]),3)])

		table_data = [t1,t2]

		col_widths = [0.08, 0.08,0.08,0.08,0.08,0.08,0.08,0.09, 0.08, 0.08]
		colors = plt.cm.hsv(np.linspace(0, 0.5, len(row_names)))[-1::-1]
		colors = ['red','white','yellow','purple']

		base_table = self.ax.table(cellText=table_data,rowLabels=row_names,colLabels=col_names,colWidths = col_widths,loc='center')
		base_table.auto_set_font_size(False)
		base_table.set_fontsize(11.0)
		base_table.scale(1,1.5)
		self.ax.axis('off') 
		self.xLoc+=1

		title_string = ' $Base Model :  \propto B_o +'+"+".join(['B_'+str(ni+1)+'*'+"\_".join(mn.split('~')[-1].split('_')) for ni,mn in enumerate(Xp) if mn != 'intercept'])+'$'

		self.ax.set_title(title_string)

		if len(title_string) < 50: self.ax.set_title(title_string,fontsize = 35, y = 1.1)
		elif len(title_string) < 100: self.ax.set_title(title_string,fontsize = 15, y = 1.2)
		else: 
			title_string = title_string[0:150]+'...'
			self.ax.set_title(title_string,fontsize = 13, y = 1.25)
		
		return self


	def add_covariate_data(self,covar_data):

		type_dict = {'binary': 'categorical', 'continuous': 'continuous'} 
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 10, colspan = 1)
		row_names,cell_data =[],[]
		col_names = ['type','size','vif','$\overline{R^{2}_{c}}$','$\overline{R^{2}_{c+m}}$']
			
		for cName,cData in covar_data.items():
			row_names.append(cName)
			c_stats,m_disc,c_disc = cData
			if len(self.dK) == 0:
				self.dK = sorted(m_disc.keys(),reverse=True)
				col_names +=['$M_{net}  < %2.2f$' % self.dK[0],'$M_{net} $'+'<%1.1e '%  self.dK[1]]
				col_names +=['$S_{net}  < %2.2f$' % self.dK[0],'$S_{net} $'+'<%1.1e '%  self.dK[1]]

			cell_data.append(c_stats + [m_disc[self.dK[0]],m_disc[self.dK[1]], c_disc[self.dK[0]],c_disc[self.dK[1]]])

		cell_colors = [['white' for x in col_names] for y in range(len(row_names))]
		
		for j,name in enumerate(col_names):
			c_srt = sorted([(cell_data[i][j],i,j) for i in range(len(cell_data))])
			if j in [1,5,6]:
				cell_colors[c_srt[0][1]][c_srt[0][2]] = 'tomato'
				if c_srt[-1][0] > 0: 
					cell_colors[c_srt[-1][1]][c_srt[-1][2]] = 'lightgreen'
			if j in [7,8]:
				if c_srt[0][0] != 0:
					cell_colors[c_srt[0][1]][c_srt[0][2]] = 'lightgreen'
				if c_srt[-1][0] > 0: 
					cell_colors[c_srt[-1][1]][c_srt[-1][2]] = 'tomato'
			
		if len(cell_data) < 6: 
			sw,mw, fs,xc, sc = 0.07, 0.1, 13.5,1,  3.0

		elif len(cell_data) < 12:
			sw, mw, fs,xc, sc = 0.055, 0.08, 12.5,1,  2.5
		elif len(cell_data) < 16:
			sw, mw, fs,xc, sc = 0.04, 0.07, 12.0,1,  2.2
		elif len(cell_data) < 22:	
			sw, mw, fs,xc, sc = 0.04, 0.075, 12.0, 1,  1.6
		elif len(cell_data) < 30:	
			sw, mw, fs,xc, sc = 0.05, 0.075, 12.0,1,  1.5
		elif len(cell_data) < 40:	
			sw, mw, fs,xc, sc = 0.05, 0.075, 12.25,1,  1.0
		elif len(cell_data) < 50:	
			sw, mw, fs,xc, sc = 0.044, 0.068, 11.0, 1,  0.89
		elif len(cell_data) < 60:	
			sw, mw, fs,xc, sc = 0.033, 0.060, 9.0,1,  0.77
		elif len(cell_data) < 75:	
			sw, mw, fs,xc, sc = 0.030, 0.060, 7.8,1,  0.66
		else:
			sw, mw, fs,xc, sc = 0.0275, 0.052, 6.5, 1.25,  0.51
		col_widths = [mw,sw,sw,mw,mw,mw,mw,mw,mw,mw,mw,mw]



		base_table = self.ax.table(cellText=cell_data,rowLabels=row_names,colLabels=col_names,cellColours = cell_colors,colWidths = col_widths,loc='center')
		base_table.auto_set_font_size(False)
		base_table.set_fontsize(fs)
		base_table.scale(xc,sc)
		self.ax.axis('off') 

		return self	












	def save_efig(self,mt,preds,covars,key={}):
		plt.subplots_adjust(left=0.02, bottom=0.01,right=0.99)
		if 'axis_off' in key: plt.axis('off') 
		
		self.fig.savefig(self.options.prefix+'_'+mt+'_covariateAnalysis_'+'_'.join(self.options.predictors)+'_covar_'+'|'.join(covars)+'.pdf', dpi=200) 
	
		if self.options.show: 
			plt.show() 







































class predictor_plot:
	def __init__(self,options,xLen,key={}):
		self.options = options 
	
		sns.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'cornflowerblue'})
		self.fig = matplotlib.pyplot.gcf()
		self.fig.set_size_inches(19.5, 10.5)

		matplotlib.rcParams['ytick.labelsize'] = 5.5
		matplotlib.rcParams['xtick.labelsize'] = 3.5
		matplotlib.rcParams['savefig.facecolor'] = options.plotColors['FACE']

		seaborn.set(rc={'axes.facecolor': options.plotColors['AXIS'],  'figure.facecolor':options.plotColors['FACE']})
		self.fig.patch.set_facecolor(options.plotColors['FACE'])
		self.fig.set_facecolor(options.plotColors['FACE']) 

		self.xLen, self.yLen = xLen*2, 5
		self.xLoc, self.yLoc = 0,0 
		self.options = options

		if self.options and self.options.prefix != None: self.prefix = options.prefix
		else:						 self.prefix = 'model_plot'
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 3)
		self.pv_key = dd(int) 



	def add_predictor_row(self,predictor,samples,pcaI,pcaR,Mi,sim_res): 

		self.predictor = predictor

		ax1 = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc),rowspan=2,colspan=1)
		ax2 = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc+1),rowspan=2,colspan=1)
		ax3 = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc+2),rowspan=1,colspan=1)
		ax4 = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc+1,self.yLoc+2),rowspan=1,colspan=1)
		ax5 = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc+3),rowspan=2,colspan=2)


		var_red = str(int(100*pcaR['total_var'] / pcaI['total_var']))+'%'
		null_red = str(int(100* np.mean(sim_res['var']) / pcaI['total_var']))+'%' 

		self.add_pca_data(ax1,samples,pcaI,LEGEND=True,TITLE='Raw Data ( '+str(len(samples))+' )') 
		self.add_pca_data(ax2,samples,pcaR,TITLE='Residualized ( remaining variance: '+var_red+' null='+null_red+' )') 

#		print Mi.pv_cnt, sim_res['pv'] 

		self.add_pv_bars(ax3,Mi.pv_key,Mi.pv_cnt,sim_res['pv'])		
		self.add_rs_bars(ax4,Mi.rs_key,Mi.rs_cnt,sim_res['rs'])		

		self.add_summary_table(ax5,Mi,samples) 


		self.xLoc +=2





	def add_pca_data(self,ax,samples,pca,LEGEND=False,TITLE=None):
		items,labels,OBS,parents,children = [],[],dd(bool),[],[]
		for i,p in enumerate(pca['pts']): 

			if samples[i].label != None: 
				ax.plot(p[0],p[1],alpha=0.66,**samples[i].label.vals)
			if not OBS[samples[i].label.id]:
				sLab = samples[i].label.id
				sParents  = [sl.split('~')[0].split('=')[0] for sl in sLab.split('|')]
				sChildren = [sl.split('~')[-1].split('=')[-1] for sl in sLab.split('|')]
				sChildren = [",".join(sc.split(',')[0:2])+'...' if len(sc.split(','))>2 else sc for sc in sChildren ]		
				parents.append(' & '.join(sParents)) 	
				labels.append('\n'.join(sChildren))
				items.append(samples[i].label.proxy)
				OBS[samples[i].label.id] = True
			else:
				ax.plot(p[0],p[1])

			ax.text(p[0],p[1],samples[i].name)

		ax.set_xlabel(pca['axes'][0]) 
		ax.set_xticks([])
		ax.set_yticks([]) 
		if TITLE: ax.set_title(TITLE,fontsize=10) 
		if LEGEND:
			SL = len("".join(labels))
			if   len(labels) < 6:
				leg = ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.15,  fontsize=9,ncol=1,bbox_to_anchor=(-0.30,0.95),loc='upper center')
			elif len(labels) < 15:  
				leg = ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.12,  fontsize=8,ncol=2,bbox_to_anchor=(-0.60,0.95),loc='upper center')
			else:
				leg = ax.legend(items,labels,handletextpad=-0.1,columnspacing=0.12,  fontsize=7,ncol=3,bbox_to_anchor=(-0.80,0.95),loc='upper center')

			leg.set_title(self.predictor,prop={'size':10,'weight': 'bold'}) 

		return self		



	def add_pv_bars(self,ax,p_key,pvs,sim_pvs=[]):


		r_width = 0.2 
		s_width = (1.0-r_width)/len(sim_pvs)
		myticks,myticklabels = [] , []
		for i,k in enumerate(p_key):

			myticks.append(i) 
			myticklabels.append(' < '+str(k))
			ax.bar(i,pvs[i],width=r_width,color='blue')
			iLoc =i+r_width
			for b in sorted([bh[i] for bh in sim_pvs],reverse=True):
				ax.bar(iLoc,b,width=s_width,color='red') 
				iLoc+=s_width
		ax.set_xticks(myticks) 	
		ax.set_xticklabels(myticklabels,horizontalalignment='left')		
		ax.set_title(self.predictor,fontsize=18,fontweight='bold',y=1.1) 
		ax.set_xlabel('Nominal P-value') 
		items = [Rect((0,0),1,1,fc='b'),Rect((0,0),1,1,fc='r')]
		ax.legend(items,['Results','Simulation'],loc = 'upper right',bbox_to_anchor=(1,1),fontsize=10)
		return self

		

	def add_rs_bars(self,ax,p_key,pvs,sim_pvs=[]):


		r_width = 0.2 
		s_width = (1.0-r_width)/len(sim_pvs)
		myticks,myticklabels = [] , []
		for i,k in enumerate(p_key):

			myticks.append(i) 
			myticklabels.append(' >= '+str(k))
			ax.bar(i,pvs[i],width=r_width,color='blue')
			iLoc =i+r_width
			for b in sorted([bh[i] for bh in sim_pvs],reverse=True):
				ax.bar(iLoc,b,width=s_width,color='red') 
				iLoc+=s_width
		ax.set_xticks(myticks) 	
		ax.set_xticklabels(myticklabels,horizontalalignment='left')	
		ax.set_title('Variance Explained') 	
		items = [Rect((0,0),1,1,fc='b'),Rect((0,0),1,1,fc='r')]
		ax.legend(items,['Results','Simulation'],loc = 'upper right',bbox_to_anchor=(1,1),fontsize=10)
		return self

		




		









	def add_summary_table(self,ax,M,S): 

		
			
		cell_data, row_names, row_keys = [], ['BIC','$R^2$','$R^2_{adj}$','Power-001'], ['bic','rsq','rsa','pwr2']
		for k in row_keys:	
			cell_data.append(['-',round(np.mean(M.out[k]),3), round(np.percentile(M.out[k],99),2), round(max(M.out[k]),2)])

		col_names = ['n','$ \overline{X}$','top 1%','top']

		for k,vals in M.pv_dict.items():
			if k == 'intercept': continue 

			if S.attribute_class[k.split('~')[0]] == 'binary': 
				sLen = len([s for s in  S if s.attributes[k.split('~')[0]] == k])
			
			else:   sLen = len(S) 
		


			row_names.append(k) 
			cell_data.append([sLen,round(np.mean(vals),1), '%1.2e ' % np.percentile(vals,1), '%1.2e ' % min(vals)])
			



			colors = plt.cm.BuPu(np.linspace(0, 0.5, len(row_names)))


			the_table = ax.table(cellText=cell_data,rowLabels=row_names,rowColours=colors,colLabels=col_names,colWidths = [0.1,0.11,0.15,0.13],bbox=(0.25,0.01,0.80,0.90))
			the_table.auto_set_font_size(False)
		if len(row_names) < 8: 				the_table.set_fontsize(11)
		elif len(row_names) < 12:		       the_table.set_fontsize(9) 
		else:		       				the_table.set_fontsize(8) 
		ax.axis('off') 

	def save(self,outname):

		plt.subplots_adjust(left=0.12, bottom=0.03, right=0.90, top=0.85,wspace=0.5,hspace=0.75)
		plt.suptitle('$'+self.options.model+':  y  \propto  (predictor) + covariates: '+",".join(self.options.covariates)+'$',fontsize=18,fontweight='bold')
		plt.savefig(outname+'.png',dpi=200) 
		if self.options.show: plt.show()



