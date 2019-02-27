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
from ..Rage_Transforms import rage_KDE
sns.set(color_codes=True)



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
        def __init__(self,xLen,yLen,options=None):

	
		sns.set(rc={'axes.facecolor':'pink', 'figure.facecolor':'cornflowerblue'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(18.5, 9.5)
		self.xLen, self.yLen = xLen,yLen
		self.xLoc, self.yLoc = 0,0 
		self.options = options

		if self.options and self.options.prefix != None: self.prefix = options.prefix
		else:						 self.prefix = 'model_plot'


		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		self.pv_key = dd(int) 

	def update(self,key={}):

		if 'clear_axes' in key: 
			self.ax.set_xticks([]) 
			self.ax.set_yticks([]) 

		


		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])
		if 'title' in key: 	self.ax.set_title(key['title'],fontweight='bold',y=1.05,x=0.0,horizontalalignment='left',fontsize=14)

		if 'xlim' in key: self.ax.set_xlim(key['xlim'][0],key['xlim'][1])


		self.yLoc += 1 
		if self.yLoc == self.yLen: 
			self.xLoc,self.yLoc = self.xLoc+1, 0 
		if self.xLoc >= self.xLen: return False
		
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		return True 



	def save(self,f_name,key={}):

		if 'axis_off' in key: plt.axis('off') 
		if 'title' in key: plt.suptitle(key['title'],fontsize=20,fontweight='bold')
		elif self.options and 'title' in vars(self.options).keys() and self.options.title != None:	plt.suptitle(self.options.title)
               	self.fig.savefig(f_name, dpi=100) 

		if self.options.show: 
			plt.show() 



	def add_hist(self,h_data):


		if type(h_data) == list:
			x,y = range(len(h_data)), sorted(h_data)
		else: 

			x,y = h_data.keys(), sorted(h_data.values())




		if len(self.ax.hist(y,bins='auto')[0]) < 2:
			self.ax.clear() 
			self.ax.hist(np.hstack(y),bins=min(4,len(cc(y)))) 
		yMin, yLim = self.ax.get_ylim()
		xMin, xLim = self.ax.get_xlim()
		HI,LO=False,False

		

		if yMin == 0: 

			h_srt = sorted(y) 
			out_bool, out_scrs, out_colors =  mad_based_outlier(np.array(h_srt)) 
			g_color = out_colors[sorted([(out_scrs[oi],oi) for oi in range(len(out_scrs))])[0][1]]

			
			if len(h_srt) < 100: s = 10 
			elif len(h_srt) < 1000: s = 5
			else: 			s = 3

			for i,h in enumerate(h_srt): 
				if not out_bool[i]:
					self.ax.scatter(h,yLim*1.025,alpha=0.7,color=g_color,s=s,clip_on=False) 
				else: 
					self.ax.scatter(h,yLim*1.025,alpha=0.7,color=out_colors[i],s=s,clip_on=False) 
					if i > 5 and (out_scrs[i]/(out_scrs[i-1]+out_scrs[i-2])) > 0.66:  HI=True 
					if HI: 			
						try: self.ax.text(h,yLim*1.030,x[i].name.split(";")[-1])
						except AttributeError: continue 

			self.ax.set_ylim(yMin,yLim) 
			self.ax.set_xlim(xMin,xLim) 
		return self




	def add_pca_data(self,pca_pts,key={}):

		sns.set(rc={'axes.facecolor':'black', 'figure.facecolor':'cornflowerblue'})
		if 'components' in key: comps = key['components']
		else:			comps = [0,1] 
		if 'type' in key:
			if key['type'] == 'tsne':	
				axis_labs,title = ['TS'+str(c+1) for c in comps], 'TSNE'
			elif key['type'] == 'kca': 
				axis_labs,title = ['KC'+str(c+1) for c in comps], 'KCA'
			elif key['type'] == 'mds': 
				axis_labs,title = ['MD'+str(c+1) for c in comps], 'MDS'
			else:
				axis_labs,title = ['X'+str(c+1) for c in comps], key['type'].upper() 
		else:
			title = 'PCA'
			if 'vars' in key: 				axis_labs = ['PC'+str(c+1)+':  '+str(round(100*key['vars'][c],2))+'%' for c in comps]
			else:						axis_labs = ['PC'+str(c+1) for c in comps]

		for i,p in enumerate(pca_pts):
			if 'colors' in key and key['colors'][i] != 'white': 
				if 'sizes' in key: self.ax.scatter(p[0],p[1],s = key['sizes'][i],color=key['colors'][i],alpha=0.6)
				else: 		   self.ax.scatter(p[0],p[1],color=key['colors'][i])
			else: 		    self.ax.scatter(p[0],p[1])

		self.ax.set_xlabel(axis_labs[0]) 	
		self.ax.set_ylabel(axis_labs[1])

		if 'zoom' in key: 
			pX = sorted([p[0] for p in pca_pts])
			pY = sorted([p[1] for p in pca_pts])
			xQ1,xQ3 = int(len(pX)*0.25), int(len(pX)*0.75)
			yQ1,yQ3 = int(len(pY)*0.25), int(len(pY)*0.75)
			iqX = pX[xQ3]-pX[xQ1]
			iqY = pY[yQ3]-pY[yQ1]  
			self.ax.set_xlim(pX[xQ1]-(iqX*2.5),pX[xQ3]+(iqX*2.5))
			self.ax.set_ylim(pY[yQ1]-(iqY*2.5),pY[yQ3]+(iqY*2.5))



		if 'title' in key: self.ax.set_title(key['title']) 
		else:		   self.ax.set_title(title) 
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


	def add_lines(self,x,y,xlabel=None,ylabel=None,clr=None):
		if not clr: self.ax.plot(x,y,zorder=1)
		else: 	    self.ax.plot(x,y,color=clr,zorder=1) 
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 
		return self

	def add_line(self,x,y,key={}):
		lw,alpha = 1.0, 1.0 
		if 'lw' in key.keys():    lw = key['lw'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'color' in key: self.ax.plot(x,y,linewidth=lw,color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.plot(x,y,linewidth=lw,zorder=1,alpha=alpha)
		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])


	def scatter_pts(self,X,Y,key={}):
		alpha,size,mark =  1.0, 20,'o'
		if 'size' in key.keys():    size = key['size'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'mark' in key: mark= key['mark'] 
		if 'yjitter' in key: Y = [np.random.normal(y,(0.0195)) for y in Y]
		

		if 'color' in key: self.ax.scatter(X,Y,s=size,marker = mark, color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.scatter(X,Y,zorder=1,alpha=alpha,s=size)


	def add_scatter_data(self,x,y,xlabel=None,ylabel=None,clr=None):
		if clr and len(clr.split(',')) == 2: 
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),scatter_kws={"color": clr.split(',')[0]}, line_kws={"color": clr.split(',')[1]})
		else:
			if clr == None: clr = 'k'
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),color=clr)
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

	
	def add_scatter_pts(self,X,Y,notes='None'):
		notes = notes.split(',') 
		if 'RED' in notes: clr = 'red' 
		elif 'PURPLE' in notes: clr = 'purple' 
		elif 'BLUE' in notes: clr = 'blue'
		else:		      clr = 'k'  
		if 'XY_JITTERS' in notes:
			xJ = [np.random.normal(x,(0.0000001+x*0.1)) for x in X]
			yJ = [np.random.normal(y,(0.0000001+y*0.2)) for y in Y]
			self.ax.scatter(xJ,yJ,color=clr)

		else:
			self.ax.scatter(X,Y,color=clr) 


	




	def add_labels(self,title,xlabel,ylabel):
		self.ax.set_title(title) 
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

	def add_legend(self,labels,colors):
		labs,items = [],[] 
		for a,b in zip(labels,colors):
			labs.append(a) 
			items.append(Rect((0,0),1,1,fc=b))

		plt.legend(items,labs,loc = 'upper left',ncol = len(labs)/3,bbox_to_anchor=(-0.1,1.16),fontsize=10)
		
	def add_skree(self,pca_vars,key={}):
		if 'color' in key: 	self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))],color=key['color'])
		else:			self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))])
		return self	

	def add_skrees(self,pca_var_list,key={}):
		for pca_vars in pca_var_list:


			if 'color' in key: 	self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))],color=key['color'])
			else:			self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))])
		return self	

	def add_skree_legend(self,sum_key,color_key):
		maxV = max(sum_key.values())
		labs,items = [],[] 
		for k in sorted(sum_key.keys()):
			myVal = str(int(100*(sum_key[k] / maxV)))

			labs.append(k+' :'+myVal)
			#items.append(Rect((0,0),1,1,fc=b))
			items.append(Line2D([],[],linestyle='-', color=color_key[k], marker=None))

		self.ax.set_xticks([]) 
		self.ax.set_yticks([]) 
		self.ax.set_xlabel('PCA Components') 
		self.ax.set_ylabel('Total Variance')
#		self.ax.set_xlim(0,min(20,self.ax.get_xlim()[1]))
		self.ax.legend(items,labs,loc = 'upper left',bbox_to_anchor=(0.7,0.4),fontsize=10)
		self.yLoc += 1 
		if self.yLoc == self.yLen: 
			self.xLoc,self.yLoc = self.xLoc+1, 0 
		if self.xLoc >= self.xLen: return False

		self.ax.set_title('Skree Plots') 
		
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		return True 


	

	def add_density(self,pts):
		kde = rage_KDE.samples(0.3)
		x,y = kde.run(pts)
		lw,alpha=1,1
		self.ax.plot(x,y,linewidth=lw,zorder=1,alpha=alpha) 
		#self.ax.plot(x,y,linewidth=lw,color=key['color'],zorder=1,alpha=alpha) 
		plt.show()
		sys.exit() 



	def set_pv_key(self,pvs):

		my_keys = self.pv_key.keys() 		
		if len(my_keys)  == 0:

			pv5 = np.percentile(pvs,50)
			my_keys = [0.05,0.01,0.005]
			if pv5 < 0.001: my_keys = []
			elif pv5<0.01:  my_keys = [0.01,0.005] 
			else:		my_key  = [0.05,0.01,005]

	#		my_keys.extend([0.6,0.4,0.2]) 
			p,pv = pvs[0],0.001
			while True: 
				if p < pv: my_keys.append(pv) 
				else:	   break 
				pv /= 10.0 
				if pv < 0.0000001: break
		self.pv_key = {k: 0 for k in my_keys}


	def add_pv_bar(self,pvs,sim_pvs=[]):

		self.set_pv_key(pvs) 
		bar_cnts = [] 
		pvk = sorted(self.pv_key.keys()) 
		for sp in [pvs]+sim_pvs: 
			sp.sort()
			pcnts,k,p = [0 for i in pvk],0,0  
			while True:
				while k<len(pvk) and sp[p] > pvk[k]: k+=1
				while k< len(pvk) and p < len(sp) and sp[p] <= pvk[k]:
					for j in range(k,len(pvk)): pcnts[j]+=1
					p+=1 
				if k == len(pvk) or p == len(sp): break 
			bar_cnts.append(pcnts) 

		myticks,myticklabels = [] , []

		r_width = 0.2 
		s_width = (1.0-r_width)/len(bar_cnts) 
		for i,k in enumerate(pvk[-1::-1]):
			myticks.append(i) 
			myticklabels.append(' < '+str(k))
			idx = -1-i
			self.ax.bar(i,bar_cnts[0][idx],width=r_width,color='blue')
			iLoc =i+r_width
			bh = sorted([b[idx] for b in bar_cnts[1::]],reverse=True) 
 
			for b in bh: 
				self.ax.bar(iLoc,b,width=s_width,color='red') 
				iLoc+=s_width
		self.ax.set_xticks(myticks) 	
		self.ax.set_xticklabels(myticklabels,horizontalalignment='left')		
		items = [Rect((0,0),1,1,fc='b'),Rect((0,0),1,1,fc='r')]
		self.ax.legend(items,['Results','Simulation'],loc = 'upper right',bbox_to_anchor=(1,1),fontsize=10)
		return self






	def add_pv_bars(self,pvs):	

		if len(self.pv_key) == 0:
			p,pv = pvs[0], 0.05 
			while True: 
				if p < pv: self.pv_key[pv] = 0		
				else:	   break 
				pv /= 10.0 
				if pv < 0.00001: break 


		self.pv_key = {k: 0 for k in self.pv_key.keys()}
		pvk = [k for k in sorted(self.pv_key.keys(),reverse=True)]
		my_pvs = dd(int) 
		for p in pvs:
			for k in pvk:
				if p < k: my_pvs[k] +=1
				else:     break 
	
		myticks,myticklabels = [] , []
		for i,k in enumerate(pvk):
			self.ax.bar(i,my_pvs[k],width=0.8)
			myticks.append(i) 
			myticklabels.append(' < '+str(k))
		self.ax.set_xticks(myticks) 	
		self.ax.set_xticklabels(myticklabels,horizontalalignment='left')		
		
		return self

	def add_model_table(self,model,key={}):


		cell_data, row_names, row_keys = [], ['BIC','$R^2$','$R^2_{adj}$','Power'], ['bic','rs','ars','pwr']
		for k in row_keys: 
			k_data = sorted(model[k],reverse=True) 
			#if k == 'bic':  k_data = sorted(model[k])
			#else:		k_data = sorted(model[k],reverse=True) 
			k_1 = int(len(k_data) / 100.0+0.99)
			k1,k10 =k_data[0:int(k_1)],k_data[0:int(k_1*10)]


			cell_data.append([round(np.mean(k_data),3),round(np.std(k_data),2),round(np.mean(k10),3),round(np.mean(k1),3)])

		

 
		col_names = ['$ \overline{X}$','$S^2$','top 10%','top 1%']
        	colors = plt.cm.BuPu(np.linspace(0, 0.5, len(row_names)))
        	the_table = self.ax.table(cellText=cell_data,rowLabels=row_names,rowColours=colors,colLabels=col_names,colWidths = [0.2, 0.2,0.2,0.2],loc='center')
        	the_table.auto_set_font_size(False); the_table.set_fontsize(15.5); the_table.scale(1,2)

		if 'ss' in key:
			col_names,row_names = sorted(key['ss'].keys(),reverse=True), ['Sensitivity','Specificity'] 
			cell_data = [[round(key['ss'][r]['se'],3) for r in col_names],[round(key['ss'][r]['sp'],3) for r in col_names]]
			colors = plt.cm.hsv(np.linspace(0, 0.5, len(row_names)))
			table2 = self.ax.table(cellText=cell_data,rowLabels=row_names,rowColours=colors,colLabels=col_names,colWidths = [0.2, 0.2,0.2,0.2],loc='bottom')
#			the_table.auto_set_font_size(False)
#			the_table.set_fontsize(15.5)

		elif 'pvs' in key:

			self.set_pv_key(key['pvs'][0]) 
			my_pvs,sim_pvs = key['pvs']
			self.set_pv_key(my_pvs)
			sn = float(len(my_pvs))/len(sim_pvs)
			col_names, row_names = sorted(self.pv_key.keys(),reverse=True),['Found','FalsePos']
			row_data = [len([p for p in my_pvs if p < col_names[i]]) for i in range(len(col_names))]
			row_scr = [] 
			for i in range(len(col_names)):
				if row_data[i] > 0: 	row_scr.append(round((sn*len([p for p in sim_pvs if p < col_names[i]]))/float(row_data[i]),3))
				else:			row_scr.append(0)


			col_widths = [0.11, 0.11,0.11,0.11,0.11,0.11,0.10,0.10, 0.10, 0.10]

			colors = plt.cm.hsv(np.linspace(0, 0.5, len(row_names)))[-1::-1]
			table2 = self.ax.table(cellText=[row_data,row_scr],rowLabels=row_names,rowColours=colors,colLabels=col_names,colWidths = col_widths,loc='bottom')


		table2.auto_set_font_size(False)
		table2.set_fontsize(12.5)


		if 'interest' in key and 'preds' in key: 
			it,preds = key['interest'], key['preds']
			subs = [p.split('=')[-1] for p in key['preds'] if p.split('=')[0] == key['interest']]
			if len(subs) == 1 and subs[0] == it: self.ax.set_title('INTEREST: '+key['interest']+' (continuous)',fontsize=20,y=0.99,fontweight='bold')
			elif len(subs) < 4: self.ax.set_title('INTEREST: '+key['interest']+'   ('+",".join(subs)+')',fontsize=15,y=0.92,fontweight='bold')
			else:	            self.ax.set_title('INTEREST: '+key['interest']+'   ('+",".join(subs)+')',fontsize=12,y=0.92,fontweight='bold')

		#if 'title' in key: 	self.ax.set_title(key['title'],fontweight='bold',y=1.05,x=0.0,horizontalalignment='left',fontsize=14)
		self.ax.axis('off') 
		return self


	def save_mfig(self,mt,preds,key={}):


		p_key, f_key, m_str, x = dd(int), {},  mt+':  y  \propto  ', 1
		for p in preds:
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


 		plt.subplots_adjust(left=0.04, bottom=0.05, right=0.96, top=0.88,wspace=0.15,hspace=0.45)
		plt.suptitle('$'+m_str+'$',fontsize=25, fontweight='bold')
		if 'axis_off' in key: plt.axis('off') 


		
		self.fig.savefig(self.options.prefix+'_'+mt+'_predictorModelAnalysis_'+'_'.join(self.options.predictors)+'.pdf', dpi=200) 
	
		if self.options.show: 
			plt.show() 







