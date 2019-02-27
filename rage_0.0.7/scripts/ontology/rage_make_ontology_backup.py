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
from sklearn.preprocessing import RobustScaler
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import Circle as Circ
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
from math import fabs

from sklearn.linear_model import LinearRegression
import random

from sklearn import manifold


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



fp = os.path.dirname(os.path.realpath(__file__))


ONTOLOGY_KEY="/home/tade/rage_0.0.7/rage_0.0.7/rage/data/GENE_ONTOLOGY.key"
ONTOLOGY_TABLE="/home/tade/rage_0.0.7/rage_0.0.7/rage/data/ONTOLOGY_TABLE.txt"


ONTOLOGY_KEY=fp+"/../../data/GENE_ONTOLOGY.key"
ONTOLOGY_TABLE=fp+"/../../data/ONTOLOGY_TABLE.txt"






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





class OntologyPlot:
        def __init__(self,options,dex_key):

                self.options = options
                self.VERBOSE = True
        
		self.DUPES = options.dupes
		self.MIN_CLOUD_DIST = 4
		self.GROUP_CLOUD_DIST = 6
		self.MAX_CLOUD_DIST = 8
		self.MATCH_DIST=15
		self.PRIMARY_DIST=20
		self.INTERMEDIATE_DIST=30
		self.ISLAND_DIST=35
		self.FCFACT=50

		self.minX,self.maxX,self.minY,self.maxY = 0,0,0,0 

        
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
		self.color_key = {'CR': 'cyan', 'MINI': 'lime','AB': 'red','EB': 'blue','MEGA': 'blue'} 


		self.ax = plt.subplot2grid((1,1), (0,0), rowspan = 1, colspan = 1)
		self.color_offset = 0
		self.x,self.y = 0,100 

		self.maxX,self.minY = 0,100 
		self.dex = dex_key 
		self.group_locs = {}

		self.primary_locations = {} 
		self.all_gene_locations = dd(list) 
 
		self.LABEL_LOC = 'TOP'


	def finish2(self,out='BASIC',title='figure',xLen=2,yLen=3,kind='RAW'):
		out_name = out+'.png'
		if kind == 'REL': 	
			plt.suptitle(title+' (RELATIVE SIZE)',fontsize=25,fontweight='bold') 
		else: plt.suptitle(title,fontsize=25,fontweight='bold') 
   		plt.subplots_adjust(left=0.05, bottom=0.075, right=0.95, top=0.92,wspace=0.10,hspace=0.45)
		
#		plt.show() 
		plt.savefig(out_name,dpi=200) 
		plt.clf() 

                self.xLen, self.yLen = xLen,yLen 
		self.xLoc,self.yLoc = 0, 0 
		self.xOffset = 0 







	def make_gene_circle(self,dist,num,origin=(0,0),start_value=None):

		if num == 0: return [] 

		xMult = 1/(num/2.0)
		rads = [math.pi*x*xMult for x in range(0,num)] 

		
		if start_value != None: 
			if start_value < 0: start_value = (2*math.pi)+start_value
			elif start_value > (2*math.pi): start_value -= (2*math.pi)
			rads = [r for r in rads if r >= start_value] + [r for r in rads if r if r<=start_value]


		pts = [] 
		for z in  [math.pi*x*(1.0/(num/2.0)) for x in range(0,num)]:					
			my_dist = dist+np.random.normal(0,fabs(dist)/6.0)


			px = origin[0]+math.cos(z)*my_dist
			py = origin[1]+math.sin(z)*my_dist
			
			pxj = px + np.random.normal(0,fabs(px)/20.0)
			pyj = py + np.random.normal(0,fabs(py)/20.0)

			#pts.append((px,py))
			pts.append((pxj,pyj))

		return pts
	

		pts = [(origin[0]+math.cos(z)*dist,origin[1]+math.sin(z)*dist) for z in  [math.pi*x*(1.0/(num/2.0)) for x in range(0,num)]]
		



		return pts



	def draw_gene(self,loc,g,CENTER=None):

		fx,fy = loc[0],loc[1]
		self.all_gene_locations[g].append((fx,fy))
		HA,VA='center','center'

		if   fx > self.maxX: self.maxX = fx 
		elif fx < self.minX: self.minX = fx 

		if   fy > self.maxY: self.maxY = fy 
		elif fy < self.minY: self.minY = fy 

		if CENTER != None:
			if fy >   CENTER[1]+1: VA='bottom'
			elif fy < CENTER[1]-1: VA='top'

			if fx >   CENTER[0]+1: HA='left'
			elif fx < CENTER[0]-1: HA='right'


		if g not in self.primary_locations: 
			self.primary_locations[g] = fx,fy 
			self.ax.text(fx,fy,g.split(';')[-1],horizontalalignment=HA,verticalalignment=VA)
			self.ax.scatter(fx,fy,marker='o',color='black',s=50,alpha=0.5) 

		elif self.DUPES:
			self.ax.text(fx,fy,g.split(';')[-1],horizontalalignment=HA,verticalalignment=VA)
			self.ax.scatter(fx,fy,marker='o',color='black',s=50,alpha=0.5) 

		if CENTER != None:	
			if self.DUPES: self.ax.plot([CENTER[0],fx],[CENTER[1],fy],zorder=0,color='k',linewidth=0.3,alpha=0.4)
			else:          self.ax.plot([CENTER[0],self.primary_locations[g][0]],[CENTER[1],self.primary_locations[g][1]],zorder=0,color='k',linewidth=0.3,alpha=0.4)
					

		if self.dex.pv[g] < 0.05 and self.dex.top[g] in self.color_key:
			self.ax.scatter(fx, fy, marker='o',color=self.color_key[self.dex.top[g]],s=self.FCFACT*self.dex.fc[g],alpha=0.5)





		return 



	def draw_group(self,loc,grp):
		fx,fy = loc[0],loc[1] 
					
		gstr = grp.split('_')
		gname = " ".join([gs for gs in gstr[0:int(len(gstr)/2)]])+'\n'+" ".join([gs for gs in gstr[int(len(gstr)/2)::]])
		self.ax.add_patch(Circ((fx, fy), 0.6,facecolor='pink',zorder=2))
		self.ax.text(fx,fy,gname,horizontalalignment='center',verticalalignment='bottom')

		self.group_locs[grp] = loc 

		return



	def draw_cluster_to_primary_genes(self,loc,genes):

		fx,fy = loc[0],loc[1]
		for gene in genes:
			if gene in self.primary_locations:
				self.ax.plot([self.primary_locations[gene][0],fx],[self.primary_locations[gene][1],fy],zorder=0,color='k',linewidth=0.3,alpha=0.4)

		return






	def iterative_island_draw(self,cluster_groups,island_pts): 

		share_key,island_genes = dd(list),[]
		for r in cluster_groups:
			pIsland = [g for g in self.group_key[r] if g in island_genes]
			r_genes = [g for g in self.group_key[r] if g not in island_genes]
			pG = [g for g in r_genes if g in self.primary_locations.keys()]
			fG = [g for g in r_genes if g not in self.primary_locations.keys()] 
			tScr,tj,tPt = sorted([(sum([np.linalg.norm(np.array(self.primary_locations[pg])-np.array(pt)) for pg in pG]),j,pt) for j,pt in enumerate(island_pts) if pt != False])[0]
			island_pts[tj] = False 
			self.draw_group(tPt,r) 
			if self.DUPES:
				my_genes = self.group_key[r] 
				island_gene_pts = self.make_gene_circle(self.MIN_CLOUD_DIST,int(len(self.group_key[r])),origin=tPt)
				for (g,(a,b)) in zip(self.group_key[r],island_gene_pts): self.draw_gene((a,b),g,CENTER=tPt) 
				island_genes.extend(fG) 

			else:
				self.draw_cluster_to_primary_genes(tPt,pG)
				island_gene_pts = self.make_gene_circle(self.MIN_CLOUD_DIST,int(len(fG)),origin=tPt)
				for (g,(a,b)) in zip(fG,island_gene_pts): self.draw_gene((a,b),g,CENTER=tPt) 
				self.draw_cluster_to_primary_genes(tPt,pIsland+fG)
				island_genes.extend(fG) 
		return 



	def draw_multi_cluster(self,cluster,cluster_pts,genes,tPt):



		for g_name,pt in zip(cluster,cluster_pts):	self.draw_group(pt,g_name) 

		gene_pts = self.make_gene_circle(self.GROUP_CLOUD_DIST,len(genes),origin=(tPt[0],tPt[1]))					
		for gene,fp in zip(genes,gene_pts):
			self.draw_gene(fp,gene)
			for group in [c for c in cluster if gene in self.group_key[c]]: self.ax.plot([self.group_locs[group][0],fp[0]],[self.group_locs[group][1],fp[1]],zorder=0,color='k',linewidth=0.3,alpha=0.4)
					 					
		return

	def draw_multi_group(self,multi_group): 

		names, locs = [mg[0] for mg in multi_group],[mg[-1] for mg in multi_group]
		gc = cc([a for b in [self.group_key[n] for n in names] for a in b]) 
		centroid = np.mean([mg[-1][0] for mg in multi_group]),np.mean([mg[-1][1] for mg in multi_group])
		cluster_pts  = self.make_gene_circle(5.0,len(names),origin=centroid) 
		shared_genes = [g for g in gc if gc[g] > 1] 

		for pt in cluster_pts:
			group  = names[sorted([(np.linalg.norm(np.array(pt)-np.array(multi_group[i][-1])),i) for i in range(len(multi_group)) if multi_group[i][0] not in self.group_locs])[0][1]]
			self.draw_group(pt,group) 
			uG  = [g for g in self.group_key[group] if self.obs[g] == 1]
			fG  = [g for g in self.group_key[group] if self.obs[g] > 1 and gc[g] == 1]
			if len(uG) > 1:	
				uniq_pts   =   self.make_gene_circle(10,len(fG),origin=pt,start_value = math.atan(pt[1]/pt[0])-(math.pi/4))
				for p,g in zip(uniq_pts,uG):	self.draw_gene(p,g,CENTER=pt)
			if len(fG) > 1:
				far_pts   =    self.make_gene_circle(30,len(fG),origin=pt,start_value = math.atan(pt[1]/pt[0])-(math.pi/6))
				for p,g in zip(far_pts,fG):	self.draw_gene(p,g,CENTER=pt)

		if self.DUPES: 
			cluster_uniq, cluster_multi = [g for g in shared_genes if self.obs[g] == gc[g]] , [g for g in shared_genes if self.obs[g] > gc[g]] 
		else:
			cluster_uniq, cluster_multi = [g for g in shared_genes if g not in self.primary_locations and self.obs[g] == gc[g]] , [g for g in shared_genes if g not in self.primary_locations and self.obs[g] > gc[g]] 

		if len(cluster_uniq) > 0:
			uniq_pts   =   self.make_gene_circle(20,len(cluster_uniq),origin=centroid,start_value = math.atan(centroid[0]/centroid[1])-(math.pi/2))
			for p,g in zip(uniq_pts,cluster_uniq):	self.draw_gene(p,g,CENTER=centroid)
		if len(cluster_multi) > 0: 
			far_pts   =    self.make_gene_circle(40,len(cluster_multi),origin=centroid,start_value = math.atan(centroid[0]/centroid[1])-(math.pi/5))
			for p,g in zip(far_pts,cluster_multi):	self.draw_gene(p,g,CENTER=centroid)

		return




	def draw_single_group(self,group,loc):
		self.draw_group(loc,group)
		genes = self.group_key[group]



		if self.DUPES: 
			cluster_uniq, cluster_multi = [g for g in genes if self.obs[g] == 1] , [g for g in genes if self.obs[g] > 1]
		else:	
			cluster_uniq, cluster_multi = [g for g in genes if self.obs[g] == 1 and g not in self.primary_locations] , [g for g in genes if self.obs[g] > 1 and g not in self.primary_locations]



		if len(cluster_uniq) > 0:
			uniq_pts   =   self.make_gene_circle(20,len(cluster_uniq),origin=loc,start_value = math.atan(float(loc[0])/(0.01+loc[1]))-(math.pi/2))
			for p,g in zip(uniq_pts,cluster_uniq):	self.draw_gene(p,g,CENTER=loc)
		if len(cluster_multi) > 0: 
			far_pts   =    self.make_gene_circle(40,len(cluster_multi),origin=loc,start_value = math.atan(float(loc[0])/(0.01+loc[1]))-(math.pi/5))
			for p,g in zip(far_pts,cluster_multi):	self.draw_gene(p,g,CENTER=loc)

		return



		

		



        def run_clustering(self): 

#		self.gene_key, self.group_key = gene_key, group_key 
		group_names,  gene_membership = self.group_key.keys(), {} 
		SHARED, grp_clusters = dd(bool), dd(list)

		group_shares = [(sum([self.obs[g]-1 for g in self.group_key[group_names[i]]]),sum([self.obs[g] for g in self.group_key[group_names[i]]]),group_names[i]) for i in range(len(group_names))]

		islands, share_sort = [g[3] for g in group_shares if g[0] == 0] , sorted([g for g in group_shares if g[0] > 0],reverse=True) 
		init_group, rem_groups= share_sort[0][-1], [s[-1] for s in share_sort[1::]]
		init_genes = sorted([(len(self.gene_key[gene]),gene) for gene in self.group_key[init_group]],reverse=True)
		pts = self.make_gene_circle(7,len(init_genes))
	
		#self.draw_group((0,0),init_group)
		#for (gL,g),p in zip(init_genes,pts):
			#gene_membership[g] = init_group
			#self.draw_gene(p,g,CENTER=(0,0)) 

		xOffset,yOffset = self.mds_locs[init_group] 
		offset_locs = {grp: (self.mds_locs[grp][0]-self.mds_locs[init_group][0],self.mds_locs[grp][1]-self.mds_locs[init_group][1]) for grp in self.mds_locs}
		
		groups = [init_group] + [grp for grp in offset_locs.keys() if grp != init_group]



		self.MIN_DIST = 2

		multi_group,multi_idxs, group_idxs = [(init_group,offset_locs[init_group])], [0] , range(1,len(groups))

		while True:
			GROWING=False 
			for i in group_idxs:
				for grp,loc in multi_group:
					if i not in multi_idxs and np.linalg.norm(np.array(loc)-np.array(offset_locs[groups[i]])) < self.MIN_DIST:
						multi_group.append((groups[i],offset_locs[groups[i]]))
						multi_idxs.append(i) 
						GROWING=True 
			 			break
					
			if not GROWING: 
				group_idxs = [ix for ix in group_idxs if ix not in multi_idxs] 

				if len(multi_group) > 1:	self.draw_multi_group(multi_group) 
				else:				self.draw_single_group(multi_group[0][0],multi_group[0][1])
				multi_group,multi_idxs = [],[]		
				if len(group_idxs) == 0: break 
				multi_group = [(groups[group_idxs[0]],offset_locs[groups[group_idxs[0]]])]
				multi_idxs  = [group_idxs[0]]





		if not self.DUPES:
			for group in groups:
				for gene in self.group_key[group]:
					print group,gene
					px,py = self.primary_locations[gene]  
					gx,gy = self.group_locs[group] 
					self.ax.plot([px,gx],[py,gy],zorder=0,color='k',linewidth=0.3,alpha=0.4)


		xLocs,yLocs = [loc[0] for loc in offset_locs.values()],[loc[1] for loc in offset_locs.values()] 

		self.ax.set_xlim(self.minX-5,self.maxX+5)
		self.ax.set_ylim(self.minY-5,self.maxY+10)



		



        def run_clustering2(self): 

#		self.gene_key, self.group_key = gene_key, group_key 
		group_names,  gene_membership = self.group_key.keys(), {} 
		SHARED, grp_clusters = dd(bool), dd(list)

		group_shares = [(sum([self.obs[g]-1 for g in self.group_key[group_names[i]]]),sum([self.obs[g] for g in self.group_key[group_names[i]]]),group_names[i]) for i in range(len(group_names))]

		islands, share_sort = [g[3] for g in group_shares if g[0] == 0] , sorted([g for g in group_shares if g[0] > 0],reverse=True) 
		init_group, rem_groups= share_sort[0][-1], [s[-1] for s in share_sort[1::]]
		init_genes = sorted([(len(self.gene_key[gene]),gene) for gene in self.group_key[init_group]],reverse=True)
		pts = self.make_gene_circle(7,len(init_genes))
	
		self.draw_group((0,0),init_group)
		for (gL,g),p in zip(init_genes,pts):
			gene_membership[g] = init_group
			self.draw_gene(p,g,CENTER=(0,0)) 

		print init_group, self.group_locs[init_group]
		xOffset,yOffset = self.mds_locs[init_group] 
		for grp in self.mds_locs:
			print grp, self.mds_locs[grp]  

		
		gUniq,gFirst,gShared,gPrev = dd(list),dd(list),dd(list),dd(list)  
		share_key = dd(lambda: dd(int)) 
		for grp in rem_groups:
			for gene in [g for g in self.group_key[grp]]:		
				if gene in self.primary_locations.keys(): gPrev[grp].append(gene)	
				elif gene in gene_membership.keys(): 	
									gShared[grp].append(gene) 
									share_key[gene_membership[gene]][grp]+=1
				else: 
					if len(self.gene_key[gene]) == 1:	gUniq[grp].append(gene) 
					else: 				gFirst[grp].append(gene) 
					gene_membership[gene] = grp

		for grp in rem_groups:
			share_cands = sorted([(share_key[grp][sg],sg) for sg in share_key[grp] if (share_key[grp][sg] > 3 and (1+len(gFirst[sg])+len(gUniq[sg])) < share_key[grp][sg])],reverse=True)[0:4]
			for s in share_cands:
				grp_clusters[grp].append(s[1]) 
				SHARED[s[1]] = True
		
		rem_islands,k = [grp for grp in rem_groups if len(gFirst[grp])==0 and SHARED[grp] == False],0
		rem_groups = [g[1] for g in sorted([(len(gFirst[grp]),grp) for grp in rem_groups if len(gFirst[grp])>0 and SHARED[grp] == False],reverse=True)]
		rem_share_pts = [p for n,p in enumerate(self.make_gene_circle(self.PRIMARY_DIST,2*int(len(rem_groups)))) if n % 2 == 0]
		rem_offset_pts = [p for n,p in enumerate(self.make_gene_circle(self.MATCH_DIST,2*int(len(rem_groups)))) if n % 2 != 0 and n > 0] 

		while True:
			while k < len(rem_groups) and rem_groups[k] ==False: k+=1 		
			if k == len(rem_groups): break 	

			grp, rem_groups[k] = rem_groups[k] , False 
			pG,uG,fG,sG = gPrev[grp],gUniq[grp],gFirst[grp],gShared[grp]
			
			if len(grp_clusters[grp]) > 0:

				cluster = grp_clusters[grp] + [grp]
				aG = list(set([a for b in [self.group_key[gp] for gp in cluster] for a in b]))
				pG = list(set([a for b in [gPrev[gp] for gp in cluster] for a in b]))
				uG = list(set([a for b in [gUniq[gp] for gp in cluster] for a in b]))+list(set([a for b in [gFirst[gp] for gp in cluster] for a in b]))
				tScr,tj,tPt = sorted([(sum([np.linalg.norm(np.array(self.primary_locations[pg])-np.array(pt)) for pg in pG]),j,pt) for j,pt in enumerate(rem_share_pts) if pt != False])[0]
				rem_share_pts[tj] = False
				cluster_pts  = self.make_gene_circle(1.0,len(cluster),origin=(tPt[0],tPt[1]))					
				
				if self.DUPES:	self.draw_multi_cluster(cluster,cluster_pts,aG,tPt)
				else:		self.draw_multi_cluster(cluster,cluster_pts,uG,tPt) 


			else:
				if len(gPrev[grp]) > 0:
					tScr,tj,tPt = sorted([(sum([np.linalg.norm(np.array(self.primary_locations[pg])-np.array(pt)) for pg in pG]),j,pt) for j,pt in enumerate(rem_share_pts) if pt != False])[0]
					rem_share_pts[tj] = False
				else:
					tj,tPt = [(j,pt) for (j,pt) in enumerate(rem_share_pts) if pt != False][0] 
				self.draw_group(tPt,grp) 
				if self.DUPES:
					island_gene_pts = self.make_gene_circle(5,int(len(self.group_key[grp])),origin=(tPt[0],tPt[1]))
					for (g,(a,b)) in zip(self.group_key[grp],island_gene_pts):	 self.draw_gene((a,b),g,CENTER=tPt)

				else:
					grp_far = self.make_gene_circle(4.0,len(uG+fG),origin=(tPt[0],tPt[1]),start_value = math.atan(tPt[1]/tPt[0])-(math.pi/6))
					grp_close = self.make_gene_circle(2.5,len(uG+fG),origin=(tPt[0],tPt[1]),start_value = math.atan(tPt[1]/tPt[0])-(math.pi/6))
					for c,p,g in zip(grp_close,grp_far,uG+fG):
						if g in uG: 	self.draw_gene(c,g,CENTER=tPt)
						else:		self.draw_gene(p,g,CENTER=tPt)
					self.draw_cluster_to_primary_genes(tPt,pG)



		if not self.DUPES: 
			size_dict = {r: len([g for g in self.group_key[r] if g not in self.primary_locations.keys()]) for r in rem_islands} 
			size_med =  np.median(size_dict.values()) 	
			rem_niches = [rn[1] for rn in sorted([(v,r) for (r,v) in size_dict.items() if v<3 or (v*2 <= size_med and v <= 10)])[0:len(rem_offset_pts)]]
			rem_islands = [r for r in size_dict.keys() if r not in rem_niches]
			self.iterative_island_draw(rem_niches,rem_offset_pts)

		island_pts = [p for n,p in enumerate(self.make_gene_circle(self.ISLAND_DIST,(2*int(len(rem_islands))))) if n % 2 != 0][0:len(rem_islands)]		
		self.iterative_island_draw(rem_islands,island_pts) 

		if self.DUPES:
			for gene,locs in self.all_gene_locations.items():
				for m in range(len(locs)):
					for n in range(m+1,len(locs)):	self.ax.plot([locs[m][0],locs[n][0]],[locs[m][1],locs[n][1]],zorder=0,color='k',linewidth=0.3,alpha=0.4)
				
					
		x_locs = sorted([pt[0] for pt in self.group_locs.values()])
		y_locs = sorted([pt[1] for pt in self.group_locs.values()])


		self.ax.set_xlim(x_locs[0]-(2*self.MIN_CLOUD_DIST),x_locs[-1]+(2*self.MIN_CLOUD_DIST))
		self.ax.set_ylim(y_locs[0]-(2*self.MIN_CLOUD_DIST),y_locs[-1]+(2*self.MIN_CLOUD_DIST))


		



	def finish(self): 



   		plt.subplots_adjust(left=0.01, bottom=0.075, right=0.95, top=0.92,wspace=0.10,hspace=0.45)
		self.ax.axis('off') 
		plt.show() 


	def set_cluster_locations(self,ontology_groups,ontology_genes,ontology_counts): # ID=FOCUS) 
		
		self.group_key, self.gene_key, self.obs = ontology_groups, ontology_genes, ontology_counts
		ontology_dists = {} 
		ontology_names = ontology_groups.keys() 
		for i in range(len(ontology_names)):
			nI = ontology_names[i] 
			gI = ontology_groups[ontology_names[i]]
			cI = float(len(gI))
			ontology_dists[(nI,nI)] = 0.0
			for j in range(i+1,len(ontology_names)):
				nJ = ontology_names[j] 
				gJ = ontology_groups[ontology_names[j]]
				dist = len([g for g in gI if g in gJ]) / float(len(list(set(gI+gJ))))
				ontology_dists[(nI,nJ)] = 1.0-dist
				ontology_dists[(nJ,nI)] = 1.0-dist

		d_matrix = np.matrix([[ontology_dists[n,m] for m in ontology_names] for n in ontology_names])

		mds = manifold.MDS(n_components=2,metric=False,max_iter=3000,dissimilarity='precomputed',n_jobs=1)
		
		mds.fit_transform(d_matrix)
		group_pts = mds.embedding_
		
		gX,gY = np.array([p[0] for p in group_pts]),np.array([p[1] for p in group_pts])
		scaler = MinMaxScaler(feature_range=(0,150))
		scaler = RobustScaler() #feature_range=(0,150))

 

                mX,mY =[p[0]*100 for p in scaler.fit_transform(np.array(gX).reshape(-1,1))] ,[p[0] * 100 for p in scaler.fit_transform(np.array(gY).reshape(-1,1))]





		group_pts = [(mX[i],mY[i]) for i in range(len(group_pts))]

	



		self.mds_locs = {ontology_names[i]:  [p[0],p[1]] for i,p in enumerate(group_pts)} 





class Ontology:
        def __init__(self,options):


		self.table = {} 

		idx, self.gene_idx, self.gene_key, self.ontology_key = 0,{}, dd(list), dd(list)
		
		for line in open(ONTOLOGY_KEY):
			line = line.split() 
			gene,opts = line[0],line[1].split(',') 
			self.gene_key[gene] = opts 
			for k in opts: self.ontology_key[k].append(idx) 
			self.gene_idx[idx] = gene
			idx += 1 

		for line in open(ONTOLOGY_TABLE): 
			line = line.strip().split('|') 
			idx,go,ont = line[0].split() 
			self.table[idx]  = [go,ont,line[-1]]	
			






class DexGenes:
        def __init__(self,group_names,options=None):


		self.groups = group_names 
		self.genes = [] 
		self.key = dd(lambda: {}) 
		self.sig, self.pv,self.fc, self.ranks, self.top = {}, {}, {}, {}  , {} 


	def add_gene_cnts(self,gene,group_counts): 

		
			
		group_data = sorted([[np.mean(group_counts[i]),len([x for x in group_counts[i] if x>0])/float(len(group_counts[i])),grp,group_counts[i]] for i,grp in enumerate(self.groups)],reverse=True)


		group_pv = stats.ttest_ind(group_data[0][-1],group_data[-1][-1])[1]
		group_fc = group_data[0][0] / (0.01+group_data[-1][0])

		self.genes.append(gene) 
		self.pv[gene] = group_pv 
		self.fc[gene] = group_fc
		self.ranks[gene] = [gd[2] for gd in group_data] 
		self.top[gene] = group_data[0][2] 
	






def run_script(options):

	if options.id == None or len(options.id.split('='))!=2:
		print 'warning: need an id'
		ID = 'CT'
		FOCUS='MEGA'
		IGNORE_NAME = 'MINI'

	else:
		ID = options.id.split('=')[0] 
		FOCUS = options.id.split('=')[1] 
		IGNORE_NAME = None

	ontology = Ontology(options) 
	key = dd(lambda: {}) 

	for line in open(options.key):
		line = line.split() 
		if line[0] == '---': headers = line  
		else:
			for i in range(1,len(line)): key[headers[i]][line[0]] = line[i]  
				

	cnts, groups = {} , dd(list)  
	for line in open(options.cnts): 
		line = line.split() 
		if line[0] == '---': 
			samples = line[1::] 
			for i,s in enumerate(samples): 
				if s in key[ID] and key[ID][s] != 'NA': 
					groups[key[ID][s]].append(i) 
			group_names = sorted(groups.keys()) 		
			dex = DexGenes(group_names) 


		else:
			gene,vals = line[0],[float(x) for x in line[1::]]
			vals = [log(v+1,2) for v in vals]
			if gene not in ontology.gene_key.keys():	continue 
			if gene in ontology.gene_key.keys(): 
				cnts[gene] = vals 
				dex.add_gene_cnts(gene,[[vals[i] for i in groups[k]] for k in group_names])

#			if len(cnts) > 2000: break 
	
	kCands = [] 
	for k,members in ontology.ontology_key.items():
#		print ontology.ontology_table[k]
		gm = [ontology.gene_idx[m] for m in members if ontology.gene_idx[m] in cnts] 
		gm_pv = [dex.pv[g] for g in gm]
		sig_genes = [g for g in gm if dex.pv[g] < 0.01]
		if len(gm) < 4 or len(sig_genes) < 5: continue 

		if len(gm) < 5: continue 
		if min(gm_pv) > 0.00000001: continue

		sig_tops = [dex.top[g] for g in sig_genes] 		
		topName,topCnt = sorted(cc(sig_tops).items(),key=lambda X: X[1],reverse=True)[0] 
		topRate = topCnt/float(len(sig_tops))
		topGlobal = topCnt / float(len(gm)) 

		if topRate < 0.5 or topGlobal < 0.2: continue 
		if len(gm) - len(sig_genes) > 25: continue 
		if topCnt < 3: continue 

		topPvs = [dex.pv[g] for g in sig_genes if dex.top[g] == topName]
		minP = min(topPvs) 	
		kCands.append([minP,topRate,topGlobal,topCnt,topName,k]) 

	kCands.sort() 
	ADDED = 0 
#	CAT_NAME = 'MINI'
	ontology_genes = dd(list) 
	ontology_groups = dd(list) 
	ontology_counts = dd(int) 
	ontology_types  = dd(lambda: dd(list))
	op =  OntologyPlot(options,dex)
	
	for i,(minP,topRate,topGlobal,topCnt,topName,k) in enumerate(kCands):
		gm = [ontology.gene_idx[m] for m in ontology.ontology_key[k] if ontology.gene_idx[m] in cnts]
		sig_genes = [g for g in gm if dex.pv[g] < 0.05]
		sig_tops = [dex.top[g] for g in sig_genes] 		


		if IGNORE_NAME == None:
			if topName != FOCUS: continue 		

		if len(gm) > 40: continue 

#		if ontology.table[k][1][0:3] in ['mit','dis','cel','cen','cil','pro']: continue

		ot_term,ot_name,ot_str = ontology.table[k]

		if ot_name.split('_')[0] in ['mitotic','mitosis']: continue 	




		for g in gm: 
			ontology_counts[g]+=1 
			ontology_genes[g].append(ot_name) 

		ontology_groups[ot_name] = gm 
		ontology_types[topName][ot_name] = gm 

		ADDED+=1
		if ADDED > 15:
			break		

	for k in ontology_types:
		op.set_cluster_locations(ontology_groups,ontology_genes,ontology_counts) #,ID=FOCUS) 
		op.run_clustering() 
		op.finish() 

	#op.add_data(ontology_groups,ontology_genes,ontology_counts,ID=FOCUS) 
#	op.add_data(ontology_groups,ontology_genes,ontology_counts,ID='MEGA') 

	sys.exit() 
	

				















if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("--id", default = None, type='string', help="Category Id")
	parser.add_option("-c", "--cnts", default = None, type='string', help="Output Filename Prefix")
 	parser.add_option("--dupes", default = False, action='store_true', help="Output Filename Prefix")



	(options, args) = parser.parse_args()

	run_script(options)	













