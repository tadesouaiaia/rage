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

ONTOLOGY_KEY=fp+"/../../data/GENE_ONTOLOGY_MAJOR.key"
ONTOLOGY_TABLE=fp+"/../../data/ONTOLOGY_TABLE_MAJOR.txt"







G_GROW=["activation_of_axon_extension","activation_of_axonogenesis","activation_of_neuron_differentiation","axonal_fasciculation","axon_development","axon_extension","axon_growth_cone","axon_guidance","basal_dendrite","dendrite_branch","dendrite_development","growth_cone","inhibition_of_axonogenesis","inhibition_of_neuron_differentiation","neural_crest_cell_migration","neuron_development","neuron_differentiation","neuron_fate_commitment","neuron_guidance","positive_regulation_of_neuron_guidance","proton_transport","somatodendritic_compartment",]

G_GROW.extend(["structural_constituent_of_myelin_sheath"])
G_GROW.extend(["asymmetric_neuroblast_division","ATP_anabolism","cerebral_cortex_GABAergic_interneuron_migration","embryonic_hemopoiesis","microglial_cell_activation","neurofibrillary_tangle","positive_regulation_of_neural_precursor_cell_proliferation"])

G_SPLICE=["AT-AC_spliceosomal_complex","catalytic_step_2_spliceosome","nuclear_import","spliceosome","spliceosome_assembly","splice_site_selection",]


G_TUBE=["activation_of_microtubule_polymerization","axon_development","cytoplasmic_microtubule","cytoplasmic_microtubule_organisation","microtubule","microtubule_associated_complex","microtubule-based_movement","microtubule-based_process","microtubule_bundling","microtubule_cytoskeleton","microtubule_dynamics","microtubule_nucleation","microtubule_plus_end","microtubule_plus-end-directed_vesicle_distribution","microtubule_rescue","structural_constituent_of_cytoskeleton",]


G_TUBE.extend(["actin_filament_branch_nucleation"])

G_SYN=["activation_of_synapse_assembly","anterograde_synaptic_vesicle_transport","calcium_channel_regulator_activity","chloride_channel_activity","chloride_channel_inhibitor_activity","GABA-A_receptor_complex","glutamate_receptor_activity","inhibitory_extracellular_ligand-gated_ion_channel_activity","integral_component_of_presynaptic_membrane","kainate_selective_glutamate_receptor_activity","kainate_selective_glutamate_receptor_complex","Kv_channel_clustering","metabotropic_glutamate_receptor,_adenylyl_cyclase_inhibiting_pathway","modulation_of_synaptic_transmission","postsynaptic_membrane_assembly","postsynaptic_neurotransmitter_receptor_endocytosis","postsynaptic_neurotransmitter_receptor_endosomal_trafficking","receptor_clustering","regulation_of_postsynaptic_density_assembly","regulation_of_postsynaptic_membrane_potential","regulation_of_presynapse_assembly","sodium_channel_inhibitor_activity","synaptic_vesicle_clustering","vesicle"]

G_SYN.extend(["synapse_maturation","anterograde_axonal_transport","dendrite_membrane"])
G_SYN.extend(["dendritic_spine_development","membrane_depolarization_during_action_potential","postsynaptic_density,_intracellular_component"])



G_DIFF=(["chemokine-mediated_signaling_pathway","glial_cell_development","glia_proliferation","inhibition_of_amyloid_fibril_assembly","microglial_cell_proliferation","negative_regulation_of_hippocampal_neuron_apoptotic_process","neuroepithelial_cell_differentiation","neurulation","noradrenergic_neuron_differentiation","spinal_cord_development"])









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
        def __init__(self,options,dex_key,dex_groups):

                self.options = options
                self.VERBOSE = True
        
		self.xOffset,self.yOffset = 0,0  
		self.DUPES = options.dupes
		self.MIN_CLOUD_DIST = 5
		self.GROUP_CLOUD_DIST = 6
		self.MAX_CLOUD_DIST = 10
		self.MATCH_DIST=10
		self.PRIMARY_DIST=15
		self.ISLAND_DIST=20

		self.minX,self.maxX,self.minY,self.maxY = 0,0,0,0 

        
                sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(4.5, 2.5)
                self.fig.set_facecolor('white') 
                self.fig.patch.set_facecolor('white')
                matplotlib.rcParams['savefig.facecolor'] = 'white'
                matplotlib.rcParams['ytick.labelsize'] = 7.5
                
                seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})

                #self.fig.patch.set_facecolor('lightgrey')
                self.fig.patch.set_facecolor('white')
		self.color_key = {'MULTI': 'orange', 'ES': 'cyan', 'CR': 'cyan', 'MINI': 'lime','AB': 'red','EB': 'blue','MEGA': 'blue'} 

		self.xLen,self.yLen = 1,1	
		if len(dex_groups) < 4:
			self.xLen,self.yLen = 1, len(dex_groups) 

		elif len(dex_groups) == 4:
			self.xLen,self.yLen = 2, 2
		elif len(dex_groups) > 4:
			self.xLen,self.yLen = 2, 3

		self.xLen,self.yLen = 2,2
	
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (0,0), rowspan = 1, colspan = 1)
		self.a1 = plt.subplot2grid((self.xLen,self.yLen), (0,1), rowspan = 1, colspan = 1)
		self.a2 = plt.subplot2grid((self.xLen,self.yLen), (1,0), rowspan = 1, colspan = 2)
		self.plt_num = 0 
		self.xLoc, self.yLoc = 0,0

		self.color_offset = 0
		self.x,self.y = 0,100 

		self.LABEL_LOC = 'TOP'
		self.maxX,self.minY = 0,100 
		self.dex = dex_key 
		self.group_locs = {}
		self.primary_locations = {} 
		self.all_gene_locations = dd(list)  
		self.clustered_groups = [] 



	def finish(self): 
		out_name = 'cluster_fig.png'
		self.ax.axis('off') 
   		#plt.subplots_adjust(left=-0.05, bottom=-0.05, right=0.99, top=0.98,wspace=0.1,hspace=0.1)
		plt.savefig(out_name,dpi=500) 
		plt.show() 





	def make_gene_circle(self,dist,num,origin=(0,0),start_value=None,OFFSET=True):

		if num == 0: return [] 

		xMult = 1/(num/2.0)
		rads = [math.pi*x*xMult for x in range(0,num)] 

		
		if start_value != None: 
			if start_value < 0: start_value = (2*math.pi)+start_value
			elif start_value > (2*math.pi): start_value -= (2*math.pi)
			rads = [r for r in rads if r >= start_value] + [r for r in rads if r if r<=start_value]


		my_pts = [] 
		for z in  [math.pi*x*(1.0/(num/2.0)) for x in range(0,num)]:					
			my_dist = dist+np.random.normal(0,fabs(dist)/40.0)


			px = origin[0]+math.cos(z)*my_dist
			py = origin[1]+math.sin(z)*my_dist
			
			pxj = px + np.random.normal(0,0.5+fabs(px)/25.0)
			pyj = py + np.random.normal(0,0.5+fabs(py)/25.0)

			#pts.append((px,py))
			my_pts.append((pxj,pyj))

	

		#my_pts = [(origin[0]+math.cos(z)*dist,origin[1]+math.sin(z)*dist) for z in  [math.pi*x*(1.0/(num/2.0)) for x in range(0,num)]]
		if OFFSET: 
			pts = [(p[0]+ self.xOffset, p[1] +self.yOffset) for p in my_pts]
			return pts
		else:
			return my_pts


	def get_size_color(self,pv,fc,top_group):

		if top_group not in self.color_key: return False,5,'grey'
		
		clr = self.color_key[top_group] 

		if   pv < 0.00001 and fc > 4:  s = 120 
		elif pv < 0.0001  and fc > 3:  s = 100 
		elif pv < 0.0001   and fc > 2: s = 80
		elif pv < 0.001:              s = 40 
		elif fc < 3: 		      return False,10,clr 
		else:			      return False,8,'grey'
	
		return True,s,clr	




	def draw_gene(self,loc,g,CENTER=None):

		fx,fy = loc[0],loc[1]
		self.all_gene_locations[g].append((fx,fy))
		PLOT,HA,VA=False,'center','center'

		if   fx > self.maxX: self.maxX = fx 
		elif fx < self.minX: self.minX = fx 

		if   fy > self.maxY: self.maxY = fy 
		elif fy < self.minY: self.minY = fy 

		if CENTER != None:
			if fy >   CENTER[1]+1: VA='bottom'
			elif fy < CENTER[1]-1: VA='top'

			if fx >   CENTER[0]+1: HA='left'
			elif fx < CENTER[0]-1: HA='right'

		if self.DUPES or g not in self.primary_locations:
			if g not in self.primary_locations: self.primary_locations[g] = fx,fy 
			PLOT=True

		if PLOT: 
			pv,fc,top_group = self.dex.pv[g], self.dex.fc[g], self.dex.top[g] 
			TEXT,size,clr = self.get_size_color(pv,fc,top_group) 
			
			self.ax.scatter(fx,fy,marker='o',color=clr,s=size,alpha=0.5) 
			if TEXT:	self.ax.text(fx,fy,g.split(';')[-1],fontsize=size/10,color=clr,horizontalalignment=HA,verticalalignment=VA)
				 



		if CENTER != None:	

			if self.DUPES: self.ax.plot([CENTER[0],fx],[CENTER[1],fy],zorder=0,color='k',linewidth=0.1,alpha=0.4)
			else:          self.ax.plot([CENTER[0],self.primary_locations[g][0]],[CENTER[1],self.primary_locations[g][1]],zorder=0,color='k',linewidth=0.1,alpha=0.4)
					




		return 



	def draw_group(self,loc,grp,FORCE=False):
		fx,fy = loc[0],loc[1] 
					
		gstr = grp.split('_')
		gname = " ".join([gs for gs in gstr[0:int(len(gstr)/2)]])+'\n'+" ".join([gs for gs in gstr[int(len(gstr)/2)::]])


		print grp,'yo'
		if grp == 'neuroepithelial_cell_differentiation': loc = (loc[0] -5,loc[1]) 

		if FORCE:
			print 'yup' 
			self.ax.add_patch(Circ((fx, fy), 1.0,facecolor='pink',zorder=2))
			self.ax.text(fx,fy,gname,horizontalalignment='center',fontsize=10,verticalalignment='bottom')

		self.group_locs[grp] = loc 

		return



	def draw_cluster_to_primary_genes(self,loc,genes):

		fx,fy = loc[0],loc[1]
		for gene in genes:
			if gene in self.primary_locations:
				self.ax.plot([self.primary_locations[gene][0],fx],[self.primary_locations[gene][1],fy],zorder=0,color='k',linewidth=0.1,alpha=0.3)

		return






	def iterative_island_draw(self,cluster_groups,island_pts,OFFSET=False): 

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
				island_gene_pts = self.make_gene_circle(self.MIN_CLOUD_DIST,int(len(self.group_key[r])),origin=tPt,OFFSET=False)
				for (g,(a,b)) in zip(self.group_key[r],island_gene_pts): self.draw_gene((a,b),g,CENTER=tPt) 
				island_genes.extend(fG) 

			else:
				self.draw_cluster_to_primary_genes(tPt,pG)
				island_gene_pts = self.make_gene_circle(self.MIN_CLOUD_DIST,int(len(fG)),origin=tPt,OFFSET=False)
				for (g,(a,b)) in zip(fG,island_gene_pts): 
					self.draw_gene((a,b),g,CENTER=tPt) 
				self.draw_cluster_to_primary_genes(tPt,pIsland+fG)
				island_genes.extend(fG) 
		return 



	def draw_multi_cluster(self,cluster,cluster_pts,genes,tPt):



		for g_name,pt in zip(cluster,cluster_pts):	self.draw_group(pt,g_name) 

		gene_pts = self.make_gene_circle(self.GROUP_CLOUD_DIST,len(genes),origin=(tPt[0],tPt[1]),OFFSET=False)					
		for gene,fp in zip(genes,gene_pts):
			self.draw_gene(fp,gene)
			for group in [c for c in cluster if gene in self.group_key[c]]: self.ax.plot([self.group_locs[group][0],fp[0]],[self.group_locs[group][1],fp[1]],zorder=0,color='k',linewidth=0.3,alpha=0.4)
					 					
		return






		

		





        def run_clustering2(self): 

#		self.gene_key, self.group_key = gene_key, group_key 
		group_names,  gene_membership = self.group_key.keys(), {} 
		SHARED, grp_clusters = dd(bool), dd(list)

		group_shares = [(sum([self.obs[g]-1 for g in self.group_key[group_names[i]]]),sum([self.obs[g] for g in self.group_key[group_names[i]]]),group_names[i]) for i in range(len(group_names))]



		share_sort = sorted([g for g in group_shares if g[0] > 0],reverse=True) 
		total_islands   = [g[2] for g in group_shares if g[0] == 0] 


		init_group, rem_groups= share_sort[0][-1], [s[-1] for s in share_sort[1::]]
		init_genes = sorted([(len(self.gene_key[gene]),gene) for gene in self.group_key[init_group]],reverse=True)
		pts = self.make_gene_circle(7,len(init_genes))
	


		self.draw_group((self.xOffset,self.yOffset),init_group)
		for (gL,g),p in zip(init_genes,pts):	
			gene_membership[g] = init_group
			self.draw_gene(p,g,CENTER=(self.xOffset,self.yOffset)) 



		
		gUniq,gFirst,gShared,gPrev = dd(list),dd(list),dd(list),dd(list)  
		share_key = dd(lambda: dd(int))



		print 'init',init_group
 
		for grp in rem_groups:

			print 'rem',grp
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
			share_cands = sorted([(share_key[grp][sg],sg) for sg in share_key[grp] if (share_key[grp][sg] > 3 and (len(gFirst[sg])+len(gUniq[sg])) < 3+share_key[grp][sg])],reverse=True)[0:4]
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
				self.clustered_groups.append(cluster)
				
				for c in cluster: print 'cluster',c


				

				aG = list(set([a for b in [self.group_key[gp] for gp in cluster] for a in b]))
				pG = list(set([a for b in [gPrev[gp] for gp in cluster] for a in b]))
				uG = list(set([a for b in [gUniq[gp] for gp in cluster] for a in b]))+list(set([a for b in [gFirst[gp] for gp in cluster] for a in b]))
				tScr,tj,tPt = sorted([(sum([np.linalg.norm(np.array(self.primary_locations[pg])-np.array(pt)) for pg in pG]),j,pt) for j,pt in enumerate(rem_share_pts) if pt != False])[0]
				rem_share_pts[tj] = False
				cluster_pts  = self.make_gene_circle(3.0,len(cluster),origin=(tPt[0],tPt[1]),OFFSET=False)					
				
				if self.DUPES:	self.draw_multi_cluster(cluster,cluster_pts,aG,tPt)
				else:		self.draw_multi_cluster(cluster,cluster_pts,uG,tPt) 


			else:
				print 'serio',grp

				if len(gPrev[grp]) > 0:
					print 'serioA',grp
					tScr,tj,tPt = sorted([(sum([np.linalg.norm(np.array(self.primary_locations[pg])-np.array(pt)) for pg in pG]),j,pt) for j,pt in enumerate(rem_share_pts) if pt != False])[0]
					rem_share_pts[tj] = False
				else:
					
					tj,tPt = [(j,pt) for (j,pt) in enumerate(rem_share_pts) if pt != False][0] 
					rem_share_pts[tj] = False


				self.draw_group(tPt,grp) 
				if tPt[0] == 0: tPt = (tPt[0] + 0.001, tPt[1]) 
				if self.DUPES:
					island_gene_pts = self.make_gene_circle(5,int(len(self.group_key[grp])),origin=(tPt[0],tPt[1]),OFFSET=False)
					for (g,(a,b)) in zip(self.group_key[grp],island_gene_pts):	 self.draw_gene((a,b),g,CENTER=tPt)

				else:
					
					
					grp_far = self.make_gene_circle(4.0,len(uG+fG),origin=(tPt[0],tPt[1]),start_value = math.atan(tPt[1]/tPt[0])-(math.pi/6),OFFSET=False)
					grp_close = self.make_gene_circle(2.5,len(uG+fG),origin=(tPt[0],tPt[1]),start_value = math.atan(tPt[1]/tPt[0])-(math.pi/6),OFFSET=False)
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
					for n in range(m+1,len(locs)):	self.ax.plot([locs[m][0],locs[n][0]],[locs[m][1],locs[n][1]],zorder=0,color='k',linewidth=0.03,alpha=0.4)
	
					
		x_locs = sorted([pt[0] for pt in self.group_locs.values()])
		y_locs = sorted([pt[1] for pt in self.group_locs.values()])


		if len(total_islands) > 0: 

			xStep = (x_locs[-1]-x_locs[0])/(float(len(total_islands)))
			xStart = x_locs[0]+(xStep/3.0)


			for ti,t in enumerate(total_islands):

				print 'wtf',ti,t
				
				t_genes = self.group_key[t]
				t_mid = xStart + (ti*xStep), y_locs[0] - self.MAX_CLOUD_DIST+1.0
				t_pts  = self.make_gene_circle(self.MIN_CLOUD_DIST,len(self.group_key[t]),origin=t_mid)
				self.draw_group((t_mid[0]+self.xOffset,t_mid[1]+self.yOffset),t) 
				for g,p in zip(t_genes,t_pts): self.draw_gene(p,g,CENTER=(t_mid[0]+self.xOffset,t_mid[1]+self.yOffset))

		x_locs = sorted([pt[0] for pt in self.group_locs.values()])
		y_locs = sorted([pt[1] for pt in self.group_locs.values()])


#		self.ax.set_xlim(x_locs[0]-(2*self.MIN_CLOUD_DIST),x_locs[-1]+(2*self.MIN_CLOUD_DIST))
#		self.ax.set_ylim(y_locs[0]-(2*self.MIN_CLOUD_DIST),y_locs[-1]+(2*self.MIN_CLOUD_DIST))



	def add_outline(self,kp):


		kpv = kp 
		x_srt = sorted(self.group_locs.values())
		y_srt = sorted(self.group_locs.values(),key=lambda X: X[1])

		x_rev = sorted(self.group_locs.values(),reverse=True)
		y_rev = sorted(self.group_locs.values(),reverse=True,key=lambda X: X[1])

		xMin,xMax = x_srt[0][0],x_srt[-1][0]
		yMin,yMax = y_srt[0][1],y_srt[-1][1]

		#print 'hmm',self.plt_num, kp, self.xOffset, self.yOffset,xMin,xMax, yMin,yMax
		FS = 20
		if kp == 'VARIOUS': return 

		elif self.plt_num == 1: self.ax.text(xMin-3,yMax-2,kpv,fontsize=FS)
		elif self.plt_num == 2:	self.ax.text(self.xOffset,yMax+10,kpv,fontsize=FS,horizontalalignment='center')
		elif self.plt_num == 3:	self.ax.text(xMax,yMax,kpv,fontsize=FS)
		elif self.plt_num == 4:	self.ax.text(xMin-5,yMax-5,kpv,fontsize=FS)
		elif self.plt_num == 5:	self.ax.text(xMax,yMax,kpv,fontsize=FS)


		return 

		p1 = [y_srt[0]] 
		for x,y in y_srt[1::]:
			if x<p1[-1][0]: p1.append((x,y))
			if x == xMin: break 

		for x,y in x_srt[1::]: 
			if y > p1[-1][1]: p1.append((x,y))
			if y == yMax: break 

		for x,y in y_rev[1::]:
			if x > p1[-1][0]: p1.append((x,y))
			if x == xMax: break 

		for x,y in x_rev[1::]:
			if y < p1[-1][1]: p1.append((x,y))
			if y == yMin: break 
	
		x_pts,y_pts = [],[] 
		for x,y in p1[0:-1]: 
			if x < 0: x_pts.append(x - 20 + np.random.normal(0,5))
			else:     x_pts.append(x + 20 + np.random.normal(0,5))
			if y < 0: y_pts.append(y - 20 + np.random.normal(0,5))
			else:     y_pts.append(y + 20 + np.random.normal(0,5))
		
		x_pts.append(x_pts[0])
		y_pts.append(y_pts[0])


		self.ax.plot(x_pts,y_pts,'--',color='k')
	
	def add_group_labels(self,kp):



		VA_LOC = {} 
		for c_group in self.clustered_groups:
			v_grps = [vg[1] for vg in sorted([(self.group_locs[grp][1],grp) for grp in c_group],reverse=True)]
			h_grps = [vg[1] for vg in sorted([(self.group_locs[grp][0],grp) for grp in c_group])]

			if len(v_grps) == 2: 
				VA_LOC[v_grps[0]] = 'top'
				VA_LOC[v_grps[1]] = 'bottom'

			elif len(v_grps) == 3: 
				VA_LOC[v_grps[0]] = 'top'
				VA_LOC[v_grps[1]] = 'center'
				VA_LOC[v_grps[2]] = 'bottom'

			elif len(v_grps) == 4:
				continue  

			else:

				for i,grp in enumerate(v_grps):

	
					if i % 3 == 0: VA_LOC[grp]='bottom'
					elif i % 2 == 0: VA_LOC[grp]='top'
					elif i == 0: 	 VA_LOC[grp] = 'bottom'


		self.add_outline(kp)



		for grp,loc in self.group_locs.items():
			#print grp,'loc',str(loc[0])+','+str(loc[1]),",".join(sorted(self.group_key[grp]))

			if grp in VA_LOC: VA = VA_LOC[grp] 
			else: 	          VA = 'center' 

			gstr = grp.split('_')
			gname = " ".join([gs for gs in gstr[0:int(len(gstr)/2)]])+'\n'+" ".join([gs for gs in gstr[int(len(gstr)/2)::]])
			gn = gstr[0] 
			for i in range(1,len(gstr)): 
				g_add = gstr[i] 
				if len(g_add)> 3 and len(gstr[i-1])+len(g_add) > 8: gn+='\n'+g_add
				else:  gn+=' '+g_add

			gname = gn

			f_scr,r_scr,t_scr,topName = self.group_scores[grp] 
			fs = 3+(f_scr*1.3)
			if fs > 13: fs = 13
			clr = self.color_key[topName] 
			#self.ax.add_patch(Circ((fx, fy), 0.6,facecolor='pink',zorder=2))
			self.ax.add_patch(Circ((loc[0], loc[1]), r_scr,facecolor=clr,zorder=2,alpha=0.3))
			#self.ax.text(loc[0],loc[1],gname,fontweight='bold',zorder=9,horizontalalignment='center',color=clr,verticalalignment=VA,fontsize=fs) 
			self.ax.text(loc[0],loc[1],gname,fontweight='bold',zorder=9,horizontalalignment='center',color='k',verticalalignment=VA,fontsize=fs) 
			
#		self.ax.axis('off') 
#		self.xLoc +=1


#		if self.xLoc == self.xLen:
#			self.yLoc+=1
#			self.xLoc = 0 





	def set_cluster_locations(self,kp,ontology_groups,ontology_scores): #,ontology_genes,ontology_counts): # ID=FOCUS) 

		#self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		#self.ax = plt.subplot2grid((self.xLen,self.yLen), (0,0), rowspan = 1, colspan = 1)

		if kp != 'VARIOUS':

			print kp,self.plt_num 
			if self.plt_num == 0: self.xOffset,self.yOffset = -40,55
			elif self.plt_num == 1: self.xOffset,self.yOffset = -3,65
			elif self.plt_num == 2: self.xOffset,self.yOffset = 35,60
			elif self.plt_num == 3: self.xOffset,self.yOffset = -53,12
			elif self.plt_num == 4: self.xOffset,self.yOffset = 38,-2
			elif self.plt_num == 5: self.xOffset,self.yOffset = 0,-100
			else: 			 self.xOffset,self.yOffset = 10,-30
			self.plt_num+=1

		else:
			self.xOffset,self.yOffset = -5,12
	
		self.MIN_CLOUD_DIST = 5
		self.GROUP_CLOUD_DIST = 6
		self.MAX_CLOUD_DIST = 10
		self.MATCH_DIST=10
		self.PRIMARY_DIST=15
		self.ISLAND_DIST=20

		if kp == 'Synapse': 
			self.xOffset,self.yOffset = 40,9
			self.MIN_CLOUD_DIST = 6
			self.GROUP_CLOUD_DIST = 7
			self.MAX_CLOUD_DIST = 9
			self.MATCH_DIST=13
			self.PRIMARY_DIST=16
			self.ISLAND_DIST=19
		elif kp == 'VARIOUS':
			self.MIN_CLOUD_DIST = 7
			self.GROUP_CLOUD_DIST = 10
			self.MAX_CLOUD_DIST = 12
			self.MATCH_DIST=14
			self.PRIMARY_DIST=25
			self.ISLAND_DIST=20
		
		elif kp == 'Spliceosome':
			self.xOffset,self.yOffset = -50,6
			self.MATCH_DIST=8
			self.PRIMARY_DIST=12


		self.group_locs = {}
		self.primary_locations = {} 
		self.all_gene_locations = dd(list)  
		self.clustered_groups = [] 
	
		self.group_scores= ontology_scores
		self.group_key = ontology_groups 
		self.gene_key, self.obs = dd(list), dd(int) 
		for group in ontology_groups: 
			for gene in ontology_groups[group]: 
				self.gene_key[gene].append(group) 
				self.obs[gene] += 1 
				


#		self.group_key, self.gene_key, self.obs = ontology_groups, ontology_genes, ontology_counts
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
		self.FCFACT=50

 

                mX,mY =[p[0]*100 for p in scaler.fit_transform(np.array(gX).reshape(-1,1))] ,[p[0] * 100 for p in scaler.fit_transform(np.array(gY).reshape(-1,1))]


		group_pts = [(mX[i],mY[i]) for i in range(len(group_pts))]

		self.mds_locs = {ontology_names[i]:  [p[0],p[1]] for i,p in enumerate(group_pts)} 
		return 
		for i in range(len(ontology_names)):
			for j in range(i+1,len(ontology_names)): 
				oI,oJ = ontology_names[i],ontology_names[j] 
				gI,gJ = ontology_groups[oI], ontology_groups[oJ] 
				pI,pJ = group_pts[i], group_pts[j] 	
				m_dist = np.linalg.norm(np.array(pI)-np.array(pJ))
				s_num = len([g for g in gI if g in gJ])
				s_dist = len([g for g in gI if g in gJ]) / float(len(list(set(gI+gJ))))

				print kp,oI,oJ,len(gI),len(gJ),pI[0],pI[1],pJ[0],pJ[1],s_num,m_dist,s_dist








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
	


def group_class(ot_name): 


	#return "VARIOUS" 
	grps = [G_GROW,G_SPLICE,G_SYN,G_TUBE,G_DIFF]
	grp_names = ['Growth_and_Guidance','Spliceosome','Synapse','Microtubles','Differentiation']

	for a,b in zip(grps,grp_names):
		for g in a: print 'SAVE',g 
		if ot_name in a: 
			return b

	

	return "VARIOUS" 

	axon_list = ["activation_of_synapse_assembly","cerebral_cortex_GABAergic_interneuron_migration","dendrite_development","dendritic_spine_development","ERK_cascade","growth_cone","neurogenesis","neuron_guidance","node_of_Ranvier"]
	rc_list   = ["activation_of_calcium_ion-dependent_exocytosis","inhibition_of_calcium_ion_transmembrane_transporter_activity","membrane_depolarization_during_action_potential","somatodendritic_compartment","vesicle"]

	glial_list = ["chemokine-mediated_signaling_pathway","inhibition_of_amyloid_fibril_assembly","negative_regulation_of_hippocampal_neuron_apoptotic_process","neurofibrillary_tangle","regulation_of_neurogenesis"]
	if 'axon' in ot_name or ot_name in axon_list: return 'AXON' 

	if "channel" in ot_name or "recept" in ot_name or "synap" in ot_name or ot_name in rc_list: return "RECEPTORS_AND_CHANNELS"

	if 'microtub' in ot_name: return "MICROTUBLES" 

	if "glia" in ot_name or ot_name in glial_list: return "glial" 



	return "VARIOUS" 




def run_script(options,DEX_CUT=0.05,MID_CUT=0.01,HI_CUT=0.0001):

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
			for i in range(1,len(line)): 
				if line[i] == 'ES': continue 
				key[headers[i]][line[0]] = line[i]  
	IGNORE=dd(bool) 	
	if options.ignore != None: 
		for line in open(options.ignore): 
			line = line.split() 
			if len(line)>0:	IGNORE[line[0]] = True 


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


	topCands = dd(list) 
	kCands = [] 
	w = open('ontology_all.txt','w') 
	for k,members in ontology.ontology_key.items():

		ot_term,ot_name,ot_str = ontology.table[k]	
		
		if IGNORE[ot_name]: continue  
#		elif invalid_group(ot_name): continue 
#		if group_class(ot_name) not in ['VARIOUS']: continue 
#		if group_class(ot_name) not in ['Microtubles']: continue 
		gm = [ontology.gene_idx[m] for m in members if ontology.gene_idx[m] in cnts] 
		gm_pv = [dex.pv[g] for g in gm]

		PASS=False 
		topNames,topCnts,r_scrs,t_scrs,sig_fcs, avg_ps = [], [] ,[], [] ,[] , [] 



		for CUT in [DEX_CUT,MID_CUT,HI_CUT]: 

			sig_tops = [dex.top[g] for g in gm if dex.pv[g] < CUT] 
		
			if len(sig_tops) < 2: 
				PASS=False 
				break 
			else:
				PASS=True
			sig_cc = cc(sig_tops) 
			topName,topCnt = sorted(sig_cc.items(),key=lambda X: X[1],reverse=True)[0] 
			topNames.append(topName) 
			topCnts.append(topCnt)

			t_scr = sig_cc[topName]/float(len(gm)) 
			r_scr = sig_cc[topName]/float(sum(sig_cc.values()))

			t_scrs.append(t_scr) 
			r_scrs.append(r_scr) 

			sig_fc  =  [dex.fc[g] for g in gm if dex.pv[g] < CUT and topName == dex.top[g] ] 
			sig_non =  [dex.fc[g] for g in gm if dex.pv[g] < CUT and topName != dex.top[g] ] 

			if len(sig_non) > 0: 
				sig_fcs.append(np.mean(sig_fc)/(0.1+np.mean(sig_non)))
			else:
				sig_fcs.append(np.mean(sig_fc))
			sig_top_scr = [dex.pv[g] for g in gm if dex.top[g] == topName]	
			avg_ps.append(np.mean(sig_top_scr))

			minP = min(sig_top_scr) 
			



		if len(topNames) == 0: continue 

		elif PASS == True and len(list(set(topNames))) == 1:
			topName = topNames[0] 
		else:
			topName = 'MULTI'

		w.write('%-30s %20s %10.3e %10.3f %10.3f %10.5f %10.3f %10.3f %10s\n' % (topName,ot_name,minP,np.mean(r_scrs),np.mean(t_scr),avg_ps[0],np.mean(topCnts),np.mean(sig_fcs),len(gm)))

		
		
		topCands[group_class(ot_name)].append([minP,avg_ps[0],np.mean(r_scrs),np.mean(t_scrs),np.mean(topCnts),topName,k])





		




	ontology_groups = dd(list) 
	ontology_scores = {} 

	op_groups = sorted(topCands.keys())
#	op_groups = ['ES','EB','AB'] 
#	op_groups = ['EB','AB'] 
#	op_groups = ['EB'] 
	op =  OntologyPlot(options,dex,op_groups)
	



	for kp in op_groups:
		ontology_groups = dd(list) 
		ontology_scores = {} 
		tCands = sorted(topCands[kp]) 
		for i,(minP,avgP,topRate,topGlobal,topCnt,topName,k) in enumerate(tCands):
			gm = [ontology.gene_idx[m] for m in ontology.ontology_key[k] if ontology.gene_idx[m] in cnts]
			sig_genes = [g for g in gm if dex.pv[g] < DEX_CUT]
			sig_tops = [dex.top[g] for g in sig_genes] 		
			ot_term,ot_name,ot_str = ontology.table[k]	
			sig_cc = cc(sig_tops) 

			t_scr = sig_cc[topName]/float(len(gm)) 
			r_scr = sig_cc[topName]/float(sum(sig_cc.values()))
			f_top = [dex.fc[g] for g in gm if dex.top[g] == topName] 
			f_not = [1.0] + [dex.fc[g] for g in gm if dex.top[g] != topName] 
			if len(gm)>30: continue

			f_scr = sum(f_top)/sum(f_not) 
			ontology_groups[ot_name] = gm 
			ontology_scores[ot_name] = (f_scr,r_scr,t_scr,topName) 

			if len(ontology_groups.keys()) > 8: break 



		op.set_cluster_locations(kp,ontology_groups,ontology_scores) 
#		op.run_clustering() 
		op.run_clustering2() 
		op.add_group_labels(kp) 
	op.finish() 
	sys.exit() 



	

				















if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-i", "--ignore", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("--id", default = None, type='string', help="Category Id")
	parser.add_option("-c", "--cnts", default = None, type='string', help="Output Filename Prefix")
 	parser.add_option("--dupes", default = False, action='store_true', help="Output Filename Prefix")



	(options, args) = parser.parse_args()

	run_script(options)	













