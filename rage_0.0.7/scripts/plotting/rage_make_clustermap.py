#!/usr/bin/env python

import os
import sys
import matplotlib
import seaborn as sns 
#matplotlib.use('Agg')
import pylab as plt
import numpy as np
import math
import matplotlib.lines as mlines
from scipy.optimize import curve_fit

from collections import Counter as cc
from matplotlib import cm 
from collections import defaultdict as dd
from matplotlib.lines import Line2D as Line
from matplotlib import cm as cm
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Rectangle as Rect
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib._png import read_png
from PIL import Image
from math import exp
from math import sin
from math import pi
import matplotlib
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import StandardScaler
#import seaborn
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
from matplotlib.lines import Line2D as Line
import seaborn
# Generate the image
import copy
import sys
import os
import pickle
from collections import defaultdict as dd 
from collections import defaultdict as dd 
from datetime import datetime
from math import fabs
from math import log
import numpy as np 
import pylab 
import scipy.stats as stats
import os
import sys
import numpy as np
from scipy.interpolate import spline
from scipy.interpolate import spline
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.collections as mcoll
from matplotlib.collections import LineCollection
#!/usr/bin/env python
from matplotlib.patches import Circle as Circ
import matplotlib
#import seaborn
import matplotlib.pyplot as plt
from random import random 
# Generate the image
import copy
import sys
import os
import pickle
from collections import defaultdict as dd 
from collections import defaultdict as dd 
from collections import Counter as cc 
from datetime import datetime
from math import fabs
#!/usr/bin/env python
from math import log
import numpy as np 
import pylab 
import scipy.stats as scistats
import os
import sys
import math 
import matplotlib.pyplot as plt
from scipy.stats import chi2
from collections import defaultdict as dd 
from random import random
from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import FancyArrowPatch as Arrow
from scipy.stats import variation as cV
#from matplotlib.colors import ListedColormap, BoundaryNorm

matplotlib.rcParams['xtick.labelsize'] = 14



COLORS_1 = [ 'indigo', 'gold', 'hotpink', 'firebrick', 'indianred', 'sage', 'yellow', 'mistyrose', 'darkolivegreen', 'olive', 'darkseagreen', 'pink', 'tomato', 'lightcoral', 'orangered', 'navajowhite', 'lime', 'palegreen', 'darkslategrey'
, 'greenyellow', 'burlywood', 'seashell', 'mediumspringgreen', 'fuchsia', 'papayawhip', 'blanchedalmond', 'chartreuse', 'dimgray', 'black', 'peachpuff', 'springgreen', 'aquamarine', 'white', 'orange', 'lightsalmon', 'darkslategray', 'brown', 'ivory', 'dodgerblue', 'peru', 'darkgrey', 'lawngreen', 'chocolate', 'crimson', 'forestgreen', 'slateblue', 'lightseagreen', 'cyan', 'mintcream', 'silver', 'antiquewhite']

COLORS_2 = [ 'mediumorchid', 'skyblue', 'gray', 'darkturquoise', 'goldenrod', 'darkgreen', 'floralwhite', 'darkviolet', 'darkgray', 'moccasin', 'saddlebrown', 'grey', 'darkslateblue', 'lightskyblue', 'lightpink', 'mediumvioletred', 'slategrey', 'red', 'deeppink', 'limegreen', 'darkmagenta', 'palegoldenrod', 'plum', 'turquoise', 'lightgrey', 'lightgoldenrodyellow', 'darkgoldenrod', 'lavender', 'maroon', 'yellowgreen', 'sandybrown', 'thistle', 'violet', 'navy', 'magenta', 'dimgrey', 'tan', 'rosybrown', 'olivedrab', 'blue', 'lightblue', 'ghostwhite', 'honeydew', 'cornflowerblue', 'linen', 'darkblue', 'powderblue', 'seagreen', 'darkkhaki']

COLORS_3 = [ 'snow', 'sienna', 'mediumblue', 'royalblue', 'lightcyan', 'green', 'mediumpurple', 'midnightblue', 'cornsilk', 'paleturquoise', 'bisque', 'slategray', 'darkcyan', 'khaki', 'wheat', 'teal', 'darkorchid', 'deepskyblue', 'salmon', 'darkred', 'steelblue', 'palevioletred', 'lightslategray', 'aliceblue', 'lightslategrey', 'lightgreen', 'orchid', 'gainsboro', 'mediumseagreen', 'lightgray', 'mediumturquoise', 'darksage', 'lemonchiffon', 'cadetblue', 'lightyellow', 'lavenderblush', 'coral', 'purple', 'aqua', 'lightsage', 'whitesmoke', 'mediumslateblue', 'darkorange', 'mediumaquamarine', 'darksalmon', 'beige', 'blueviolet', 'azure', 'lightsteelblue', 'oldlace']

COMMON = ['red','blue','green','yellow','orange','purple','lime','cyan','k']
COLORS = COMMON+list(set([x for x in [c for c in COLORS_1+COLORS_2+COLORS_3] if x not in COMMON]))



























































class SS_COUNTS:
        def __init__(self,options): 
		self.options = options 
		self.genes, self.cnts, self.key = [], [], dd(lambda: {}) 


		self.multi_colors = {} 
		self.multi_colors['REG']={'CSP': 'blue','ISVZ':'purple','MZ': 'cyan', 'ES': 'lime','TP': 'crimson','HP': 'pink'} 

		self.multi_colors['LOC']= {'CP': 'blue','SP': 'purple','IZ': 'darkgreen','SVZ':'lightgreen','MZ':'cyan','IZ_SVZ': 'green'}

		self.multi_colors['CT'] = {'CR':'cyan'} #'MINI': 'brown','CR': 'cyan', 'MN': 'red','SN': 'blue', 'ITN': 'crimson', 'PYR': 'orange'}
		self.multi_colors['CTX'] = {'MEGA': 'lime','CR':'cyan'} #'MINI': 'brown','CR': 'cyan', 'MN': 'red','SN': 'blue', 'ITN': 'crimson', 'PYR': 'orange'}
		self.multi_colors['CTY'] = {'MEGA': 'lime','CR':'cyan'} #'MINI': 'brown','CR': 'cyan', 'MN': 'red','SN': 'blue', 'ITN': 'crimson', 'PYR': 'orange'}
		self.multi_colors['CT_ALT'] = {'MEGA': 'lime','MINI': 'brown','CR': 'cyan', 'MN': 'red','SN': 'blue', 'ITN': 'crimson', 'PYR': 'orange'}
		self.multi_colors['SURE'] = {'YES': 'green','NO': 'k'}
		self.multi_colors['SRC'] = {'AB': 'red','EB': 'blue','ES': 'orange'}

		self.multi_colors['LOC-CT'] = {'MZ-CR': 'cyan'} 
		self.multi_colors['REG-CT'] = {'MZ-CR': 'cyan'} 
		self.multi_colors['PBIO'] = {'UNIQ': 'snow'} 
		self.multi_colors['FS'] = {'NO': 'snow'} 
		self.multi_colors['HBB1'] = {'HBB_NO': 'snow'} 
		self.multi_colors['HBB2'] = {'HBB_NO': 'snow'} 
		self.multi_colors['FS'] = {'NO': 'snow'} 
		self.multi_colors['NBIO'] = {'UNIQ': 'snow'} 
		self.multi_colors['XBIO'] = {'UNIQ': 'snow'} 
		self.multi_colors['BAMPS'] = {'A2': 'green','A3': 'red'} 
		self.multi_colors['BARCODE'] = {'X': 'green'} 
		self.multi_colors['SPAN'] = {'AB': 'red','FT': 'blue'} 
		self.multi_colors['CT-BGW'] = {'MINI-R17': 'crimson'}
		self.multi_colors['CT-NBIO'] = {'MINI-B1D67': 'crimson'}
		self.multi_colors['CT-EI'] = {'MINI-UNSURE': 'snow'}
		self.multi_colors['ITN'] = {'PROJECTION': 'purple','UNK': 'white'} 
		self.multi_colors['PJN'] = {'INTERNEURON': 'olive','UNK': 'white'} 
		self.multi_colors['EI'] = {'EXCITATORY': 'lime','UNSURE': 'grey','INHIBITORY': 'red'}
		#self.multi_colors['BGW'] = {'UNIQ ': 'snow'} 
		self.multi_colors['BGW'] = {'R10': 'red','R12': 'orange','R14':'yellow','R15':'lightgreen','R16':'green','R17':'lightskyblue','R18':'blue','R19':'indigo','R20':'violet' } 

		self.multi_legend = dd(bool) 

		self.multi_color_index = dd(int) 


		self.gene_colors = dd(list) 
		self.SD = 60 
		self.scr_key = dd(list) 
		self.cnt_file, self.key_file = options.cnts, options.key 
		self.add_key() 


		self.mname = self.cnt_file.split('/')[-1]+'.'+options.key.split('/')[-1].split('.key')[0]+"_"+"-".join(self.IDS)+".cm"+str(self.SD)
#		if not options.scale: 
#			self.mname+='.UNSCALED'

		self.ids = [] 
		self.add_cnts() 
		self.make_plot() 
		self.finish() 


	def get_color(self,x,ID=None,SAMPLE=None):



		if x == 'NA': return 'snow' 
		elif x == 'UNK': return 'grey' 


		if ID != None: 
			if ID in self.multi_colors: 
				self.multi_legend[ID] = True 
				if x not in self.multi_colors[ID]:
					self.multi_colors[ID][x] = COLORS[self.multi_color_index[ID]]
					self.multi_color_index[ID]+=1
				return self.multi_colors[ID][x] 
			else:
				print ID
				sys.exit() 


		return 'snow'


	def add_key(self):

		if self.options.verbose: sys.stderr.write('\nAdding key...') 
	
		for line in open(self.key_file): 
			line = line.split()
			if line[0] == '---': headers = line
			else:
				for i in range(1,len(line)):
					self.key[headers[i]][line[0]] = line[i] 

		if self.options.id != None: 
			self.IDS = self.options.id.split(',') 
			self.color_ids = self.IDS 
		else: 			   
			self.IDS = ['SRC','REG','LOC','CT','CT_ALT',"SURE"] 
			self.color_ids =  ['SRC','SURE','REG','CT']

		 



	def score_ids(self,id_list,gene,my_cv): 


		TOPNAME = False 	
		for i,ID in enumerate(self.IDS): 

			vt = cc([x[i] for x in id_list])
			scores = sorted([(cv/float(self.exp_ids[ID][cx]),cx) for cx,cv in vt.items() if cv >10],reverse=True)

			if len(scores) == 0: 
				self.scr_key[gene].append((ID,'NA',1,1)) 
				continue  
			else:
	
				topScr,topName,nextScr,nextName = scores[0][0],scores[0][1],0.1,'None'
				if len(scores) > 1: nextScr,nextName = scores[1][0],scores[1][1] 

				self.wRes.write('%s %s %s %s %s | %s %s %s\n' % (gene,ID,topName,topScr,topScr/nextScr,nextName,nextScr,my_cv))
				self.scr_key[gene].append((ID,topName,round(topScr,3),round(topScr/nextScr,3)))

				if ID == 'CTX': 
					if topScr > 1.5: 
						TOPNAME = topName 

		return TOPNAME 


	def add_cnts(self): 
		
		DT = 'RAW'
		SCORE_IDS = True
		if self.options.verbose: sys.stderr.write('\nAdding counts...') 
		kk= 0 
		self.genes,self.gene_types,self.gene_colors, self.cnts = [] ,[], [], [] 
		scaler = MinMaxScaler(feature_range=(0,1))
		for line in open(self.cnt_file):
			line = line.split() 
			if line[0] == '---' or line[0] == 'Module': 
				#DT='RAW'
#				self.samples = line[1::][0:55] 
				self.samples = line[1::]
				ids = sorted([([self.key[ID][x] if x in self.key[ID] else 'MISSING' for ID in self.IDS],x,i)  for i,x in enumerate(self.samples)])


				ids = [x for x in ids if 'MISSING' not in x[0]]


						
				ordered_idxs = [idx[-1] for idx in ids] 
				ordered_samples = [idx[1] for idx in ids] 
				self.exp_ids = {} 	

				for i,ID in enumerate(self.IDS): 
					
					ordered_ids  = [idx[0][i] for idx in ids] 
					rc = cc(ordered_ids)  
					rr = {cx: cv/float(sum(rc.values())) for cx,cv in rc.items()} 
					self.exp_ids[ID] = {cx: cv*self.SD for cx,cv in rr.items()} 



				self.wRes=open(self.mname+'.res','w')
 
			else:
				gene,vals,gene_type,gene_nick,scores = line[0],[float(x) for x in line[1::]],'?',line[0].split(';')[-1],[]
				vals = [vals[j] for j in ordered_idxs]
				
				if sum(vals) == 0: continue 


				if not self.options.raw: 
					try: vals = [log(z+1,2) for z in vals]
					except ValueError: DT='RESID'
				

#				kk+=1
#				if kk > 100: break

			
				topS = None 
				my_cv = cV(vals) 
				if self.options.scale: vals = [z[0] for z in scaler.fit_transform(np.array(vals).reshape(-1,1))]
				if SCORE_IDS: 

					id_vals = [x[1] for x in sorted([(vals[j],ids[j][0]) for j in range(len(vals))],reverse=True)[0:self.SD]]
					topS = self.score_ids(id_vals,gene,my_cv) 

				if topS != False:

					if topS == 'MEGA': 		self.gene_colors.append('lime') 
					elif topS == 'CR': 		self.gene_colors.append('cyan')
					elif topS == 'MINI-X' or topS == 'MINI_X':   self.gene_colors.append('orange') 
					elif topS == 'MINI-Y':	   self.gene_colors.append('purple') 
					else: 		   		self.gene_colors.append('white') 
				else:
					self.gene_colors.append('grey') 



				val_mean, val_std = np.mean(vals), np.std(vals) 
				centered_vals = [(v-val_mean) / val_std for v in vals]


				if self.options.method in ['ward','centroid']:
					self.cnts.append(centered_vals)
			
				else:
					self.cnts.append(vals) 			


				self.genes.append(gene) 
				self.gene_types.append(gene_type) 
				#if gene_type == 'EB': self.gene_colors.append('darkblue') 
				#else: 		      self.gene_colors.append('darkred') 
			

		self.wRes.close() 

		self.ids.extend(ids) 









	def pull_clusters(self,c_data,TYPE='GENES'):



		dn = c_data.dendrogram
		Z  = c_data.linkage

		if TYPE == 'GENES':	g_swap = {j: self.genes[r] for j,r in enumerate(c_data.reordered_ind)}
		else:			g_swap = {j: self.ids[r][1] for j,r in enumerate(c_data.reordered_ind)}
		
		cluster_list = [] 
		for aa,bb,cc in zip(Z,dn['dcoord'],dn['icoord']): 	
			if cc[0] != cc[1] or cc[2] != cc[3]: 
				print 'wtf'
				sys.exit() 
			span = cc[1],cc[2]
			height,left_height,right_height = bb[2],bb[0],bb[3]
			length = aa[-1] 

			cluster_list.append((height,length,(left_height,right_height),span)) 
		cluster_list.sort(reverse=True)
		height,length,(left_height,right_height),span = cluster_list[0] 
		current_level = [[((span[0],'L'),(span[1],'R'))]]
		cluster_list[0] = False 

		next_level = [] 
		g_clusters = {} 

		while True: 
			if len(next_level) > 0: break
			for sp in current_level[-1]:
				for (sPt,sDir) in sp:
					for j in range(len(cluster_list)): 
						if cluster_list[j] == False: continue 	
						ht,length,(left_height,right_height),jSpan = cluster_list[j] 
						jMid = np.mean(jSpan)
						if jMid == sPt:
							next_level.append(((jSpan[0],sDir+'L'),(jSpan[1],sDir+'R')))		
							cluster_list[j] = False						

			if len(next_level) == 0: break 
			current_level.append([nx for nx in next_level])
			next_level = [] 

			
		current_rev = [a for b in current_level[-1::-1] for a in b]
		cluster_key = {} 
		for cR in current_rev:
			for cLoc,cName in cR:
				cIdx, cFloat= int(cLoc/10.0 - 0.5),(cLoc/10.0 - 0.5)
				if float(cIdx) == cFloat: 
					g_name = g_swap[cIdx] 		
					if g_name not in cluster_key:
						if TYPE == 'GENES': cName = "".join(['U' if x == 'L' else 'D' for x in cName])
						g_clusters[g_name] = cName


		return g_clusters






	def score_clusters(self,s_data,ID,TYPE='SAMPLES',LEN=4):




		if TYPE == 'SAMPLES':	
			xkey = self.sample_cluster_ids
			w_scrs=open(self.mname+'-score-'+ID+'.sample_clusters','w')
		else:			
			xkey = self.gene_cluster_ids 
			w_scrs=open(self.mname+'-score-'+ID+'.gene_clusters','w')

		save_key = {} 

	


		cLens =  sorted([len(x) for x in list(set(xkey.values()))])
		cMin,cMax = cLens[0],cLens[-1]

		grp_key = dd(lambda: dd(int))
		member_key = dd(list) 
		for i in range(1,cMax+1):
			for sd in s_data:

				if sd[0] in xkey:				
					s_id,s_type,sk = sd[0],sd[1], xkey[sd[0]][0:i]
				else:
					s_id,s_type,sk = sd[0],sd[1], 'NNNNNNN'


				if i == LEN:
					save_key[s_id] = sk 


				if len(sk) < i: continue 

				grp_key[sk][s_type] += 1 
				member_key[sk].append(s_id)
					


		for sk,sNames in member_key.items():

			if len(sNames) <3: continue 
			xtotal = sum(grp_key[sk].values()) 
			tupes = [(x[0],x[1]/float(xtotal)) for x in sorted(grp_key[sk].items(),reverse=True,key=lambda X: X[1])]
			
			ns = ",".join([sn.split(';')[-1] for sn in sNames])
			w_scrs.write('%-20s %20s %6d %10s %8.3f' % (sk,ID,xtotal,tupes[0][0],tupes[0][1]))	
			if len(tupes) > 1: w_scrs.write(' %6s %5.3f' % (tupes[1][0],tupes[1][1]))
			if len(tupes) > 2: w_scrs.write(' %6s %5.3f' % (tupes[2][0],tupes[2][1]))


			w_scrs.write(' | %50s\n' % (ns))

		w_scrs.close() 
		return save_key
			

			



	def make_plot(self):

		
		if self.options.verbose: sys.stderr.write('\nMaking Plot...')

		self.col_colors = [] 



		for ID in self.color_ids: 
			iLoc = self.IDS.index(ID) 
			ID_COLORS = [self.get_color(x[0][iLoc],ID,x[1]) for x in self.ids]
			self.col_colors.append(ID_COLORS) 


		#print self.col_colors
 
		#col_colors = ['blue' if x[0:2] == 'EB' else 'red' for x in samples] 
		#g = sns.clustermap(cnts,col_cluster=False,row_colors=gene_colors, col_colors=col_colors)
		#g = sns.clustermap(cnts,method='ward',col_cluster=False, row_cluster=False,row_colors=gene_colors, col_colors=col_colors)
#		g = sns.clustermap(cnts,method='ward',col_cluster=True, row_cluster=True,row_colors=gene_colors, col_colors=col_colors)
		#g = sns.clustermap(cnts,method='centroid',col_cluster=True, row_cluster=True,row_colors=gene_colors, col_colors=col_colors)
		#g = sns.clustermap(cnts,method='weighted',col_cluster=True, row_cluster=True,row_colors=gene_colors, col_colors=col_colors)
		#g = sns.clustermap(cnts,cmap=camp,method='ward',metric='correlation',col_cluster=True, row_cluster=True,row_colors=gene_colors, col_colors=col_colors)
		#g = sns.clustermap(cnts,cmap=camp,method='ward',metric='correlation',col_cluster=False, row_cluster=True,row_colors=gene_colors, col_colors=col_colors)
#		g = sns.clustermap(cnts,cmap=camp,method='ward',metric='correlation',col_cluster=True, row_cluster=True,row_colors=gene_colors, col_colors=col_colors)
#		g = sns.clustermap(cnts,cmap=camp,method='weighted',metric='correlation',col_cluster=True, row_cluster=True,row_colors=gene_colors, col_colors=col_colors)


		METHOD = self.options.method 		

#		METHOD = 'weighted'
#		METHOD = 'centroid'
	
		
		OPSTR, ROBUST, CENTER, SCALE = '',False, False, False 
		if 'robust' in self.options.settings.split(','): 
			ROBUST = True
			OPSTR += 'R'
		if 'center' in self.options.settings.split(','): 
			CENTER = True 
			OPSTR += 'C'
		elif 'scale'  in self.options.settings.split(','): 
			SCALE = True
			OPSTR += 'S'

		if len(OPSTR) > 0:	self.mname += '.'+METHOD+'-'+OPSTR
		else:			self.mname += '.'+METHOD

#		self.gene_colors = [COLORS[x] for x in range(len(self.cnts))]

	


		GC,CC = self.gene_colors, self.col_colors
		if METHOD in ['ward','centroid']:
			if CENTER:	g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=METHOD,metric='correlation',z_score=1,   robust=ROBUST,col_cluster=True, row_cluster=True,row_colors=GC, col_colors=CC, zorder=5)
			elif SCALE:	g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=METHOD,metric='correlation',standard_scale=1,   robust=ROBUST,col_cluster=True, row_cluster=True,row_colors=GC, col_colors=CC, zorder=5)
			else:	
				g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=METHOD,metric='correlation',robust=ROBUST,col_cluster=True, row_cluster=True,row_colors=GC, col_colors=CC, zorder=5)
		elif METHOD in ['weighted','single']:
			if CENTER:	g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=METHOD,z_score=1,   robust=ROBUST,col_cluster=True, row_cluster=True,row_colors=GC, col_colors=CC, zorder=5)
			elif SCALE:	g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=METHOD,standard_scale=1,   robust=ROBUST,col_cluster=True, row_cluster=True,row_colors=GC, col_colors=CC, zorder=5)
			else:		g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=METHOD,robust=ROBUST,col_cluster=True, row_cluster=True,row_colors=GC, col_colors=CC, zorder=5)
				
		else:
			g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=method,robust=True,col_cluster=True, row_cluster=True,row_colors=self.gene_colors, col_colors=self.col_colors, zorder=5,z_score=0)
			#g = sns.clustermap(self.cnts,cmap=self.options.cmap,method=method,robust=True,col_cluster=True, row_cluster=True,row_colors=self.gene_colors, col_colors=self.col_colors, zorder=5,standard_scale=1)
		plt.setp(g.ax_heatmap.set_yticklabels([]), rotation=0)
		plt.setp(g.ax_heatmap.set_xticklabels([]), rotation=0)

#		self.gk = {self.genes[r]: r for r in g.dendrogram_row.reordered_ind}
#		self.sk = {self.samples[r]: r for r in g.dendrogram_col.reordered_ind}




		self.gene_cluster_ids = self.pull_clusters(g.dendrogram_row)
		self.sample_cluster_ids = self.pull_clusters(g.dendrogram_col,TYPE='SAMPLES')



		g_full_order = [(self.genes[r],self.scr_key[self.genes[r]]) for r in g.dendrogram_row.reordered_ind]
		s_full_order = [self.ids[r] for r in g.dendrogram_col.reordered_ind]
		
		self.gene_range = 25	
		self.sample_range = 8

		if self.gene_range > len(self.cnts): 	  self.gene_range = int(len(self.cnts)/2.0) - 1
		if self.sample_range > len(self.samples): self.sample_range = int(len(self.samples)/2.0) - 1
		
#		self.wNay=open(self.mname+'-'+str(self.sample_range)+'.sample_groups','w')
#		self.gNay=open(self.mname+'-'+str(self.gene_range)+'.gene_groups','w')
		xMin,xMax = g.ax_heatmap.get_xlim()
		yMin,yMax = g.ax_heatmap.get_ylim()


		xLeft = xMin - (xMax*0.15) 

		if self.options.verbose: sys.stderr.write('\nChecking Neighbors...')

		S_JUMP, G_JUMP = 0,0 
		if len(s_full_order) > 50: S_JUMP = int(len(s_full_order) / 10.0)
		if len(g_full_order) > 50: G_JUMP = int(len(g_full_order) / 10.0)


		spec_samps, spec_genes,IN = [],[] ,False 
		for i,ID in enumerate(self.IDS): 
			if ID in ['SURE']: continue 

			s_order = [(s[0][i],s[1],s[2]) for s in s_full_order]
		        g_order = [(gi[0],gi[1][i]) for gi in g_full_order] 
		        g_order = [(gi[0],gi[1][i]) for gi in g_full_order] 
#		        g_order = [(gi[0],gi[1][i]) for gi in g_full_order] 

		
			s_data = [(s[1],s[0][i],2,2) for s in s_full_order]
		       	g_data = [(gi[0],gi[1][i][1],gi[1][i][2],gi[1][i][3]) for gi in g_full_order] 
			
			sample_choice = self.score_clusters(s_data,ID,TYPE='SAMPLES') 	
			gene_choice   = self.score_clusters(g_data,ID,TYPE='GENES') 
			continue
			for j in range(len(g_data)):
				try:
					if G_JUMP == 0:
						g.ax_heatmap.text(xLeft-3,len(self.cnts)-j-0.5,gene_choice[g_data[j][0]],color='purple',fontsize=14,zorder=10,verticalalignment='center') 
						g.ax_heatmap.text(xLeft,len(self.cnts)-j-0.5,self.gene_cluster_ids[g_data[j][0]],color='purple',fontsize=9,zorder=10,verticalalignment='center') 
					elif j % G_JUMP == 0: 
						g.ax_heatmap.text(xLeft,len(self.cnts)-j-0.5,gene_choice[g_data[j][0]],color='purple',fontsize=14,zorder=10,verticalalignment='center') 

				except KeyError: 
					continue 

			for j in range(len(s_data)):
				try:
					if S_JUMP == 0:	
						g.ax_heatmap.text(j+0.5,yMax*1.05,self.sample_cluster_ids[s_data[j][0]],color='purple',fontsize=12,rotation=80,zorder=10,verticalalignment='center',horizontalalignment='center') 
						g.ax_heatmap.text(j+0.5,yMax*1.2,sample_choice[s_data[j][0]],color='purple',fontsize=10,rotation=80,zorder=10,verticalalignment='center',horizontalalignment='center') 
					elif  j % S_JUMP == 0: 	
						g.ax_heatmap.text(j+0.5,yMax*1.2,sample_choice[s_data[j][0]],color='purple',fontsize=10,rotation=80,zorder=10,verticalalignment='center',horizontalalignment='center') 
				except KeyError:
					continue 
			continue

			self.check_neighbors(s_order,ID) 
			self.check_gene_neighbors(g_order,ID) 


		
			
		
#		self.wNay.close() 


	def finish(self): 

		if self.options.verbose: sys.stderr.write('\nAdding Legend...')
		plt.axis('off')
		items,labels=[],[] 
		leg_items,leg_labels,leg_ids = [],[],[] 
		for ID in self.multi_legend: 
			items.append([])
			labels.append([])  


			for a,b in self.multi_colors[ID].items(): 
				items[-1].append(Rect((0,0),1,1,fc=b))
				labels[-1].append(a) 
			leg_ids.append(ID) 
			
		maxLen= max([len(x) for x in labels])
	
		leg_items,leg_labels = [[] for x in range(maxLen)],[[] for x in range(maxLen)] 



		for I,L in zip(items,labels): 

		
			while len(I) < maxLen:
				I.append(Rect((0,0),1,1,fc='white',ec='white'))
				L.append('')

			

			for k in range(len(L)): 
				leg_items[k].append(I[k])
				leg_labels[k].append(L[k])

		leg_items = [[Rect((0,0),1,1,fc='white',ec='white') for x in leg_ids]]+leg_items
		leg_labels = [leg_ids]+leg_labels


		ncols = len(leg_ids) 

		ncols = maxLen+1

		leg_items = [a for b in leg_items for a in b]
		leg_labels = [a for b in leg_labels for a in b]
			

		#hp,cs = 0.1, 0.1
		hp,cs = 0.13, 0.3
		plt.legend(leg_items,leg_labels,ncol=ncols,bbox_to_anchor=(1.2,1.7),loc='upper left',handletextpad=hp,columnspacing=cs,fontsize=12) 
		plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.82) # ,wspace=0.25,hspace=1.05)


		plt.savefig(self.mname+'.png',dpi=500) 

		if options.show:
			plt.show() 

		if self.options.verbose: sys.stderr.write('\n') 
		



	def check_neighbors(self,s_order,ID):


		FULL_RANGE, HALF_RANGE,CHECK_RANGE = self.sample_range,self.sample_range/2,(self.sample_range/2)-1
		for i in range(len(s_order)): 

			if i < HALF_RANGE: 
				left=0 
				right=i+(FULL_RANGE-i)+1 

			

			elif i + HALF_RANGE > len(s_order):
				left = i - (FULL_RANGE - (len(s_order) - i) )-1 
				right = len(s_order) 

			else: 
				left = i - HALF_RANGE 
				right = 1+i+HALF_RANGE


			prev = [p[0] for p in s_order[left:i]]
			post = [p[0] for p in s_order[i+1:right]] 

			preCC,postCC,bothCC = cc(prev), cc(post) , cc(prev+post) 

			
			preTups = sorted([(x,float(y)/sum(preCC.values())) for x,y in preCC.items()],reverse=True,key=lambda X: X[1]) 
			postTups = sorted([(x,float(y)/sum(postCC.values())) for x,y in postCC.items()],reverse=True,key=lambda X: X[1]) 
			bothTups = sorted([(x,float(y)/sum(bothCC.values())) for x,y in bothCC.items()],reverse=True,key=lambda X: X[1]) 


			topBoth,nextBoth = bothTups[0],('NA',0.01)
			if len(bothTups)>1: nextBoth = bothTups[1] 



			self.wNay.write('%s %s %s %s ' % (s_order[i][1],i,ID,s_order[i][0]))
			self.wNay.write('%s %4.4f %3.2f %s %4.4f | ' % (topBoth[0],topBoth[1],topBoth[1]/nextBoth[1],nextBoth[0],nextBoth[1]))

			topPrev,topPost = ('NA',0.01),('NA',0.01) 
			if len(prev) > CHECK_RANGE: topPrev = preTups[0]
			if len(post) > CHECK_RANGE: topPost = postTups[0]
			allTops = sorted([topPrev,topPost,topBoth],reverse=True,key=lambda X: X[1]) 
			self.wNay.write('%s %4.4f \n' % (allTops[0][0],allTops[0][1]))
			
			









	def check_gene_neighbors(self,s_order,ID):


		FULL_RANGE, HALF_RANGE,CHECK_RANGE = self.gene_range,self.gene_range/2,(self.gene_range/2)-1

		for i in range(len(s_order)): 

			if i < HALF_RANGE: 
				left=0 
				right=i+(FULL_RANGE-i)+1 

			

			elif i + HALF_RANGE > len(s_order):
				left = i - (FULL_RANGE - (len(s_order) - i) )-1 
				right = len(s_order) 

			else: 
				left = i - HALF_RANGE 
				right = 1+i+HALF_RANGE




			left_key,right_key,both_key = dd(list),dd(list),dd(list)
			nay_key = dd(list) 
			for p in s_order[left:i]:	
				left_key[p[1][1]].append(p[1][2]) 
				both_key[p[1][1]].append(p[1][2]) 
			for p in s_order[i+1:right]:	
				right_key[p[1][1]].append(p[1][2]) 
				both_key[p[1][1]].append(p[1][2]) 



			left_list  = sorted([(len(vals),len(vals)/float((i-left)),np.mean(vals),p) for p,vals in left_key.items()],reverse=True)
			right_list =  sorted([(len(vals),len(vals)/(right-i-1.0),np.mean(vals),p) for p,vals in right_key.items()],reverse=True)

			

			both_list =  sorted([(len(vals),len(vals)/((right-i-1.0)+(i-left)),np.mean(vals),p) for p,vals in both_key.items()],reverse=True)
			if len(both_list) == 0: continue 

			
			self.gNay.write('%s %d %s %s %3.3f %3.3f |' % (s_order[i][0],i,s_order[i][1][0],s_order[i][1][1],s_order[i][1][2],s_order[i][1][3]))
			

			self.gNay.write(' %s %3.3f %3.3f | ' % (both_list[0][-1],both_list[0][1],both_list[0][2]))

			xCands = [both_list[0]]
			if  left-i> CHECK_RANGE: xCands.append(left_list[0]) 
			if (right-i-1)> CHECK_RANGE: xCands.append(right_list[0]) 
			xTop =  sorted(xCands,reverse=True)[0] 

			self.gNay.write(' %s %3.3f %3.3f \n' % (xTop[-1],xTop[1],xTop[2]))

				

		
			


			









































































if __name__ == '__main__':

	import sys
	import os
	from optparse import OptionParser

	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)



	parser.add_option('--key',  dest= "key", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--cmap',  dest= "cmap", type = 'str' , default = 'seismic', help = "horizontal data")
	parser.add_option('--method',  dest= "method", type = 'str' , default = 'ward', help = "horizontal data")
	parser.add_option('--settings',  dest= "settings", type = 'str' , default = '', help = "horizontal data")
	parser.add_option('--id',  dest= "id", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('-c','--cnts',  dest= "cnts", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('-b','--bcnts',  dest= "bcnts", type = 'str' , default = None, help = "horizontal data")
	parser.add_option('--scale',  dest= "scale", action = 'store_true', default = False, help = "horizontal data")
	parser.add_option('--show',  dest= "show", action = 'store_true', default = False, help = "horizontal data")
	parser.add_option('--verbose',  dest= "verbose", action = 'store_true', default = False, help = "horizontal data")
	parser.add_option('--raw',  dest= "raw", action = 'store_true', default = False, help = "horizontal data")


	##############################################################################
	######################   STEP 1: LOAD FILE  ################################## 
	##############################################################################

	# colors = ['black','darkblue','limegreen','olive','purple','pink','Chartreuse','brown','cyan','black','Brown','DarkGray']
	(options, args) = parser.parse_args()

	color_key = {'singleCell~AB':'red',  'singleCell~NEW': 'gold',  'singleCell~EB': 'blue',  'brainSpan~EB': 'cyan',  'brainSpan~AB': 'purple'}
	color_key = {'SS': {'AB': 'red', 'EB': 'blue'}, 'BS': {'AB': 'magenta', 'EB': 'cyan'},'EI': {'AB': 'red', 'EB': 'blue'}}

	xLen, yLen = 2,2 
	xLoc, yLoc = 1,1 

	CAMP='YlOrBr'
	CAMP='hot'
	CAMP='seismic'


	s_data = SS_COUNTS(options) 
#	if options.cnts: 
#		s_cnts = SS_COUNTS(options.cnts,options,camp=CAMP) 






	sys.exit() 


	plt.subplots_adjust(left=0.07,right=0.95,wspace=0.05,hspace=-0.05,top=0.92,bottom=0.05)
#	plt.savefig('DEX_FIG.png',dpi=200) 
#	plt.savefig('DEX_FIG.pdf',dpi=200) 
	plt.show() 
#	plt.show() 
	sys.exit() 
    	 














