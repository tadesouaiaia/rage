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



print 



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



class CellKey:
        def __init__(self,key,options):

		self.options,self.key = options, key 



	def identify(self,s):
		self.name = s 
		self.marker = 'o'

		cs,cl,clx,ct,cf = self.key['CS'][s],self.key['CL'][s],self.key['CLX'][s],self.key['CT'][s],self.key['CF'][s]
		self.obs,self.tot,self.gw,self.bio = self.key['R_OBS'][s],self.key['R_TOT'][s],self.key['GW'][s],self.key['BIO'][s]	



		fs,maxspikes = self.key['FS'][s],self.key['MAXSPIKES'][s]
		marker = 'o'

		if fs != 'NA': self.spikes = int(maxspikes) 
		else:          self.spikes, self.fs, self.fire_color, = None, 'NA','white' 

		if fs == 'REPI' and self.spikes > 30: 		    self.fs, self.fire_color = 'SUSTAINED', 'lime'
		elif fs in ['REPI','DUBS','SINGLE']: 		    self.fs, self.fire_color = 'ACTIVE', 'blue'
		elif fs != 'NA':				    self.fs, self.fire_color = 'ABORTIVE','gray'
	


				
		types = sorted(list(set([x for x in [cs,cl,clx,ct] if x not in ['NA','UNK']])))
		

		if len(types) == 0: self.type,self.type2,self.marker, self.type_color = 'NA','NA','o','white'
		elif 'CSP' in types:
			self.type, self.type2,self.marker,self.type_color = 'CANONICAL','CSP','^','yellow'
		elif 'MEGA' in types:	
			self.type, self.type2,self.marker,self.type_color = 'NOVEL','ISVZ','v','purple'
		elif 'CR' in types: 
			self.type, self.type2, self.marker, self.type_color = 'CR','MZ','D','cyan' 	
		elif 'ISVZ' in types: 
			self.type, self.type2,self.marker,self.type_color = 'CANONICAL','ISVZ','v','yellow'
		
		elif 'MZ' in types: 
			self.type, self.type2, self.marker, self.type_color = 'MZ','MZ','D','orange' 	
		elif 'MINI' in types:
			self.type, self.type2,self.marker,self.type_color = 'CANONICAL','CANONICAL','o','tan'
		else:
			
			print types 
			sys.exit() 

			
		if self.type == 'MZ': 
			if self.obs < 2000 and self.obs < 4000: self.type,self.type_color = 'NA', 'white'
			else: 		    self.type, self.type_color = 'CANONICAL', 'yellow'	

			



		return self















		


class ScatterPlot:
        def __init__(self,options,xLen=3,yLen=4,key={}):

                self.options = options
                self.VERBOSE = True
                
		self.xOffset = 0 
        
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
                self.xLen, self.yLen = xLen,yLen 
                self.xLoc, self.yLoc = 0,0
		self.color_key = {'NA': 'k', 'UNK': 'grey','+': 'cyan','-': 'olive','16': 'lightblue','17': 'lime','18': 'red','14':'yellow','19':'purple'} 

		self.color_key['SVZ'] = 'crimson'
		self.color_key['SP'] = 'firebrick'
		self.color_key['CP'] = 'gold'
		self.color_key['MZ'] = 'olive'
		self.cnt_index = 0 
		self.axes= [] 		
		self.color_offset = 0


		MZ_Mark =  Line2D([], [], color='white', marker='D',markeredgecolor='k',markeredgewidth=0.5,linestyle='None',markersize=10, label='Marginal Zone')
		UP_Mark =  Line2D([], [], color='white', marker='^',markeredgecolor='k',markeredgewidth=0.5,linestyle='None',markersize=10, label='Upper Layers (SP/CP)')
		DOWN_Mark=Line2D([], [], color='white', marker='v',markeredgecolor='k',markeredgewidth=0.5,linestyle='None',markersize=10, label='Lower Layers (SVZ/IZ)')

		self.loc_labs = [MZ_Mark,UP_Mark,DOWN_Mark]









	def finish(self,out='figure',title=None,xLen=2,yLen=3,kind='RAW'):
		out_name = out+'.png'
		#if title != None:
		#	plt.suptitle(title,fontsize=25,fontweight='bold') 
		#else:   
		#	plt.suptitle(out,fontsize=25,fontweight='bold') 
 

  		plt.subplots_adjust(left=0.03, bottom=0.02, right=0.98, top=0.96,wspace=0.10,hspace=0.03)		
		plt.savefig(out_name,dpi=250) 
		plt.show() 
		plt.clf() 

                self.xLen, self.yLen = xLen,yLen 
		self.xLoc,self.yLoc = 0, 0 
		self.xOffset = 0 


	def get_color(self,x): 

		if len(x.split(',')) == 2: 
			a,b = x.split(',') 
			if a not in self.color_key: 
				self.color_key[a] = COLORS[self.color_offset]
				self.color_offset += 1
				
			if b not in self.color_key: 
				self.color_key[b] = COLORS[self.color_offset]
				self.color_offset += 1

			return self.color_key[a],self.color_key[b] 
	
		else:
			try: return self.color_key[x] 

			except KeyError: 
				self.color_key[x] = COLORS[self.color_offset]
				self.color_offset +=1
			return self.color_key[x] 








	

	def plot_custom_pts(self,pts,key,cnts=None,cnt_idx =0,genes=[],neg_genes=[],color=None,VERSION='LOC',TITLE='HI',BAR_MSG=None,AXES=[]):
		if len(AXES) == 0: 
			self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))
		else:
			xl,yl,rs,cs = AXES 
			self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (xl,yl), rowspan = rs, colspan = cs))



		cell_key = CellKey(key,self.options)  
		p_data,p_genes,p_names,label_key = [],[],[],{}

		if len(genes)>0:	
			TITLE = ",".join(genes) 	
			if len(neg_genes)>0: 
				TITLE +='\n'+",".join(neg_genes)
				genes+=neg_genes
			p_genes = [[] for p in genes]


		for s in pts:
			x,y = pts[s]
			try:  
				cell = cell_key.identify(s) 
			except KeyError:
				continue


			if cell.type == 'NA': continue 
			if VERSION == 'FIRING':	
				clr,clr_label = cell.fire_color, cell.fs 
			else:			
				
				clr,clr_label = cell.type_color, cell.type
				if cell.type == 'MINI': cell.type = 'canonical'
				#print s,x,y,cell.type, cell.type2, cell.marker
						
				#if cell.type == 'MZ' and cell.type2 == 'MZ': continue 

			

			if clr_label not in label_key: 	
				self.color_key[clr_label] = clr 

				label_key[clr_label] = Rect((0,0),1,1,fc=clr)	


			if cnts == None:
				if cell.marker == 'o': continue 
				if VERSION == 'FIRING':		self.axes[-1].scatter(x,y,marker=cell.marker,c=cell.fire_color,s=70,alpha=1) 
				else:				self.axes[-1].scatter(x,y,marker=cell.marker,c=clr,s=70,alpha=1) 



			else: 
				for j in range(len(genes)): 
					p_genes[j].append(log(cnts[genes[j]][s]+1,2)) 
					#print j,p_genes[j],cnts[genes[j]]



				p_data.append([np.array([log(cnts[g][s]+1,2) for g in genes]),x,y,clr,cell])


		if cnts != None:

			#print len(p_data) 
			#print len(p_data[0])
			#print p_data[0]  
			#sys.exit() 



			if len(genes) > 0: 
				scaler = MinMaxScaler() 
				t_data = scaler.fit_transform(np.array([p[0] for p in p_data]))
				for j in range(len(p_data)): 
					pos_cont = sum([t_data[j][x] for x in range(len(genes)) if genes[x] not in neg_genes])
					neg_cont = sum([t_data[j][x] for x in range(len(genes)) if genes[x] in neg_genes])


					#for x in range(len(genes)):
					#	print x,genes[x],t_data[j][x],p_data[j][0]


					p_data[j][0] = pos_cont - neg_cont 

	
				
			p_data.sort(reverse=True) 
			p_data.sort() #reverse=True) 
			X = np.array([p[1] for p in p_data])
			Y = np.array([p[2] for p in p_data])
			I = np.array([p[0] for p in p_data]) 


			#print iI


			#z = gaussian_kde(I)(I) 



			Xs = scaler.fit_transform(X.reshape(-1,1))
			Ys = scaler.fit_transform(Y.reshape(-1,1))


			scatterplot = self.axes[-1].scatter(X,Y,c=I,cmap='seismic',s=40,edgecolor='') 
			plt.colorbar(scatterplot,ax=self.axes[-1], orientation='horizontal', fraction = 0.04, pad=0.1, ticks=[]) 
			self.axes[-1].scatter(X,Y,c=I,cmap='seismic',s=40,edgecolor='') 

			if BAR_MSG == 'RELN':
				self.axes[-1].text(0.6,-7.5,'      CR Markers\n(RELN,CXCR4,PCP4)',fontsize=11,verticalalignment='bottom') 
			elif BAR_MSG[0:3] == 'CHL':
				self.axes[-1].text(0.3,-7.5,'   Choride Ratio\nSLC12A2:SLC12A5',fontsize=11,verticalalignment='bottom') 
			elif BAR_MSG[0:3] == 'GLU':
				self.axes[-1].text(0.0,-7.5,'      Glutamatergic:GABA \n(SLC17A7,SLC1A2,GRIK3:ABAT,GAD)',fontsize=9,verticalalignment='bottom') 
			elif BAR_MSG[0:3] == 'VEN':
				self.axes[-1].text(0.3,-7.5,'    VEN Markers\n(FEZF2,BCL11B)',fontsize=11,verticalalignment='bottom') 
			else:
				self.axes[-1].text(0.3,-7.5,BAR_MSG,fontsize=11,verticalalignment='bottom') 

			#self.axes[-1].scatter(X,Y,c=z,cmap='jet',s=50,edgecolor='') 		
			#self.axes[-1].scatter(X,Y,c=C,s=100,alpha=0.5) 
			#self.axes[-1].scatter(X,Y,c=I,cmap='hot',s=50,edgecolor='') 
			#labels = label_key.keys()
			#items = [label_key[x] for x in labels]		
			#nc,fs,cs,hp,bb1,bb2 = 2,12 ,0.5,0.5,-0.01,1.1
			#leg = self.axes[-1].legend(items,labels,title=TITLE,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')
	
		else: 
			
			nc,fs,cs,hp,bb1,bb2 = 2,9.0,0.62,0.62,0.465,1.0

			PROX = [Line2D([],[],color=self.color_key[k],marker='s',mec='k',markeredgewidth=0.2,linestyle='None',markersize=10,label=k) for k in label_key if k != 'NA']
			leg = self.axes[-1].legend(handles=PROX+self.loc_labs,title=TITLE,handletextpad=hp,columnspacing=cs,fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')
			
			



		self.axes[-1].axis('off') 


		self.yLoc += 1
		if self.yLoc == self.yLen: 
			self.yLoc =0 
			self.xLoc +=1































def run_script(data_file,options):
	k=0
	key = dd(lambda: {}) 

	prefix = data_file.split('.pts')[0] 	
	anno_key = dd(list) 
	genes,cnts,annos = [],[],[]
	pt_key = {} 
	for line in open(data_file): 
		line = line.split() 
		if line[0] == '---': 
 			continue 
		else:
			pt_key[line[0]] = [float(x) for x in line[1:3]]

	samples = pt_key.keys() 

	


	for line in open(options.key):
		line = line.split() 
		if line[0] == '---': 
 			headers = line  
		else:
			s = line[0] 
			if s not in samples: continue 
			for i in range(1,len(line)): 
				try: 
					x = float(line[i])
					


					if headers[i] in ['MAXSPIKE','AVGSPIKE','TOTSPIKE']:
						if x < 3: x = str(x)  	
						elif x < 10: x = '5' 
						elif x < 30: x = '30' 
						elif x < 60: x = '50' 
						else:        x = '100' 
					elif headers[i] in ['RELN','CXCR4','CHGA','GRIK3','SLC17A7']:
						x = str(x) 
					elif headers[i] in ['R_OBS','OBS']:
						x = str(round(x,-1))
					elif headers[i] in ['R_TOT','TOT']:
						x = str(round(x,-1))
					elif x < 100: 
						x = str(int(line[i]))
					else:
						x = str(round(x,-4))
					key[headers[i]][s] = x
				except ValueError:
					key[headers[i]][s] = line[i]  

	if options.cnts != None:
		cnt_key = dd(lambda: dd(float)) 
		for line in open(options.cnts):
			line = line.split()
			if line[0] == '---': headers = line 
			else:
				s = line[0] 
				for i in range(1,len(headers)):
					cnt_key[headers[i]][s] = float(line[i]) 
				


			

	CR_GENES = ['CXCR4','PCP4','RELN','OLFM1','CLU','CALB2'] 
	CR_EXACT = ['CXCR4','PCP4','RELN','CALB2','OLFM1']
	CR_BIG = ['DBI', 'PCP4', 'CALB2', 'CBLN1', 'OLFM1', 'RELN', 'CXCR4', 'CLU', 'NR2F2']
	#CR_A = ['CXCR4','PCP4','RELN','CLU','CALB2'] 
	#CR_B = ['DBI', 'PCP4', 'CALB2', 'CBLN1', 'OLFM1']
	#CR_C = ['RELN', 'CXCR4', 'CLU', 'NR2F2']
	#CR_D = ['RELN', 'CXCR4', 'PCP4']

	MF_GENES = ['SLC17A7', 'NGEF', 'RAP1GAP', 'VGF', 'GRIK3', 'CYP26A1']
	MF_HI = ['SLC17A7', 'NGEF',  'GRIK3']
	INDIV= True


	MF_0 = ['SLC17A7','SLC1A2','ATP1B1','CHL1']
	MF_1 = [ 'FOXG1', 'CELF4']
	MF_2 = ['GRIA2', 'FOXG1', 'CELF4', 'VCAN', 'STXBP1', 'HS3ST4', 'CHL1']
	MF_3 = ['SSBP3', 'SYT1', 'PDE1A', 'ZBTB18', 'ADRA2A', 'SOBP', 'CCDN2']
	MF_4 = ['SERPINI1', 'RAP2A', 'RTN1', 'SCD5', 'CRYM', 'NELL2', 'ATP1B1']
	MF_5 = ['PBX1', 'NEUROD2', 'SEMA3C', 'NEUROD6', 'ADCY1', 'EIF4G2', 'GRIK2', 'LMO3', 'ENC1', 'PRKACB', 'CAMK2N1', 'ARPP21', 'SLC1A2']

	MINI_1 = ["SEMA3D","GAD1","GAD2"] 
	MINI_1 = ["GAD1","GAD2"] 
	MINI_2 = ["CLMN","ITGA1","XIRP2","DNAH5","WWC1"]


	sp = ScatterPlot(options,xLen=3,yLen=6) 

#	sp.plot_custom_pts(pt_key,key) 
	sp.plot_custom_pts(pt_key,key,AXES=[0,0,2,2],VERSION='FIRING',TITLE='Firing Activity') 
#	sp.plot_custom_pts(pt_key,key) 
	#sp.plot_custom_pts(pt_key,key,cnt_key,genes=['BCL11B','FEZF2','GRIK2','GRIK3','SLC17A7','ATP1B1','NGEF'],neg_genes = ['CXCR4','PCP4','RELN','CALB2','OLFM1','CLU','NR2F2'])
#	sp.plot_custom_pts(pt_key,key,cnt_key,genes=['GRIK2','GRIK3','SLC17A7','ATP1B1','NGEF'],neg_genes = ['CXCR4','PCP4','RELN','CALB2','OLFM1','CLU','NR2F2'])
#	sp.plot_custom_pts(pt_key,key,cnt_key,genes=['GRIK2','GRIK3','SLC17A7','SLC1A2','ATP1B1'],neg_genes = ['GAD1','GAD2','GLS','ABAT'])
#	sp.plot_custom_pts(pt_key,key,cnt_key,genes=['GRIK2','GRIK3','SLC17A7','SLC1A2','ATP1B1'],neg_genes = ['GAD1','GAD2','GLS','ABAT'])
#	sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[1,2,1,1],genes = ['SLC12A2'],neg_genes=['SLC12A5'],BAR_MSG='CHLORIDE SLC12A2 vs SCL12A5')
	#sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[0,3,1,1],BAR_MSG='GLUTE',genes=['SLC17A7','GRIK3','GRIK2','SLC1A2','ATP1B1'],neg_genes = ['GAD1','GAD2','GLS','ABAT'])
	#sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[1,3,1,1],BAR_MSG='VEN',genes=['BCL11B','FEZF2'])

#	sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[0,2,1,1],genes = ['RELN','CXCR4','PCP4','CALB2','OLFM1','CLU','NR2F2'],BAR_MSG='RELN')
#	sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[1,2,1,1],BAR_MSG='SARAF',genes=['SARAF']) #,neg_genes=['RBFOX3'])
	sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[0,3,1,1],BAR_MSG='CAMK2N1',genes=['CAMK2N1'])
	sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[1,3,1,1],BAR_MSG='KCNMB4',genes=['KCNMB4']) #,neg_genes=['RBFOX3'])
	#sp.plot_custom_pts(pt_key,key,cnt_key,AXES=[1,3,1,1],BAR_MSG='MAP2',genes=['MAP2'],neg_genes=['RBFOX3'])


	sp.plot_custom_pts(pt_key,key,AXES=[0,4,2,2],VERSION='LOC',TITLE='Identified CellType') 
	#sp.plot_custom_pts(pt_key,key,cnt_key,genes=['GRIK2','GRIK3','SLC17A7','SLC1A2','ATP1B1','NGEF'],neg_genes = ['GAD1','GAD2','GLS','ABAT'])
	#sp.plot_custom_pts(pt_key,key,cnt_key,genes=['SAT2B','LMO4'])
	#sp.plot_custom_pts(pt_key,key,cnt_key,genes=MINI_2)
#	sp.plot_custom_pts(pt_key,key,cnt_key,genes=MINI_1)
#	sp.plot_custom_pts(pt_key,key,cnt_key,genes=MINI_2)
	sp.finish(prefix) 
	sys.exit() 



	sys.exit() 

			




	 
				















if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-c", "--cnts", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args[0],options)	













