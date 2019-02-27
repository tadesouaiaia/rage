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

COMMON = ['red','blue','green','yellow','orange','purple','lime','cyan','k']
COLORS = COMMON+list(set([x for x in [c for c in COLORS_1+COLORS_2+COLORS_3] if x not in COMMON]))

def rage_error(estr):
	sys.stderr.write(estr+'\n')
	sys.exit() 



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





class dot_log:
        def __init__(self,args,prefix=None,command_line=None,WELCOME=None):

		self.active = True
                self.prefix = prefix
                self.out = sys.stderr
		self.status = None
		self.intro = 'Objective: '
		self.events = []
                if WELCOME != None: self.out.write(WELCOME) 




        def start_major(self,topic,msg=None,subtopics=[]): 


                if self.status in ['minor','major']: self.end()

                self.status = 'major'


                self.blank = ''.join([' ' for x in range(len(self.intro))])

                if msg: 
                	self.out.write(self.intro+topic+'\n') 
                        self.out.write(self.blank+msg+' '+', '.join([str(s) for s in subtopics])+'\n') 
                        self.sub_blank = ''.join([' ' for x in range(len(self.intro+msg))])
                else:
                	self.out.write(self.intro+topic+'...') 
                        self.sub_blank = self.blank 

		self.counter = 0 

	def mark(self,dotrate=0):
                self.counter +=1 
                if dotrate >= 1:   self.out_write('.') 
                elif dotrate > 0 and random() < dotrate: self.out_write('.') 



        def end(self): 

                self.out.write('Complete\n') 
                        
                #self.events[self.topic][self.subtopic] = [self.counter]
                self.subtopic = None 
		self.status = None



	def comment(self,comment): 
                if self.status in ['minor','major']: self.end()
		self.status = None 
		self.out.write('Result: '+comment+'\n') 



	def warn(self,WARNING):

		self.out.write('WARNING: '+WARNING+'\n') 










class GenePlot:
        def __init__(self,options,color_key): #,xLen=3,yLen=4): #,xLen=3,yLen=4,key={}):

                self.options = options
                self.VERBOSE = True
		self.xOffset = 0 
        
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(19.5, 10.5)
                #matplotlib.rcParams['ytick.labelsize'] = 7.5
                
#COLORS_2 = [ 'mediumorchid', 'skyblue', 'gray', 'darkturquoise', 'goldenrod', 'darkgreen', 'floralwhite', 'darkviolet', 'darkgray', 'moccasin', 'saddlebrown', 'grey', 'darkslateblue', 'lightskyblue', 'lightpink', 'mediumvioletred', 'slategrey', 'red', 'deeppink', 'limegreen', 'darkmagenta', 'palegoldenrod', 'plum', 'turquoise', 'lightgrey', 'lightgoldenrodyellow', 'darkgoldenrod', 'lavender', 'maroon', 'yellowgreen', 'sandybrown', 'thistle', 'violet', 'navy', 'magenta', 'dimgrey', 'tan', 'rosybrown', 'olivedrab', 'blue', 'lightblue', 'ghostwhite', 'honeydew', 'cornflowerblue', 'linen', 'darkblue', 'powderblue', 'seagreen', 'darkkhaki']
                
		FACECOLOR = args.facecolor 

		matplotlib.rcParams['savefig.facecolor'] = FACECOLOR
                self.fig.patch.set_facecolor(FACECOLOR)
                self.fig.set_facecolor(FACECOLOR) 
                seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
                sns.set(rc={'axes.facecolor':FACECOLOR, 'figure.facecolor':FACECOLOR})

                self.xLoc, self.yLoc = 0,0 
		self.iKey = dd(int) 

		self.color_key = color_key 

		self.iKey['MEGA'],self.iKey['ES'],self.iKey['CR'],self.iKey['ES'],self.iKey['ADULT'] = 1,2,3,4,5
		
		self.color_offset = 0
                self.yLen =  1
		self.xLen =  len(args.genes) 
		self.axes = [] 
		self.legend_items = {} 



	def finish(self,out='BASIC',title=None,xLen=2,yLen=3,kind='RAW'):

		#if self.options.title != None:
		#	plt.suptitle(self.options.title,fontsize=40,fontweight='bold') 

  		plt.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.90,wspace=0.10,hspace=0.05)
		
		if self.options.output !=None: plt.savefig(self.options.output+'.png',dpi=200) 
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



		



	def add_legend(self):
		nc,fs,cs,hp,bb1,bb2 = len(self.legend_items),20 ,1,1,0.5,-0.1
		leg_data = sorted([[self.iKey[x[0]],x[0],x[1]] for x in self.legend_items.items()])
		labels,items = [ld[1] for ld in leg_data],[ld[2] for ld in leg_data]			
		leg = self.axes[-1].legend(items,labels,handletextpad=hp,columnspacing=cs,  fontsize=fs,ncol=nc,bbox_to_anchor=(bb1,bb2),loc='upper center')



	def plot_boxes(self,box_data,SCALE=5): 

		
		self.scale_key, self.xOffset, self.yTop = dd(list) , 0 , 0 
		
		self.axes.append(plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1))

		box_names = [bd[0] for bd in box_data]
		self.scale_key = dd(list) 
		for gene_name,gene_vals in box_data:
			self.add_multi_box(gene_name,gene_vals,SCALE=SCALE)

                scale_sums  = {opt: [sum([vals[x][y] for x in range(len(vals))]) for y in range(len(vals[0]))] for opt,vals in self.scale_key.items()}

		if not self.options.noscale:
                	self.add_multi_box('SCALE TOTAL',{opt: [sum([vals[x][y] for x in range(len(vals))]) for y in range(len(vals[0]))] for opt,vals in self.scale_key.items()})
		


		self.axes[-1].set_xlim(0-self.width,self.xOffset)  
		self.axes[-1].axis('off') 
		self.xLoc +=1 













	def add_multi_box(self,TITLE,h_data,BW=False,SCALE=0):


		data = [d for d in sorted([[self.iKey[d[0]],d[0],d[1]] for d in h_data.items()]) if d[1] not in ['MA','UN']]

		opts,vals,items,labels,R,pv = [d[1] for d in data],[[v for v in  d[2]] for d in data],[],[],'NA','NA'

#		colors = [self.color_key[X] for X in opts] 

#		colors,global_sigs,local_sigs,jitters = [self.get_color(X) for X in opts], [1 for X in opts],[1 for X in opts],[] 
		colors,global_sigs,local_sigs,jitters = [self.color_key[X] for X in opts], [1 for X in opts],[1 for X in opts],[] 

		if type(colors[0]) == tuple:	colors,alt_colors = [c[0] for c in colors], [c[1] for c in colors]
		else:				alt_colors = [c for c in colors]
		
                self.pos, self.width = [self.xOffset+0.20*x for x in range(len(opts))], 0.16
		all_vals = [a for b in vals for a in b] 

		valMin,valMean,valMax, valMeans = min(all_vals), np.mean(all_vals), max(all_vals)*1.25, [np.mean(v) for v in vals]
		for i in range(len(vals)):

			other_vals = [a for b in [vals[k] for k in range(len(vals)) if k != i] for a in b] 
			pv = stats.ttest_ind(vals[i],other_vals)[-1] 
			if valMeans[i] > valMean and pv < global_sigs[i]: global_sigs[i] = pv
	
                	jitters.append([[self.pos[i] + np.random.normal(0,0.02) for m in range(len(vals[i]))],[v*np.random.normal(1,0.005) for v in vals[i]]])
			for j in range(i+1,len(vals)): 
				pv = stats.ttest_ind(vals[i],vals[j])[-1] 
				if valMeans[i] > valMeans[j] and pv < local_sigs[i]: local_sigs[i] = pv 
				elif valMeans[i] < valMeans[j] and pv < local_sigs[j]: local_sigs[j] = pv 

			if SCALE > 0:
				scaler = MinMaxScaler(feature_range=(0,SCALE))
				self.scale_key[opts[i]].append([vs for vs in scaler.fit_transform(np.array(vals[i],dtype=float).reshape(-1,1)).reshape(1,-1)[0]])
		self.fill_boxes(opts,vals,colors,alt_colors,JITTERS=jitters,SIG_STARS=(global_sigs,local_sigs))
		if valMax > self.yTop: self.yTop = valMax
		self.xOffset = self.pos[-1]+2*self.width             
		
    
		if TITLE:	
			self.axes[-1].text(np.mean(self.pos),0-(self.yTop)*0.1,TITLE.split(';')[-1],fontsize=16,fontweight='bold',horizontalalignment='center',verticalalignment='top') 
		
		return 



	def fill_boxes(self,opts,vals,colors,alt_colors,JITTERS=None,SIG_STARS=None):

		G_THRES,L_THRES = 0.05, 0.05
                self.bp = self.axes[-1].boxplot(vals,positions=self.pos, widths=self.width, patch_artist=True, showmeans=True,whis=0.7)
                for i,opt in enumerate(opts):

                        clr,ac = colors[i], alt_colors[i] 
                        self.bp['boxes'][i].set_edgecolor(clr) 
                        self.bp['boxes'][i].set_linewidth(1) 
                        plt.setp(self.bp['medians'][i], color=clr, linewidth=3)
                        plt.setp(self.bp['means'][i], marker='h',markersize=9,markerfacecolor=clr) 
                        plt.setp(self.bp['caps'][(i*2)+1], color=clr,linewidth=1) 
                        plt.setp(self.bp['caps'][(i*2)], color=clr,linewidth=1) 
                        plt.setp(self.bp['whiskers'][i*2], color=clr,linewidth=1) 
                        plt.setp(self.bp['whiskers'][1+(i*2)], color=clr,linewidth=1)   
                        plt.setp(self.bp['fliers'][i], markerfacecolor=clr, markeredgecolor = clr, marker='s',markersize=2.0)           



			if clr == ac:	self.legend_items[opt] = Rect((0,0),1,1,fc=clr)
			else:   	self.legend_items[opt] = Rect((0,0),1,1,fc=a,ec=b,hatch='*')
			
                for patch, clevel in zip(self.bp['boxes'], colors):
                        patch.set_facecolor(clevel); patch.set_alpha(0.5) 
		if JITTERS != None:
			for i,(xJ,yJ) in enumerate(JITTERS):
                        	plt.scatter(xJ, yJ, c=colors[i], alpha=0.7,s=4,zorder=9)
		if SIG_STARS != None:
			for i,(g,l) in enumerate(zip(SIG_STARS[0],SIG_STARS[1])): 
				if l < L_THRES: plt.scatter(self.pos[i],max(vals[i])*1.1,s=pv_2_size(l),marker='*',color='silver',edgecolor='black',linewidth=1)
				if g < G_THRES: plt.scatter(self.pos[i],(max(vals[i])*1.1)+2,s=1.5*pv_2_size(l),marker='*',color='gold',edgecolor='black',linewidth=1)
			



def pv_2_size(s):
	MAX,LESS = 200,20  

	if   s < 0.0000001: return MAX
	elif s < 0.000001:  return MAX-(LESS)
	elif s < 0.00001:   return MAX-(LESS*2)
	elif s < 0.0001:    return MAX-(LESS*3) 
	elif s < 0.001:     return MAX-(LESS*4) 
	elif s < 0.01:      return MAX-(LESS*5)
	return MAX-(LESS*6) 














def run_script(args):
	k=0
	TOPICS = ['CT'] 

	key,gene_key = dd(bool), dd(bool) 
	progress = dot_log(args,WELCOME='\nStarting Rage Plot\n') 
	progress.start_major('Reading Key File') 

	for line in args.key:
		line = line.split() 
		if line[0] == '---': 
 			headers = line  
			if args.id == None: args.id = headers[1] 
			id_idx = headers.index(args.id) 
		else:
			key[line[0]] = line[id_idx] 

	
	color_key = {kv: False for kv in list(set(key.values()))}
	g_str = " ".join([x if i%2==0 else '('+str(x)+'),' for (i,x) in enumerate([a for b in cc(key.values()).items() for a in b])])[0:-1]
	progress.comment('Key ID: '+args.id+': Groups Found: '+g_str)


	gene_groups = [[g.upper() for g in gene_group.split(',')] for gene_group in args.genes] 
	genes       = [a for b in gene_groups for a in b] 
	for gene in genes: gene_key[gene] = True 

	group_key = dd(lambda: dd(list)) 
	mean_key = dd(float) 
	if len(args.colors) > 0: 
		for ac in args.colors: 
			a,b = ac.split('=')  
			if len(a.split(',')) > 1: 
				grp_id,grp_ids = a,[x.upper() for x in a.split(',')] 
				color_key[grp_id] = b 
				for x in a.split(','): 
					for s in key:
						if key[s] in grp_ids: key[s] = grp_id 
			else:
				if a.upper() in color_key:	color_key[a.upper()] = b 
				elif a in color_key: 		color_key[a]         = b 
				else: 	progress.warn('ID '+a+' not found in key file') 
		ignored =  [x for x in color_key.keys() if not color_key[x]]
		if len(ignored)>0:
			progress.warn('IDS '+",".join(ignored)+' do not have supplied colors and will be ignored')
	else: 
		for i,x in enumerate(color_key.keys()): color_key[x] = COLORS[i] 


	progress.start_major('Reading Count File') 
	for line in args.counts:
		line = line.split() 
		if line[0] == '---': samples = line[1::] 	
		elif not gene_key[line[0].split(';')[-1]]: continue  
		else:		
			gene = line[0].split(';')[-1] 
			cnts,all_cnts = [log(float(x)+1,2) for x in line[1::]] ,[]
			for i,c in enumerate(cnts):
				if color_key[key[samples[i]]]: 
					group_key[gene][key[samples[i]]].append(c)
					all_cnts.append(c) 
			mean_key[gene] = np.mean(all_cnts) 



	missing_genes = [gene for gene in genes if gene not in group_key] 
	if len(missing_genes) > 0: rage_error('ERROR: The following genes are not in counts file: '+",".join(missing_genes))

	progress.start_major('Creating Plots') 
	bp = GenePlot(args,color_key)
	for gene_group in gene_groups:
		group_mean = int(1+np.mean([mean_key[gene] for gene in gene_group]))
		bp.plot_boxes([(gene,group_key[gene]) for gene in gene_group],SCALE=group_mean) 

	bp.add_legend() 
	bp.finish() 

	progress.end() 

	sys.exit() 
				







if __name__ == '__main__':

        import sys,os,argparse,shutil
                                
        class MyParser(argparse.ArgumentParser):
                formatter_class=argparse.RawTextHelpFormatter
                def error(self, detail=None,values=None,choice=None):
                        
                        if not detail: self.print_help()
                        elif detail and values and choice:
                                sys.stderr.write(str('InvalidOption:\n'))
                                sys.stderr.write(str(values+" "+" ".join(choice)+'\n'))
                                sys.stderr.write(str(detail)+'\n\n')
                                self.print_help()
			
                        elif detail:
                                sys.stderr.write(str(detail)+'\n\n')
                                self.print_help()	
			self.show_examples() 
			self.show_colors() 
			sys.exit() 

			
		def show_colors(self):
			sys.stderr.write('\n') 
			sys.stderr.write('A list of possible colors:\n')
			for i,c in enumerate(COLORS): 
				if i == 0: sys.stderr.write(c) 
				elif i % 17 == 0: sys.stderr.write('\n'+c) 		
				else: sys.stderr.write(','+c) 		
				

			sys.stderr.write('\n') 
			sys.stderr.write('\n') 




		def show_examples(self):
			sys.stderr.write('Examples: \n\n')
			ex1='rage_plot_bars.py --key KEY.txt --counts CNTS.txt --genes CRMP1,ACTB MAP2,SLC12A5 --id CELLTYPE  --facecolor mintcream --noscale'
			ex2='rage_plot_bars.py --key KEY.txt --counts CNTS.txt --genes CRMP1,ACTB MAP2,SLC12A5 --id CELLTYPE --colors CR=red ADULT=blue --facecolor orchid'
			ex3='rage_plot_bars.py --key KEY.txt --counts CNTS.txt --genes CRMP1,ACTB --id CELLTYPE --colors CSP,ISVZ=red ADULT=blue --facecolor orange'
			sys.stderr.write('EXAMPLE 1: USE ALL MEMBERS OF CATEGORY CELLTYPE (NO SCALE):\n') 
			sys.stderr.write(ex1+'\n\n') 
			sys.stderr.write('EXAMPLE 2: USE ONLY MEMBERS CR/ADULT by assigning colors to them\n') 
			sys.stderr.write(ex2+'\n\n') 
			sys.stderr.write('EXAMPLE 3: COMBINE THE MEMBERS CSP/ISVZ by assigning them to the same color\n')
			sys.stderr.write(ex3+'\n\n') 


        parser=MyParser(formatter_class=argparse.RawTextHelpFormatter)
        parser.set_defaults(dataset=None)
        help_format="%-40s %40s"
        fc = argparse.RawTextHelpFormatter

      	parser.add_argument('--key',type=argparse.FileType('r'),default=None,help='Key File (Reqired)')
      	parser.add_argument('--counts',type=argparse.FileType('r'),default=None,required=True,help='Counts File (Required)')
     	parser.add_argument('--id',type=str,default=None,help='Key Column ID For Groups')
 #    	parser.add_argument('--title',type=str,default=None,help='Plot Title')
     	parser.add_argument('--output',type=str,default=None,help='Output File Prefix to save Figure') 
        parser.add_argument('--colors',type=str,nargs='+',default=[],help='Colors for Groups')
        parser.add_argument('--genes',type=str,nargs='+',required=True,help='Genes for Plot')
#	parser.add_argument('--silent',action='store_true', dest='silent',default=False,help='Silence the progress log') 
	parser.add_argument('--noscale',action='store_true', dest='noscale',default=False,help='Dont show the scale total') 
     	parser.add_argument('--facecolor',type=str,default='snow',help='Color for plot face') 

	args = parser.parse_args()

	run_script(args)	













