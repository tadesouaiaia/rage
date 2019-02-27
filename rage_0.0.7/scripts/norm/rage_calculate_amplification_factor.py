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
from scipy.stats import variation as CV


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










class GeneAmps:
        def __init__(self,samples,gene_key): #,xLen=3,yLen=4): #,xLen=3,yLen=4,key={}):

       
		self.samples = samples
		self.genes, self.cnts    = gene_key.keys(), gene_key 
 
		self.amp_key = {} 


	def calculate_stats(self,amp_cnts): 

		self.amp_cnts = amp_cnts
		self.alt_key = {a: b for a,b in amp_cnts.items()}

		self.calculate_sample_stats() 

		inflation = sorted([(self.verdicts[s]['inflation'],s) for s in self.samples]) 
		self.calculate_gene_stats() 
	
		for i,(scr,s) in enumerate(inflation): 
			self.calculate_altered_inflation(s) 

		
		sys.exit() 

		self.calculate_gene_stats() 

	def calculate_altered_inflation(self,s): 


		for gene in self.genes: 

			t,a,n = self.cnts[gene][s],self.amp_cnts[gene][s],self.alt_key[gene][s] 
			print t,a,n				
				


	def calculate_gene_stats(self):

		self.gene_verdicts = {} 
		drops = [self.verdicts[s]['dropouts'] for s in self.samples] 
		for gene in self.genes: 
			cnts =  [self.cnts[gene][s] for s in self.samples] 
			amps =  [self.amp_cnts[gene][s] for s in self.samples] 
			R,pv = stats.pearsonr(cnts,drops) 
			verdicts = dd(float) 
			verdicts['total'] = sum(cnts) 
			for s in self.samples: 
				c,a = self.cnts[gene][s], self.amp_cnts[gene][s] 
				if a > c: 
					if c == 0: verdicts['dropouts'] +=1 
					verdicts['deflate'] += 1
				elif c > a: 	verdicts['inflate'] += 1
				else:		verdicts['match'] += 1 	

				verdicts['diff'] =  math.fabs(c-a) 
			#print gene,verdicts['diff'],verdicts['total'],verdicts['match'], CV(cnts)/CV(amps) 
			self.gene_verdicts[gene] = verdicts	
 





	def calculate_sample_stats(self):
		print '---','total','etotal','obs','eobs','inflation-factor','eflation-factor'
		self.verdicts = {} 	
		for s in self.samples: 
			dropouts = 0
			verdicts = dd(int)  
			for gene in self.genes: 
				mod_cnts = dd(float)
				try:  
					tc, ec  = self.cnts[gene][s], self.amp_cnts[gene][s]
				except KeyError:
					continue 

				if tc > 0: verdicts['logs'] += math.log(tc,2)  
				if ec > 0: verdicts['elogs'] += math.log(ec,2)  
				if tc > 0: verdicts['obs']  +=1
				if ec > 0: verdicts['eobs'] +=1
				if tc < ec:
 					verdicts['lost'] += ec - tc 
					if tc == 0: 	verdicts['dropouts']  += 1
					else:		verdicts['deflation'] += 1 						
				elif tc > ec: 	
					verdicts['inflation'] += 1	
 					verdicts['found'] += tc - ec 
				
				verdicts['total'] += tc 		
				verdicts['etotal'] += ec	
				verdicts['net'] += tc - ec 

			verdicts['inflation'] =  verdicts['logs'] / math.log(verdicts['total'],2)
			verdicts['eflation'] =  verdicts['elogs'] / math.log(verdicts['etotal'],2)
			self.verdicts[s] = verdicts 


			print s, verdicts['total'], verdicts['etotal'], verdicts['obs'],verdicts['eobs'], verdicts['inflation'],verdicts['eflation'] 

	
			#print s, math.log(verdicts['total'],2), verdicts['logs'], verdicts['logs'] / math.log(verdicts['total'],2)
		#sys.exit() 
		sys.exit() 		


def run_script(args):
	k=0

	key,gene_key = dd(bool), dd(bool) 
	#progress = dot_log(args,WELCOME='\nStarting Rage Plot\n') 
	#progress.start_major('Reading Key File') 

	raw_samples = args.raw.readline().split()[1::] 
	raw_cnts = {} 
	for line in args.raw:
		line = line.split() 
		cnts = [float(x) for x in line[1::]] 
		raw_cnts[line[0]]  = {s: cnts[i] for i,s in enumerate(raw_samples)} 


	for (a,b) in [('inferred',args.inferred),('imputed',args.imputed)]:
		if b != None: 
			geneAmps = GeneAmps(raw_samples,raw_cnts) 
			i_cnts = {} 
			i_samples = b.readline().split()[1::] 
			for line in b:
				line = line.split() 
				cnts = [float(x) for x in line[1::]] 
				i_cnts[line[0]]  = {s: cnts[i] for i,s in enumerate(i_samples)} 
			
			geneAmps.calculate_stats(i_cnts) 
		



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

      	parser.add_argument('--raw',type=argparse.FileType('r'),default=None,help='Raw Counts File (Reqired)')
      	parser.add_argument('--imputed',type=argparse.FileType('r'),default=None,required=False,help='Imputed Counts File (Required)')
      	parser.add_argument('--inferred',type=argparse.FileType('r'),default=None,required=False,help='Inferred Counts File (Required)')
 #    	parser.add_argument('--title',type=str,default=None,help='Plot Title')

	args = parser.parse_args()

	run_script(args)	













