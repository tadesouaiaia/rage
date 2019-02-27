#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from collections import Counter as cc
from scipy import stats
from math import log 
from scipy.stats import chisquare
import numpy as np 
fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}




def scan(cnts_file,options):
	ktype = {} 
	kidx  = {} 
	for line in open(cnts_file):
		line = line.split()

		if line[0] == '---':
			cats = line[1::] 
			key  = dd(list)
		else:
			data = line[1::]

	    		for i in range(len(data)): 
				key[cats[i]].append(data[i]) 
	    #R=stats.pearsonr(cntsI,cntsJ)[0]
	    #S=stats.spearmanr(cntsI,cntsJ)[0]
	    #R=stats.pearsonr(iV,jV)[0]
		#iS=stats.spearmanr(iV,jV)[0]
	
	for i in range(len(cats)):
		ktype[cats[i]],cd  = 'continuous', key[cats[i]] 
		try: 			key[cats[i]] = [float(x) if x!= 'NA' else 'NA' for x in key[cats[i]]] 
		except ValueError: 		ktype[cats[i]] = 'binary' 

		kidx[cats[i]] = [j for j in range(len(cd)) if cd[j] != 'NA'] 
		


	for i in range(len(cats)): 
		i_type = ktype[cats[i]] 
		i_all  = key[cats[i]] 
		
		for j in range(i+1,len(cats)): 			
			j_type = ktype[cats[j]] 
			j_all  = key[cats[j]] 

			b_idx = [a for a in kidx[cats[i]] if a in kidx[cats[j]]] 

			if len(b_idx) < 25: continue

			im = [i_all[k] for k in b_idx]
			jm = [j_all[k] for k in b_idx]


			if ktype[cats[j]] == 'continuous' and ktype[cats[i]] == 'continuous': 

	    			R,pv =stats.pearsonr(im,jm)
				print cats[i],cats[j],len(b_idx), 'continuous','continuous',R,pv 

	
			elif ktype[cats[j]] == 'binary' and ktype[cats[i]] == 'binary': 

				jcc =  cc(jm).items() 
				j_total = sum([jc[1] for jc in jcc]) 
				jpercs = {jk: jv/float(j_total) for jk,jv in jcc}
				dubs = [(im[z],jm[z]) for z in range(len(im))]  
				for x in list(set(im)): 
					jx = [d[1] for d in dubs if d[0] == x]  
					if len(jx) < 10: continue 
					j_obs = cc(jx) 
					j_tot = sum(j_obs.values()) 
					my_obs = [] 
					my_exp = [] 
					my_ids = [] 
					for jk,jv in jpercs.items(): 
						if jk in j_obs: my_obs.append(j_obs[jk]) 
						else:           my_obs.append(0) 
						my_exp.append(jv*j_tot) 
						my_ids.append(jk) 

					maxE = sorted([(jO/float(jE),jN,jO,jE) for jO,jE,jN in zip(my_obs,my_exp,my_ids)],reverse=True)[0] 

					print cats[i], cats[j], 'bin/bin assoc:', x, sum(my_obs), '|', maxE[1],maxE[0],maxE[2],maxE[3],'|', chisquare(my_obs, f_exp = my_exp)[1]  


			else:
				if ktype[cats[i]] == 'binary': 
					bT,cT,bN,cN = im, jm, cats[i], cats[j] 
				else: 
					bT,cT,bN,cN = jm, im, cats[j], cats[i] 


				 
				dubs = [(bT[z],cT[z]) for z in range(len(im))]  
				bList = list(set(bT)) 				
				b_groups = {bK: [d[1] for d in dubs if d[0] == bK] for bK in bList}
				b_means = sorted([(np.mean(b_groups[bK]),bK) for bK in b_groups if len(b_groups[bK])>10]) 
				if len(b_means) < 2: continue  

				b_min,b_max = b_means[0][1], b_means[-1][1]
				pv = stats.ttest_ind(b_groups[b_min],b_groups[b_max])[1]  
				print cats[i], cats[j], 'ttest', b_min, b_max, len(b_groups[b_min]), len(b_groups[b_max]),'means', np.mean(b_groups[b_min]), np.mean(b_groups[b_max]), pv 



				






if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()



	scan(args[0],options)













