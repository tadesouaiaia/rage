#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as mlt 
import statsmodels.sandbox.stats.multicomp as mpt 
#import statsmodels.sandbox.stats.multicomp.fdrcorrection0 as fdr

#statsmodels.sandbox.stats.multicomp.multipletests(p,alpha=0.05,method='fdr_bh')
#statsmodels.sandbox.stats.multicomp.fdrcorrection0(p,alpha=0.05)


#def parse_out_file(line): 

#			--- RS CV obs len parent maxG maxMeans maxObs maxChi maxP | params

#			['ENSG00000000971;chr1;CFH', '0.007', '3.537', '276', '2455', 'FULLGROUP', 'FULLGROUP~ES', '0.34', '0.19', '2.6e-05', '1.3e-03', '|', 'FULLGROUP~AB', 'False', '2.51e-03', '2.51e-03', '|', 'FULLGROUP~EB', 'False', '3.00e-05', '3.00e-05']




IGNORE_LIST = ['TOT','MINTYPE','MT','PTN','PTU','BGW','|','SIN',"MULT",'BSPIKES','BSPIKE','ZT','MZTYPE','MAX_LATENCY']

XK = {'PLOC': 'LOC', 'LOC': 'LOC','vBIO': 'BIO','BIO': 'BIO','PBIO': 'BIO', 'ULOC': 'LOC','MLOC': 'LOC','FS': 'FS','XLOC': 'LOC','P_100': 'PASS', 'P_125': 'PASS'} 

FK = {'GW': 'GW', 'vGW': 'GW', 'BGW': 'GW', 'OBS': 'OBS', 'TOT': 'TOT', 'RELSIZE': 'RELSIZE', 'AREA': 'AREA', 'RSIZE': 'RSIZE','R_OBS':'OBS','R_TOT': 'TOT'} 

XK['LONGSTRING'] = 'LONGSTRING'
XK['MATURE'] = 'MATURE'
XK['FIRING_STYLE'] = 'FS'
XK['SURE_LOC'] = 'LOC'
XK['bio'] = 'BIO'
XK['loc'] = 'LOC'
XK['sureloc'] = 'LOC'
XK['gw'] = 'GW'


for x in ['MAXSPIKES', 'MPIKE','MAX_SPIKES']: 
	FK[x] = 'MAXSPIKES'

for x in ['mega','MINTYPE','FK','MO','PTU','PTN','MEGA','PTC','GICU','GIN','GICU','GCN','VT','VA','ST','UTYPE','GICN','P1','P2','P3','P4','P5','CT','MT','XTYPE','TYPE']: 
	XK[x] = 'MEGA'

for x in ['TL', 'UTL','XTL']: 
	XK[x] = 'SIZE-TYPE'

FK['relative_size'] = 'RSIZE'
FK['meansize'] = 'MSIZE'
XK['orig_size'] = 'OSIZE'
FK['sizes']     = 'SIZES' 


def run_script(data_files,options):

	grp_key = dd(list) 
	h_key = dd(lambda: {})  	
	samples = [] 


	merge_key = dd(lambda: dd(list)) 

	Xs,Fs = set([]), set([]) 

	size_key = dd(lambda: dd(int)) 
	bio_key = dd(lambda: dd(int)) 
	loc_key = dd(lambda: dd(int)) 
	gw_key = dd(lambda: dd(int))
	rel_key = dd(lambda: dd(int)) 
	area_key = dd(lambda: dd(int))  
	obs_key = dd(lambda: dd(int))  
	fire_key = dd(lambda: dd(int))  

	swap_dict = {'UNCLEAR': 'UN', 'VZ': 'SVZ','SVZ_VZ': 'SVZ', 'CP': 'CP','MI': 'MINI','MG': 'MEGA','IZ_SVZ': 'SVZ_IZ','SP_IZ':'IZ_SP','EBSC': 'ES'} 

	K = dd(lambda: dd(lambda: dd(int)))
	IGNORE_LIST = ['|','SIN','MULT'] 
	V_IGNORE = ['MIA','MB'] 
	AX_LIST = ['MINTYPE','MT','PTN','PTU','ZT','MZTYPE']

	for f in data_files:
		for line in open(f): 
			line = line.split() 
			if line[0] == '---': headers = line 
			else:
				sample = line[0] 
				for i in range(1,len(line)):
					h,v = headers[i],line[i] 
					if v in swap_dict: 	v = swap_dict[v]
					if h in IGNORE_LIST: continue 
					if v in V_IGNORE: continue 
					if h in AX_LIST and v != 'NA':
						if v in ['MINI','MEGA']: 
							K[sample]['ax-size'][v]+=1 
						elif v in ['CR']: 
							K[sample]['ax-type'][v]+=1 
						elif v in ['CRMZ','SG']: 
							K[sample]['cr-type'][v]+=1 
						elif v in ['ES','EBSC']: 
							K[sample]['LOC'][v] +=1 
						elif v in ['UNK']: continue 
						else:
							print 'yo',h,v
						continue 
					#if v in ['NA']: continue


					
					try:			v = float(v) 
					except ValueError:	v = v	 

					if h in ['mega']: K[sample]['OLD-SIZE'][v]+=1 
					elif h in ['orig_size']: K[sample]['DETAIL-SIZE'][v]+=1 
					elif v in ['MEGA','MINI','UN']: 	K[sample]['SIZE'][v]+=1 
					elif v in ['CR','CRMZ']: 	K[sample]['CR']['TRUE']+=1 
					elif h in ['BIO','vBIO']: 	K[sample]['BIOPSY'][v]+=1 
					elif h in ['PBIO']: 	K[sample]['BIO-GROUP'][v]+=1 
					elif v in ['MEGAF']: 	K[sample]['MF'][v]+=1 
					elif v in ['CP_SP','SVZ_IZ','CP','SP','IZ','SVZ','MZ','IZ_SP','ES','SP_IZ_SVZ']: 		K[sample]['LOC'][v]+=1 
					elif h in ['vGW','GW']:			K[sample]['GW'][v]+=1 
					elif h in ['BGW']:			K[sample]['GW-GROUP'][v]+=1 
					elif h in ['BSPIKE','BSPIKES']:			K[sample]['SPIKE-GROUP'][v]+=1 
					elif h in ['RSIZE','RELSIZE']:	K[sample]['RELSIZE-A'][v]+=1 
					elif h in ['relative_size']:	K[sample]['RELSIZE-B'][v]+=1 
					elif h in ['AREA']:		K[sample]['AREA'][v]+=1
					elif h in ['OBS','R_OBS']:		K[sample]['GENE_OBS'][v]+=1 
					elif h in ['R_TOT','TOT']:		K[sample]['TOTAL_READS'][v]+=1 
					elif h in ['FS']:		K[sample]['FIRING_TYPE'][v]+=1
					elif h in ['MAX_LATENCY']:		K[sample]['MAX_LATENCY'][v]+=1
					elif h in ['P_100','P_125']:		K[sample][h][v]+=1
					elif h in ['MAXSPIKES','MPIKE','MAX_SPIKES']:		K[sample]['MAXSPIKES'][v]+=1
					elif h in ['FIRING_STYLE']:		K[sample]['FIRING_STYLE'][v]+=1
					elif h in ['MATURE']:		K[sample]['MATURE'][v]+=1
					elif v in ['MATURE']:		K[sample]['MATURE']['YES']+=1
					elif v in ['SG']:		K[sample]['SUB-GRAN']['YES']+=1
					else:
						if v in ['UNK','UNKNOWN']: 
							if h in ['LOC','SURE_LOC','ULOC']: K[sample]['LOC']['UNK']+=1 
						
							
						else:
							if len(str(v).split('-')) == 2: # and v != '-':		
								a,b = v.split('-') 
								K[sample]['GENERAL_LOC'][b] += 1
								K[sample]['GENERAL_SIZE'][a] += 1
								K[sample]['SIZE_LOC'][v] += 1

							elif v in ['CSP','ISVZ']: K[sample]['GEN-LOC'][v] +=1 

							elif v in ['NA','UN','EB']: continue 
							else:

								if h == 'LONGSTRING': 
									for x in v.split('@'): 
										if x in ['NA','UNKNOWN','UNK']: continue
										if x in swap_dict: x = swap_dict[x] 
										if x in ['CP_SP','SVZ_IZ','CP','SP','IZ','SVZ','MZ','IZ_SP','ES','SP_IZ_SVZ']: 	K[sample]['LOC'][x]+=1 
										elif x in ['MEGA','MINI']: K[sample]['SIZE'][x] += 1 
										elif x in ['CSP','ISVZ']: K[sample]['GEN-LOC'][x] +=1 
										elif x in ['CR','CRMZ']: 	K[sample]['CR']['TRUE']+=1 
										elif x in ['MF']:		K[sample]['FIRING_TYPE'][x]+=1
										else:
											print h,v

								

					
					continue 
					if h in IGNORE_LIST: continue 
					elif h in XK:
						merge_key[sample][XK[h]].append(v)
						Xs.add(XK[h]) 

					elif h in FK:
						try: merge_key[sample][FK[h]].append(float(v)) 
						except ValueError:	merge_key[sample][FK[h]].append(v)
					
						Fs.add(FK[h]) 
					elif v != 'NA' and v != 'UNK':
						print h,v				


	SIZING = ['RELSIZE-A','RELSIZE-B'] 
	BIOPSY = ['BIOPSY','bio'] 
	BIOPSY_GROUP = ['BIO-GROUP']  
	GW = ['GW'] #'GW-GROUP'] 
#	GW = ['GW','GW-GROUP'] 
	MEGA = ['SIZE','GENERAL_SIZE','DETAIL-SIZE','OLD-SIZE','ax-size'] 
	SIZE_LOC = ['SIZE_LOC'] 
	LOC = ['GEN-LOC','LOC','GENERAL_LOC'] 


	NOTES = ['SUB-GRAN','MATURE','CR','MF','ax-type','cr-type'] 
#	RNA_INFO = ['TOTAL_READS','GENE_OBS','P_100','P_125','AREA'] 
#	EFIZ = ['MAXSPIKES','SPIKE-GROUP','MAX_LATENCY']
#	EFIZ_BIN = ['FIRING_STYLE','FIRING_TYPE'] 
	EFIZ = ['MAXSPIKES','SPIKE-GROUP','MAX_LATENCY','FIRING_STYLE','FIRING_TYPE'] 

	my_list = [] 	
	SK = ['MEGA','LOC', 'GLOC', 'PAIRED-LOC'] 

	NK = ['SUB-G', 'MATURE', 'MFIRE', 'CR']
#	RK = ['GW', 'BIO','BIO-GROUP', 'QC', 'AREA','REL-SIZE', 'RSIZE', 'total', 'obs']
	RK = ['GW', 'BIOP','B-GRP','AREA','REL-A','TOT-READS', 'TOT-GENES','QC']
	RZ = ['GW', 'BIOP','B-GRP','AREA','REL-A','TOT-READS', 'TOT-GENES','QC_1,QC_2']
	print "\t".join(['---']+SK+NK+RZ+['FTYPE','FNUM','FGRP']) 



	for s in K.keys(): 
		tmp,res,eres = {},{},{} 
		notes,ms,km = dd(bool),dd(list),dd(lambda: dd(int)) 
		for tl,t in [(MEGA,'SIZE'),(BIOPSY_GROUP,'BIO-GROUP'),(GW,'GW'),(BIOPSY,'BIO'),(SIZING,'REL-SIZING'),(LOC,'LOC'),(SIZE_LOC,'SIZE_LOC')]:
			for k in tl:
				for n,v in K[s][k].items(): km[t][n]+=v 
			ms[t]  = sorted(km[t].items(),key=lambda X: X[1],reverse=True) 			

		res['REL-SIZE'] = sorted([x[0] for x in ms['REL-SIZING']])[0]
		
		if res['REL-SIZE'] == 'NA': 
			res['RSIZE'] = 1.0 
		else:
			res['RSIZE'] = float(res['REL-SIZE'])



		for tv in ['GW','BIO','BIO-GROUP']: 
			if  len(ms[tv]) == 0: res[tv] = 'NA'
			else:		      res[tv] =  ms[tv][0][0]  




		res['AREA'] = [x for x in K[s]['AREA'].keys()+['NA']][0] 
		tmp['obs'] = sorted([x for x in K[s]['GENE_OBS'].keys() if x != 'NA'])
		tmp['total']  = sorted([log(y) if y > 100 else y for y in [x for x in K[s]['TOTAL_READS'].keys() if x != 'NA']])
		
		for tv in tmp:
			if len(tmp[tv]) == 0: res[tv] = ['NA','NA']
			else: 		      res[tv] = [tmp[tv][0],tmp[tv][-1]]

		if 'FAIL' in K[s]['P_100'].keys() and 'FAIL' in K[s]['P_125']:  res['QC'] =  'FAIL,FAIL' 
		elif 'FAIL' in K[s]['P_100'].keys() or 'FAIL' in K[s]['P_125']: res['QC'] = 'PASS,FAIL' 
		else:								res['QC'] = 'PASS,PASS' 	

		for tv in EFIZ:	tmp[tv] =  [x for x in K[s][tv].keys() if x != 'NA'] 
			
		for tv in ['MAXSPIKES','MAX_LATENCY']: 
			if len(tmp[tv]) == 0: eres[tv] = 'NA' 
			else: 		      eres[tv] = int(np.mean(tmp[tv]) )



		tmp['FIRING_TYPE'] =  [tx for tx in tmp['FIRING_TYPE'] if tx != 'MF'] 
		for ex in ['SPIKE-GROUP','FIRING_TYPE']: 
			if len(tmp[ex]) == 0: 	 eres[ex] = 'NA' 
			else: 			 eres[ex] = str(tmp[ex][0]) 

		notes['MEGAFIRE'] = K[s]['MF']['MEGAF'] + K[s]['FIRING_TYPE']['MF'] 
		for tv in ['MATURE','SUB-GRAN','CR']: notes[tv] = K[s][tv]['YES']   

		notes['CR'] += K[s]['ax-type']['CR']+K[s]['cr-type']['CRMZ'] 
		notes['SUB-GRAN'] +=  K[s]['cr-type']['SG']  	


		sn = [x[0] for x in ms['SIZE']] 
		ln = [x[0] for x in ms['LOC']] 
		zn = [x[0] for x in ms['LOC'] if x[0] not in ['UN','UNK']] 
		bn = [x[0] for x in ms['SIZE_LOC']] 

		dec = dd(float) 
		if len(sn) > 1: wt = ms['SIZE'][0][1] / float(ms['SIZE'][1][1])




		if s[0:2] == 'ES': dec['MEGAPHENO'] = -1.0 
		elif sn[0] == 'MEGA' and 'MINI' not in sn and '-' not in sn and 'UN' not in sn: 	dec['MEGAPHENO'] = 1.0
		elif sn[0] == 'CR' and 'MEGA' not in sn and ln == ['MZ']:	dec['MEGAPHENO'] = -1.0  
		elif sn[0] == 'MINI' and 'MEGA' not in sn and notes['MEGAFIRE'] == 0 and 'MZ' not in ln: dec['MEGAPHENO'] = -1.0
		elif sn[0] == 'MEGA' and eres['FIRING_TYPE'] == 'REPI' and wt >3: dec['MEGAPHENO'] = 0.75 
		elif sn[0] == 'MINI' and eres['FIRING_TYPE'] != 'REPI' and wt > 2:
			if 'MEGA' not in sn: dec['MEGAPHENO'] = -1.0 
			else: 		     dec['MEGAPHENO'] = -0.75 
		elif eres['FIRING_TYPE'] == 'NA':
			if 'MEGA' not in sn: dec['MEGAPHENO'] = -1.0
			elif sn[0] == 'MEGA':
				if wt > 2.5 and res['RSIZE'] > 2.0 and 'MZ' not in ln and ln != ['SP']: dec['MEGAPHENO'] = 0.75	
				elif zn == ['MZ'] or 'CSP' in zn: dec['MEGAPHENO'] = 0.0 
				elif zn == ['SP']: dec['MEGAPHENO'] = 0.5 
				elif wt > 1.5 and ('BIG' in sn or 'MED' in sn): dec['MEGAPHENO'] = 0.75 
				elif wt > 2.0 or res['RSIZE'] > 2.0 or 'BIG' in sn: dec['MEGAPHENO'] = 0.5
				else:						    dec['MEGAPHENO'] = 0.25 
			else:
				if 'BIG' in sn: dec['MEGAPHENO'] = 0.0 
				else:		dec['MEGAPHENO'] = -0.25 
		elif eres['FIRING_TYPE'] != 'REPI': 

			if zn == ['MZ']: 
				if sn[0] == 'MEGA': dec['MEGAPHENO'] = 0.25 
				else:		    dec['MEGAPHENO'] = -0.75 
			else:
				if sn[0] == 'MEGA':
					if wt >= 2: dec['MEGAPHENO'] = 0.50
					else:	    dec['MEGAPHENO'] = 0.25 

				else:		    dec['MEGAPHENO'] = -0.50 
		else:

			if 'CSP' not in zn and res['RSIZE'] > 1.5 and sn[0] == 'MEGA': dec['MEGAPHENO'] = 0.75 
			elif sn[0] == 'MEGA':  dec['MEGAPHENO'] = 0.25 
			elif 'MEGA' not in sn: dec['MEGAPHENO'] = -0.5 
			else: 		       dec['MEGAPHENO'] = -0.25 


		sl = [x for x in zn if x not in ['CSP','ISVZ']]
		gl = [x for x in zn if x not in sl]

		if len(zn) == 0: 
			dec['LOC'] = 'UNK' 
			dec['GLOC'] = 'UNK'
		else:
			if len(sl) > 0: dec['LOC'] = sl[0] 
			if len(gl) > 0: dec['GLOC'] = gl[0] 
			else:
				if sl[0] in ['MZ','ES']: dec['GLOC'] = sl[0] 
				elif sl[0] in ['SP','CP']: dec['GLOC'] = 'CSP'
				elif sl[0] in ['IZ','SVZ','IZ_SP','SVZ_IZ','SP_IZ_SVZ']: dec['GLOC'] = 'ISVZ'

		if len(bn) == 1: 
			dec['SLOC'] = bn[0] 
		else:
			if dec['LOC'] == 'ES': dec['SLOC'] = 'MINI-ESP'
			elif dec['LOC'] != 'MZ': 
				if dec['MEGAPHENO'] == -1.0: 
					dec['SLOC'] = 'MINI-'+dec['GLOC'] 	
				elif dec['MEGAPHENO'] == 1.0:
					dec['SLOC'] = 'MEGA-'+dec['GLOC'] 	
				else:
					dec['SLOC'] = 'UNSURE-'+dec['GLOC'] 	
			else:
					dec['SLOC'] = 'CR-MZ'

		if dec['SLOC'] == 'CR-MZ': dec['SLOC'] = 'CRET-MZONE'

		dec['PAIRED-LOC'] = dec['SLOC'] 
		dec['MEGA'] = dec['MEGAPHENO'] 

		notes['MFIRE'] = notes['MEGAFIRE'] 


		res['B-GRP'] = res['BIO-GROUP'] 		
		res['BIOP']  = res['BIO'] 

		res['REL-A'] = str(res['REL-SIZE'])[0:5] 
		res['AREA'] = str(res['AREA']).split('.')[0]

		res['TOT-READS'] =  ",".join([str(x)[0:5] for x in res['total']])
		try: res['TOT-GENES'] =  ",".join([str(int(x))[0:6] for x in res['obs']])
		except:  res['TOT-GENES'] = 'UNKNW,UNKNW'
		if res['TOT-READS'] == 'NA,NA':
			res['TOT-READS'] = 'UNKNW,UNKNW'
#		['GW', 'BIOP', 'AREA', 'B-GRP','REL-SIZE','READS', 'GENES','QC']


		edata = [eres['FIRING_TYPE'][0:4], str(eres['MAXSPIKES']), eres['SPIKE-GROUP']]
		#['MAXSPIKES', 'SPIKE-GROUP', 'FIRING_TYPE', 'MAX_LATENCY']



		print "\t".join([s]+[str(dec[k]) for k in SK]+['True ' if notes[k]>0 else 'False ' for k in NK]+[str(res[k]) for k in RK]+edata)
#		print notes,res,dec

#		print x,s,sorted(list(set(merge_key[s]['MEGA'])))

#		for x in Xs:
#			print x,s,sorted(list(set(merge_key[s][x])))













if __name__ == '__main__':

	from optparse import OptionParser
	usage = "usage: ./%prog [options] data_file"
	parser = OptionParser(usage=usage)


	parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
	parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
	parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

	(options, args) = parser.parse_args()

	run_script(args,options)	














