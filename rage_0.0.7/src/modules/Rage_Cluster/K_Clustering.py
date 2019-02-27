#!/usr/bin/env python

import random
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict as dd
from collections import Counter as cc
import sys
import os
import scipy.stats as stats
from scipy.stats import variation as coVar 

from random import random
import statsmodels.api as sm
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import fdrcorrection as fdr
import random
from math import fabs
#from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
import pickle
from math import log
import math
import numpy as np 
import pylab 
from matplotlib.patches import Rectangle as Rect
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib._png import read_png
from matplotlib.patches import Rectangle as Rect
from matplotlib.patches import Circle as Circ

from sklearn.decomposition import PCA

from sklearn.manifold import TSNE

from random import shuffle
				
from sklearn.cluster import KMeans	


import seaborn
from sklearn.cluster import KMeans

from sklearn.neighbors import KernelDensity

from sklearn.preprocessing import MinMaxScaler

from scipy.spatial import distance_matrix

def dist(d1,d2):
        
        return np.linalg.norm(d1-d2)*np.linalg.norm(d1-d2)
        print sum([(d1[i]-d2[i])*(d1[i]-d2[i]) for i in range(len(d1))])**0.5
        print np.linalg.norm(d1-d2)





def run_pcak(vals,my_feats,comps=2):

	fitMat = np.matrix(vals)
   	fitMat = fitMat - fitMat.mean(axis=0) 

	if len(vals) < 3 and comps == 0: 	return False, [[vals[0][n],vals[-1][n]] for n in range(len(vals[0]))], False , [1,1]
	elif len(vals) > 2:			v_run = PCA(n_components = max(comps,2)).fit(fitMat.getT())
	else:					v_run = PCA(n_components = 1).fit(fitMat.getT())
	
       	vPts = v_run.transform(fitMat.getT())
	v_coefs = sorted([(c*c,c,f) for (c,f) in zip(v_run.components_[0],my_feats)],reverse=True) # for comps in v_run.components_]
	var_rates = v_run.explained_variance_ratio_ 
	
	return v_run, vPts, v_coefs, var_rates+[var_rates[0]]



def scale_data(data):
	scaler = MinMaxScaler()
	d = np.array(data,dtype=float)


	s_vals =  [scaler.fit_transform(d[x].reshape(-1,1)).reshape(1,-1)[0] for x in range(len(d))]
	return s_vals
	return scaler.fit_transform(d.astype(float).reshape(d.shape[1],-1)).reshape(d.shape[0],-1) 

def scale_vals(vals):
	scaler = MinMaxScaler() 
	return scaler.fit_transform(np.array(vals,dtype=float).reshape(-1,1)).reshape(1,-1)[0] 


	
class K_Means_On_Genes:
        def __init__(self,data,fit=True):
		self.data      = data 	
		self.gene_vals = scale_data(data.vals) 

		xc =  [self.gene_vals[x][5] for x in range(len(self.gene_vals))]



		self.total_samples, self.total_genes = len(self.gene_vals[0]), len(self.gene_vals) 
		ko = self.data.args 
		self.k_start,self.k_end = [int(k) for k in ko.krange.split(',')]

		self.neighborRatio,self.cluster_ratio = float(ko.cluster_ratios.split(',')[0]),float(ko.cluster_ratios.split(',')[1])
		self.sample_outlier_rate = 1.0 - ko.outlier_rate


		self.min_cluster_size ={x: max(10,min((ko.min_cluster_size),int((self.total_samples/float(x))*self.data.args.min_cluster_ratio))) for x in range(self.k_start-1,self.k_end+1)}


		self.within_dist, self.between_dist = [float(k) for k in ko.cluster_dists.split(',')]
		#self.self.cluster_ratio = ko.cluster_ratio
		self.gene_results = dd(lambda: dd(bool))

		self.gene_results = {feat: dd(list) for feat in self.data.feats} 


	def run(self):
		for i,(feat,gene) in enumerate(zip(self.data.feats,self.gene_vals)):
			#if i > 3: break 
			s_vals = gene.reshape(-1,1)  
			self.test_1D(feat,s_vals)
		return self


	def test_1D(self,feat,vals): 

		self.km = [KMeans(n_clusters=p) for p in range(self.k_start,self.k_end)]
		self.run = [self.km[p].fit(vals) for p in range(len(self.km))]	
		k_centers, k_labels = [k.cluster_centers_ for k in self.run], [k.labels_ for k in self.run]
		for i,centers,labels in zip(range(len(self.run)),k_centers,k_labels):
			pass_clusters = [] 
			cluster_pass,cLen, center_key, val_key, sample_to_center, center_loc, center_samples  = [],i+self.k_start, dd(list), dd(lambda: dd(list)), dd(lambda: {}),{}, dd(list) 
			if sorted(cc(labels).values(), reverse=True)[1] < self.min_cluster_size[cLen]: continue 
                        for k,c in enumerate(centers):	center_key[k] = {m: dist(centers[m],c) for m in range(len(centers))}
			for j in range(len(vals)): 	center_samples[labels[j]].append(j) 

		
			pass_centers =  [c for (c,S) in center_samples.items() if len(S)*self.sample_outlier_rate >= self.min_cluster_size[cLen]]

#			if len(pass_centers) < 2: continue 

			for c,S in [(c,S) for (c,S) in center_samples.items() if len(S)*self.sample_outlier_rate >= self.min_cluster_size[cLen]]:
			
				print c,centers[c]
				s_major =  sorted([(dist(vals[j],centers[c]),j) for j in S])[0:int(len(S)*self.sample_outlier_rate)]
				s_items = [x[1] for x in s_major] 
				#print len(pass_centers) 

				s_scr =    sum([x[0] for x in s_major])/len(s_major) 
				s_out =    sorted([sum([dist(vals[j],centers[n]) for j in s_items])/len(s_major) for n in pass_centers if n != c])[0] 
				#print s_scr, i, len(S), centers[c], sorted([sum([dist(vals[j],centers[n]) for j in s_items])/len(s_major) for n in pass_centers if n != c])[0] 
				if s_out / s_scr <  self.neighborRatio: continue 
				#print 'cool'
				#print s_scr, i, len(S), centers[c], sorted([sum([dist(vals[j],centers[n]) for j in s_items])/len(s_major) for n in pass_centers if n != c])[0] 
				if s_scr < self.within_dist:	cluster_pass.append([s_scr,c,len(S),[x[1] for x in s_major]])
			

#			if feat.split(';')[1] not in ['CYP26A1','RELN']: continue 
			if len(cluster_pass) < 2: continue
			for x in range(len(cluster_pass)-1):
				x_scr,x_id,x_size, x_members = cluster_pass[x] 
				zDist,yIdx  = sorted([(center_key[x_id][cluster_pass[y][1]],y) for y in range(x+1,len(cluster_pass))],reverse=True)[0]
				y_scr,y_id,y_size,y_members = cluster_pass[yIdx]

				#print zDist,self.between_dist,zDist/x_scr, zDist/y_scr, self.cluster_ratio 
				#print zDist, min(zDist/x_scr,zDist/y_scr), self.between_dist, self.cluster_ratio
				if (zDist < self.between_dist) and  (min(zDist/x_scr,zDist/y_scr) < self.cluster_ratio): continue 
				x_types,y_types = [self.data.key[self.data.samples[nX]] for nX in x_members], [self.data.key[self.data.samples[nY]] for nY in y_members]
				dk = dd(lambda: [0,0]) 
				for a,b in cc(x_types).items():dk[a][0]=b 
				for a,b in cc(y_types).items():dk[a][1]=b 
				pass_clusters.append((x_id,x_members,[vals[x][0] for x in x_members]))
				pass_clusters.append((y_id,y_members,[vals[x][0] for x in y_members])) 
				self.gene_results[feat][cLen].append([(x_id,centers[x_id],x_size,x_scr),(y_id,centers[y_id],y_size,y_scr),zDist,dk])
			## GENERALIZE TO VOLUMES ## 
			if len(pass_clusters)>1:
				if len(centers[0]) == 1: 
					p_tmp =  [[vals[x][0] for x in center_samples[p[0]]] for p in pass_clusters]
					p_spots = sorted(list(set([(min(p),max(p)) for p in p_tmp])))
					if p_spots[0][1] == 0.0 and p_spots[1][0] > 0: 	p_spots[0] = (p_spots[0][0],p_spots[1][0]/2.0)
					all_spots = []
					p_volume = sum([p[1]-p[0] for p in p_spots])
					for p in range(len(p_spots)):
						if p == 0: 
							if p_spots[p][0] > 0:
								all_spots.append((0.0,p_spots[p][0]))
							all_spots.append(p_spots[p])
							if p_spots[p][1] < p_spots[p+1][0]:
								all_spots.append((p_spots[p][1],p_spots[p+1][0]))
							else:
								print 'wow'
								print p_spots
						 
						else:
							all_spots.append(p_spots[p]) 
						 	if p + 1 < len(p_spots):
								if p_spots[p][1] < p_spots[p+1][0]: 
									all_spots.append((p_spots[p][1],p_spots[p+1][0]))
								else:
									print 'weird'
									print p_spots
									print all_spots
									print "" 
							elif p_spots[p][1] < 1.0:
								all_spots.append((p_spots[p][1],1.001))

					val_locs = sorted([(vals[j][0],j) for j in range(len(vals))])
					k = 0 
					k_dict = dd(list) 
					k_grouped, k_total = 0, 0  
					for v,j in val_locs:
						v_data = [j,v,self.data.key[self.data.samples[j]]]

						while k < len(all_spots) and v > all_spots[k][1]: k+=1

						if v == all_spots[k][0] and v == 0:
							k_dict[all_spots[k]].append(v_data) 
							continue 
						if v >= all_spots[k][0] and v <= all_spots[k][1] and all_spots[k] in p_spots:
							k_dict[all_spots[k]].append(v_data) 
							k_grouped +=1

						elif v > all_spots[k][0] and v < all_spots[k][1] and all_spots[k] not in p_spots:
							k_dict[all_spots[k]].append(v_data)
							if v < p_spots[-1][0]:
								k_total += 1							
						elif v == all_spots[k][1]:
							k+=1 
							k_grouped+=1 
							k_dict[all_spots[k]].append(v_data) 	
						else:
							print 'wtf'
					pos_space =  sum([b-a for (a,b) in p_spots if a > 0])
					neg_space =  sum([b-a for (a,b) in all_spots if (a,b) not in p_spots and b <= p_spots[-1][0]])
					pos_density = k_grouped / pos_space
					neg_density = k_total/neg_space
					if neg_density > pos_density: continue 

					print feat,i,cLen,k_grouped,k_total,pos_density,neg_density
					for (a,b) in sorted(all_spots):
						if (a,b) in p_spots: lw = 3
						else: lw = 0 
						k_mems = k_dict[(a,b)]	
						k_height = len(k_mems) / (b-a) 
						if a == 0: k_height = len(k_mems) 
						kx = cc([km[-1] for km in k_mems]).items()
						hs = dd(int)
						for x,y in kx:
							if x[-1] == '+': hs[('r','BIG')] += y 
							elif x in ['TP','HP']: hs[('purple','AD')] += y 
							elif x in ['SVZ','IZ','SP']: hs[('blue','REG')]+=y 
							elif x == 'MZ': hs[('lime','MZ')] += y 
							elif x == 'ES': hs[('cyan','ES')] += y 
						#	else: 		hs[('k','NA')] += y 
						hx =  sorted(hs.items())
						bv = 0 
						for (c,t),lv in hx:
 							lh = (float(lv) / len(k_mems) ) * k_height 
							plt.bar(a,lh,b-a,bottom=bv,color=c,linewidth=lw)
							bv += lh  #bottom=base,color = color_key[name],linewidth=0)

						

					plt.title(feat+' '+str(cLen)) 
					plt.show()

	def write(self,out=sys.stdout):

		sys.exit() 	
		out.write('%-15s %4s %6s %8s %12s %6s %6s %8s %12s %6s ' % ('---','k#','c1_id','c1_loc','c1_scr','c1#','c2_id','c2_loc','c2_scr','c2#'))
		out.write(' | %10s | %9s %6s %4s %4s |\n' % ('dist','id','rate','c1#','c2#'))
	 	for f in self.gene_results.keys():
			for k in range(self.k_start,self.k_end):
				if len(self.gene_results[f][k]) == 0:
					continue 
					out.write('%-15s %4s  %4s\n' % (f,k,'None'))
				else:
#					for k_res in self.gene_results[f][k]: 
#						print k_res
					for x_data,y_data,z_dist,dk in self.gene_results[f][k]: 
						out.write('%-15s %4s ' % (f,k))
						#x_data,y_data,z_dist,dk = self.gene_results[f][k] 
						for kd in [x_data,y_data]:	
							loc = ",".join([str(round(s,2)) for s in kd[1]])
							if kd[3] < 0.01:	out.write('%6s %8s %12.3e %6s ' % ('c'+str(kd[0]),loc,kd[3],kd[2]))
							else:			out.write('%6s %8s %12.3f %6s ' % ('c'+str(kd[0]),loc,kd[3],kd[2]))
						out.write(' | %10.6f | ' % z_dist)
						id_data = [] 
						for x,y in dk.items(): 
							if y[0] > y[1]: id_data.append([float(y[0])/(y[0]+y[1]),y[0],y[1],x])
							else: 		id_data.append([float(y[1])/(y[0]+y[1]),y[0],y[1],x])
						id_data.sort(reverse=True) 
						for rate,x1,x2,name in id_data:
							out.write('%9s %6.3f %4d %4d |' % (name,rate,x1,x2))
						out.write('\n')



	
#					print dk.items() 
					
					#out.write('%s %s %s %s %s %s\n' % (f,k,x_data[1],x_data[2],x_data[3]))



		                        #k_result = K_Means_On_Genes(gene_data).run()
                        #kRange = range(int(args.krange.split(',')[0]),int(args.krange.split(',')[1]))
                        #for f in k_result.keys():
                         #       for k in kRange:
                          #              print f,k,k_result[f][k]










class K_Means_On_Matrix:
        #def __init__(self,samples,colors,ids,bw=0.2):
        def __init__(self,data,fit=True):

		scaler = MinMaxScaler()	
		self.data = data
		self.SAMPLES = True 
		self.SAMPLES = False 

		d = np.array(data.vals,dtype=float)

		if self.SAMPLES:
			self.vals = scaler.fit_transform(d.astype(float).reshape(d.shape[1],-1))
		else: 
			self.vals = scaler.fit_transform(d.astype(float).reshape(d.shape[1],-1)).reshape(d.shape[0],-1) 
		self.key = data.key
		self.samples = data.samples 

 		#sys.exit() 
		#v_run,vPts,v_coefs,v_vars = run_pcak(data.vals,data.feats,100) 
		#self.vals = vPts

#		print cc(self.data.key.values())

#		self.vals = d.reshape(d.shape[1],-1)

	def create_groups(self,test_num=19):

		self.group_cutoff = 0.1  

		### FIRST YOU GOTTA DO A SAMPLE PRUNING ### ---- ### THEN YOU CAN DO A FEATURE PRUNING ### 

		self.km = [KMeans(n_clusters=p) for p in range(1,test_num)]
		self.run = [self.km[p].fit(self.vals) for p in range(1,len(self.km))]	
		k_centers, k_labels = [k.cluster_centers_ for k in self.run], [k.labels_ for k in self.run]
		k_inertia = [k.inertia_ for k in self.run]

		
		for i,centers,labels in zip(range(len(self.run)),k_centers,k_labels):
			center_key, val_key, sample_to_center, center_loc, center_samples  = dd(list), dd(lambda: dd(list)), dd(lambda: {}),{}, dd(list) 
			cLen = i+2.0


			if max(sorted(cc(labels).values(),reverse=True)[1::])/(len(labels)/cLen) < 0.05: continue 
#			print max(sorted(cc(labels).values(),reverse=True)[1::])/(len(labels)/cLen),i+2	
#			print len(labels)/cLen,len(labels),cc(labels).values()  

			for k,c in enumerate(centers): center_loc[k] = c 
			for k,c in enumerate(centers): center_key[k] = sorted([(dist(centers[m],c),m) for m in range(len(centers)) if m != k])
			for j in range(len(self.vals)): 
				center_samples[labels[j]].append(j) 
				sample_to_center[j] = {k: dist(self.vals[j], c) for k,c in enumerate(centers)}
			
			for c,S in center_samples.items():
				s_dists = [(dist(self.vals[j],center_loc[c]),j) for j in S]
#				print sorted([x[0] for x in s_dists],reverse=True)[0:5]
#				print i,len(centers),centers
#				print 'to',len(s_dists) 

			
			for c,S in center_samples.items():
				if self.SAMPLES: 
					sIDs = [(self.key[self.samples[s]],self.samples[s]) for s in S]
					sTypes = sorted([(x[1],x[0]) for x in cc([s[0] for s in sIDs]).items()],reverse=True)
					print i+2.0,c,len(S),sTypes[0][1],sTypes[0][0]/float(len(S))
				else:
					sIDs = [(self.data.feats[s],s) for s in S]
					print i+2.0,c,len(S),sIDs
		 

		# ideas are to do condensation based distanaces, and do set parameters for the grouping, we need a min size rate for distance clusters, min support (number close) and no outliers, and consistency	

		sys.exit() 








