#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd
from scipy import stats
from math import log 


fKey = {'EXONIC': 'EX', 'INTRONIC': 'IR','INTERGENIC': 'IT','KJXN': 'KJ','CJXN': 'CJ','EXONIC-OVERLAP': 'EXO','EXONIC-SPECIFIC': 'EXS'}


def revcomp(seq):
    C = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([C[x] for x in seq[-1::-1]])


def scan(cnts_file,options):
    genes = dd(list)

    minC = 8 
    minT = 4 
 
    for line in open(cnts_file):
        line = line.split()
        if line[0] == '---':
            names,totals = line[1::],[[] for i in range(len(line)-1)]

        else:



	    gene,Tcnts = line[0],[float(x) for x in line[1::]] 
	    if len([t for t in Tcnts if t > 0]) < minC: continue
	    if len([t for t in Tcnts if t > 2]) < minT: continue 
 	    cnts = [c for c in Tcnts if c != 'NA'] 
	    if len(cnts) < 5 or sum(cnts) == 0: continue
            mc = sum(cnts)/float(len(cnts))
            if mc == 0: continue
            obs = len([x for x in cnts if x != 'NA' and x>0])/float(len(cnts))
            std = (round((sum([(t-mc)*(t-mc) for t in cnts]) / float(len(cnts)-1))**0.5,3))
            genes[gene] = [round(mc,2),round(obs,2),round(std,2),round(std/mc,2),[log(c+1,2) for c in Tcnts]]

    gk = sorted(genes.keys())
    for i in range(len(gk)):
        mI,oI,sI,cvI,cntsI = genes[gk[i]]
        for j in range(i+1,len(gk)):
            mJ,oJ,sJ,cvJ,cntsJ = genes[gk[j]]
		
	    
            iOnly = [z for z in range(len(cntsI)) if cntsI[z] > 0 and cntsJ[z] == 0] 
            jOnly = [z for z in range(len(cntsI)) if cntsI[z] == 0 and cntsJ[z] > 0] 
            cNone = [z for z in range(len(cntsI)) if cntsI[z] == 0 and cntsJ[z] == 0]
            cBoth = [z for z in range(len(cntsI)) if cntsI[z] > 0 and cntsJ[z] > 0] 
             

	    iOne = [cntsI[z] for z in range(len(cntsI)) if (cntsI[z]+cntsJ[z]>0)] 
	    jOne = [cntsJ[z] for z in range(len(cntsI)) if (cntsI[z]+cntsJ[z]>0)] 
            #R=stats.pearsonr(cntsI,cntsJ)[0]
            #S=stats.spearmanr(cntsI,cntsJ)[0]
            R=stats.pearsonr(cntsI,cntsJ)[0]
            S=stats.spearmanr(cntsI,cntsJ)[0]
	    #if R*R > 0.5 or S*S > 0.5: 	
            R2=stats.pearsonr(iOne,jOne)[0]
            S2=stats.spearmanr(iOne,jOne)[0]


	    
            print gk[i],gk[j],"m/o/s/cv",mI,oI,sI,cvI,"|",mJ,oJ,sJ,cvJ,"| R/S",R,S
            print gk[i],gk[j],"Both,A,B,None",len(cBoth),len(iOnly),len(jOnly),len(cNone),'|',R2,S2












if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)


    parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
    parser.add_option("-k", "--key", default = None, type='string', help="Output Filename Prefix")

    (options, args) = parser.parse_args()



    scan(args[0],options)













