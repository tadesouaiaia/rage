# 2019.02.06 16:19:20 EST
#Embedded file name: /home/tade/RAGE_PACKAGE/code/rage_0.0.7/rage_0.0.7/src/modules/Rage_Regression/rage_regmodels.py
from collections import Counter as cc
from collections import defaultdict as dd
import random
from sklearn.preprocessing import MinMaxScaler
from statsmodels.stats import power as smp
from statsmodels.stats.multitest import fdrcorrection as fdr
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif
import math
import numpy as np
import os
import sys
import statsmodels.stats.multitest as mlt
import rage_regtests
import warnings
warnings.filterwarnings('ignore')

def scale_vals(vals, f1 = 0, f2 = 1):
    scaler = MinMaxScaler(feature_range=(f1, f2))
    return scaler.fit_transform(np.array(vals, dtype=float).reshape(-1, 1)).reshape(1, -1)[0]


def vif_test(X):
    vd, vd_out = dd(list), {}
    for i, n in enumerate(X.names):
        try:
            vd_out[n] = round(vif(X.array, i), 4)
        except:
            vd_out[n] = 0.0

    return vd_out


class RegModel():
    def __init__(self, X, dist, options, progress = False, FULL = True):
        self.X, self.dist, self.options, self.progress, self.vif = X,dist,options,progress,vif_test(X) 
        self.alp1, self.alp2 = (0.05, 0.001)
        self.regress = rage_regtests.RegTest(X, self.dist, alphas=[self.alp1, self.alp2])
        self.out = dd(list)
        self.rs_key = [0.01,0.03,0.05,0.1,0.25]
        self.pv_key = [0.05,0.01,0.001,0.0001,1e-05,1e-07]
        self.resid = []

    def run(self, Y, Ynames = None):
        t_var, d_var, self.y_names = 0, 0, Ynames
        if self.dist[-3::].upper() == 'LOG':	self.Y = [ [ math.log(yc + 1.0, 2) for yc in y ] for y in Y ]
        else:	self.Y = Y
        if self.progress:	self.progress.start_minor('Running ' + self.dist + ' Regression', len(self.Y))
        for yi, y in enumerate(self.Y):
            if self.progress:	self.progress.mark()
            t = self.regress.test(y)
            self.out['params'].append(t.output)
            if len(self.X.predictor_idx) > 0:	self.out['pvs'].append([ p[0] for p in sorted(t.output) if p[-1] ][0])
            for r, k in zip([t.valid,t.history,t.zero_infl,t.v_explained,t.rsq,t.rsa,t.bic,t.pwr[self.alp1],t.pwr[self.alp2]], ['VALID','history','zero_infl','v_exp','rsq','rsa','bic','pwr1','pwr2']): self.out[k].append(r)

            self.out['run'].append(t.res)

        fdrs = [ mlt.fdrcorrection([ self.out['params'][i][n][0] for i in range(len(self.Y)) ])[1] for n in range(len(self.X.names)) ]
        self.results = [ [ [fdrs[j][i]] + list(self.out['params'][i][j]) for j in range(len(fdrs)) ] for i in range(len(self.Y)) ]
        self.param_key = ['fdr',
         'pval',
         'tval',
         'bw',
         'name',
         'pred_bool']
        return self

    def aggregate(self, PARAMS = False):
        if len(self.X.predictor_idx) > 0: pv_list = [ sorted([ p[0] for p in P if p[-1] ])[0] for P in self.out['params'] ]
        else:				  pv_list = []
        self.pv_cnt = [ len([ p for p in pv_list if p < self.pv_key[j] ]) for j in range(len(self.pv_key)) ]
        self.rs_cnt = [ len([ p for p in self.out['rsq'] if p > self.rs_key[j] ]) for j in range(len(self.rs_key)) ]
        if PARAMS:	self.pv_dict = {self.X.names[i]:[ p[i][0] for p in self.out['params'] ] for i in range(len(self.X.names))}
        return self

    def get_resids(self, COVAR = True, SCALE = True):
        self.resids, self.p_resids, self.c_resids = [], [], []
        rIdx, cIdx, pIdx = range(len(self.X.names)), [ i for i in range(len(self.X.names)) if self.X.COVARIATE[self.X.names[i]] ], [ i for i in range(len(self.X.names)) if self.X.PREDICTOR[self.X.names[i]] ]
        if self.dist[0:3].upper() == 'OLS':
            for i, y in enumerate(self.Y):
                print y
                sys.exit()
                mR = [ y[s] - sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in rIdx ]) for s in range(len(y)) ]
                cR = [ y[s] - sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in cIdx ]) for s in range(len(y)) ]
                if SCALE:
                    sc = MinMaxScaler(feature_range=(0, max(y)))
                    mR = [ x[0] for x in sc.fit_transform(np.array(mR).reshape(-1, 1)) ]
                    cR = [ x[0] for x in sc.fit_transform(np.array(cR).reshape(-1, 1)) ]
                self.resids.append(mR)
                self.c_resids.append(cR)

            return (self.resids, self.c_resids)
        SCALE = True
        for i, y in enumerate(self.Y):
            if not self.out['VALID'][i]:
                self.resids.append(y)
                self.p_resids.append(y)
                self.c_resids.append(y)
                continue
            r_pred = [ np.exp(sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in rIdx ])) for s in range(len(y)) ]
            p_pred = [ np.exp(sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in pIdx ])) for s in range(len(y)) ]
            c_pred = [ np.exp(sum([ self.out['params'][i][j][1] * self.X.array[s][j] for j in cIdx ])) for s in range(len(y)) ]
            if SCALE:
                try:
                    sc = MinMaxScaler(feature_range=(0, max(y)))
                    r_pred = [ x[0] for x in sc.fit_transform(np.array(r_pred).reshape(-1, 1)) ]
                    p_pred = [ x[0] for x in sc.fit_transform(np.array(p_pred).reshape(-1, 1)) ]
                    c_pred = [ x[0] for x in sc.fit_transform(np.array(c_pred).reshape(-1, 1)) ]
                except ValueError:
                    print 'FAIL'
                    print self.y_names[i]
                    print max(y), min(y)
                    sys.exit()

            r_resid = [ ys - r_pred[s] for s, ys in enumerate(y) ]
            p_resid = [ ys - r_pred[s] for s, ys in enumerate(y) ]
            c_resid = [ ys - r_pred[s] for s, ys in enumerate(y) ]
            if min(r_resid) < 0:
                r_resid = [ rr + -1 * min(r_resid) for rr in r_resid ]
            if min(p_resid) < 0:
                p_resid = [ rr + -1 * min(p_resid) for rr in p_resid ]
            if min(c_resid) < 0:
                c_resid = [ rr + -1 * min(p_resid) for rr in c_resid ]
            self.resids.append(r_resid)
            self.p_resids.append(p_resid)
            self.c_resids.append(c_resid)

        return (self.resids, self.c_resids)

    def run_permutation_tests(self, V):
        maxT = self.options.maxtrials
        if maxT < 100:
            maxT = 100
        if maxT > 1000:
            minT = 100
            midT = 200
        elif maxT > 499:
            minT = 50
            midT = 250
        else:
            minT = 10
            midT = 50
        minCut, midCut = 0.05 * minT, 0.01 * (midT + minT)
        self.permutations = dd(list)
        for i, y in enumerate(self.Y):
            p_cands, self.p_key = [ [a[3], a[0]] for a in [ r for r in self.out['params'][i] if r[-1] ] ], {}
            for c, pv in p_cands:
                if pv > 0.1:
                    self.permutations[c].append((round(pv, 2), 1))
                else:
                    self.p_key[c] = [pv, 0, 0]

            if len(self.p_key) > 0:
                self.regress.permute(y, [ V.select_variables(V.variables, permute=V.predictors) for j in range(minT) ], self.p_key)
                for c, (pv, L, G) in self.p_key.items():
                    if L > minCut:
                        self.permutations[c].append((self.p_key.pop(c)[1] / float(L + G), L + G))

                if len(self.p_key) > 0:
                    self.regress.permute(y, [ V.select_variables(V.variables, permute=V.predictors) for j in range(midT) ], self.p_key)
                    for c, (pv, L, G) in self.p_key.items():
                        if L > midCut:
                            self.permutations[c].append((self.p_key.pop(c)[1] / float(L + G), L + G))

                    if len(self.p_key) > 0:
                        self.regress.permute(y, [ V.select_variables(V.variables, permute=V.predictors) for j in range(maxT) ], self.p_key)
                        for c, (pv, L, G) in self.p_key.items():
                            self.permutations[c].append((self.p_key.pop(c)[1] / float(L + G), L + G))


