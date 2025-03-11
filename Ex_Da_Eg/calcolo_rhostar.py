#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 21:45:26 2023

@author: egio
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt

def rhostar(Ti, Bt, a):
    return 2.04e-4*np.sqrt(Ti)/a/Bt

def nustar(ne, R, q, Ti, eps):
    return 520*ne*1e-19*R*q/eps**1.5/Ti**2

def betan(rhos, nus):
    return 550*rhos**1.02 * nus**0.19

w = ppfs(104522)
R = 3.0
a = 1.0
Bt = 3.85
#Esempio JET
R = 3.0
a = 1.0
Bt = 0.97
ne = 2.6e19
Ti = 6000

_,Ti = w.hrts.te.get(r=0)
t,ne = w.hrts.ne.get(r=0)
tq, q95 = w.efit.q95.get(r=0)
_, Bvac = w.efit.bvac.get(r=0)
q = np.interp(t, tq, q95)
Bt = -np.interp(t, tq, Bvac)
#q = 3
eps = a/R

rhos = rhostar(Ti, Bt, a)
nus = nustar(ne, R, q, Ti, eps)
betn = betan(rhos, nus)

plt.figure('Ciao', clear=True)
plt.plot(t, betn)
w.efit.btnm.plot()
plt.legend(betan)

