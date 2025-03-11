#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 21:45:26 2023
author: egio
Modified 09/10/2023 - lsenni
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt

shot1 = 104522
w = ppfs(shot1)
R = 3.0
a = 1.0
Bt = 3.85
#Esempio JET
R = 3.0
a = 1.0
Bt = 0.97
ne = 2.6e19
Ti = 6000

def rhostar(Ti, Bt, a):  # Normalized Larmor Radius
    return 2.04e-4*np.sqrt(Ti)/a/Bt

def nustar(ne, R, q, Ti, eps):  # Collisionality
    return 520*ne*1e-19*R*q/eps**1.5/Ti**2

def betan(rhos, nus):  # Normalized Beta: Betan
    return 550*rhos**1.02 * nus**0.19

_,Ti = w.hrts.te.get(r=0)
t,ne = w.hrts.ne.get(r=0)
tq, q95 = w.efit.q95.get(r=0)
_, Bvac = w.efit.bvac.get(r=0)
q = np.interp(t, tq, q95)
Bt = -np.interp(t, tq, Bvac)
#q = 3
eps = a/R
Te = w.hrts.te

rhos = rhostar(Ti, Bt, a)
nus = nustar(ne, R, q, Ti, eps)
betn = betan(rhos, nus)

plt.figure('Ciao', clear=True)
plt.plot(t, betn,label='BETAn')
w.efit.btnm.plot(label='EFIT')
plt.legend()

plt.figure('Te')
plt.plot(Te.t,Te.v[1,:],label='Te - HRTS')
plt.xlabel('time(s)')
plt.ylabel('Electron temperature (keV)')
plt.title(shot1)
plt.legend()