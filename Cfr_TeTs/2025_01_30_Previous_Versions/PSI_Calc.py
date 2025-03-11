#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:00:09 2024

@author: lsenni
Calcolo e controllo degli andamenti di PSI sulle linee di vista
"""
import numpy as np
import my_flush
from ppfeg import ppfs
import matplotlib.pyplot as plt

plt.close('all')  

shot = 96994 
tim = 52


w = ppfs(shot)   

tTs = w.hrts.te       # Chan Te HRTS  
psiTs = w.hrts.psi    # Channel psi hrts

tEce = w.ecm1.prfl 
zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y
time_ece = tEce.t
time_ts = tTs.t
psiKk1 = np.zeros(tEce.v.shape)
radii = np.zeros(tEce.r.shape)

for i,time in enumerate(time_ece):
    r, te = tEce.get(t=time)  # Profilo a T=time
    z = np.full_like(r, zKk1) # Retta linea di vista parallela asse r
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
    psiKk1[:,i] = psi
    
# Calcolo del profilo a dato tempo 

rEce = r
idx= np.argmin(abs(time_ece - tim))
psi_ece = psiKk1[:,idx]

idx = np.argmin(abs(time_ts-tim))
psi_ts = psiTs.v[:,idx]
idd = psi_ts>=0
psi_ts = psi_ts[psi_ts>=0]

tempEce = tEce.slice(t=tim)
tempTs = tTs.slice(t=tim)
temp_ts = tempTs.v[idd]
rad_ts = psiTs.r[idd]

fig, ax = plt.subplots(nrows=1,  num = f'JPN {shot} PSI profiles')
ax.plot(rEce,psi_ece,label='psi ece')
ax.plot(rad_ts,psi_ts, label='psi hrts')
ax.set_xlabel('R (m)')
ax.set_ylabel('PSI')
ax.legend()
ax.set_title(f'JPN {shot} PSI(R) at t={tim} sec')

fig2,ax2 = plt.subplots(nrows=1, num = f'JPN {shot} Te vs PSI')
ax2.scatter(psi_ece, tempEce.v/1000, label='Te ece - kk1')
ax2.scatter(psi_ts,temp_ts/1000,label='Te hrts')
ax2.set_xlabel('PSI')
ax2.set_ylabel('Electron Temperature (keV)')
ax2.legend()
ax2.set_title(f'JPN {shot} Te(PSI) at t={tim} sec')


