#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:01:05 2024

@author: lsenni
"""
from ppfeg import ppfs
import matplotlib.pyplot as plt
import numpy as np
import my_flush

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse

# Creo dizionario dei parametri
d = {
    'shot' : 104522,  # 99950, 99971
    'tlim1' : 47.5,     # limite inferiore selezione tempi - in secondi
    'tlim2' : 52.5, 
    'delta' : 3,      # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.0,      # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.008,   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.015,   
    'rho1' : 0.01,    
    'rho2' : 0.15,
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : '/home/lsenni/Python_LS/Cfr_TeTs/prova/'  # folder to save plot
    }    

shot = d['shot']     
tlim1 = d['tlim1']   
tlim2 = d['tlim2'] 
w = ppfs(shot) 

vars = {
    'tlim' : (tlim1+tlim2)/2,          # tempo medio al quale vengono calcolati i profili per plot di controllo
    'w' : w,                           # canale dati del ppfeg con tutti i dati relativi allo sparo in oggetto
    'tTs' : w.hrts.te,                 # Chann Te HRTS  
    'errTs' : w.hrts.dte,              # Chann Errors HRTS.dte
    'tEce' : w.ecm1.prfl,              # Chann Te Ece-Kk1
    'errEce' : (w.ecm1.prfl)*d['eP'],  # Error associated to the Kk1 meas: ab 2%
    'psiTs' : w.hrts.psi,              # Chann PSI HRTS 
    'psiKk1' : None }   

# vars['tlim'] = (tlim1+tlim2)/2  

print('JPN = ', shot)
print('t lim1 = ', tlim1)
print('t lim2 = ', tlim2)

# mye.psicalc(d,vars)  
# mye.rhocalc(d, vars)
##################################################   
tlim1 = d['tlim1']
tlim2 = d['tlim2']
delta = d['delta']
w= vars['w']

tTs = vars['tTs']       # Chan. Te HRTS  
psiTs = vars['psiTs']    # Channel psi hrts
tEce = vars['tEce']
####################
psiTs = w.hrts.psi
idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
time_ts = tTs.t[idt]
temp_ts = tTs.v[:,idt]
psi_ts = psiTs.v[:,idt]
tlim = vars['tlim']
######################################################

# Prove di calcolo e plots
tEce = w.ecm1.prfl 

zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista (in metri), 
zTs = w.hrts.z

# la seconda è l'intersezione con l'asse y. 'Position of antenna aperture centre (z,majR) '
rEce = tEce.r
rTs = tTs.r

z = np.full_like(rEce, zKk1) # Retta linea di vista parallela asse r
# zz = np.full_like(zTs) # Retta linea di vista parallela asse r

timeEce = tEce.t
timeTs = tTs.t

psiKk1 = np.zeros(tEce.v.shape)
psiTh = np.zeros(tTs.v.shape)

for i,time in enumerate(timeTs):
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi, _ = my_flush.Flush_getFlux(rTs*100, zTs.v[:,0]*100) # Le coordinate vanne messe in cm!
    psiTh[:,i] = psi
   
# vars['psiKk1'] = psiKk1      

idts = np.argmin(abs(timeTs - tlim))
psiTs_s = np.absolute(psiTs.v[:,idts]) # Abs del profilo psi hrts al tempo t=tlim  
psiTh_s = np.absolute(psiTs.v[:,idts]) # Abs del profilo psi hrts al tempo t=tlim  
    
fig01,ax01 = plt.subplots(nrows=1, num = f'JPN {shot} PSI profiles')
ax01.plot(rTs, psiTs_s, linewidth=0.5, color='green', label= 'psi hrts')
ax01.plot(rTs, psiTh_s, linewidth=0.5, color='blue', label='psi ece')
ax01.set_ylabel('PSI')
ax01.legend()
# ax01.axhline(y = psi1, c='r',ls='--',lw=.4)
# ax01.axhline(y = psi2,c='r',ls='--',lw=.4)
# ax01.plot(rTs, psiTs_s, linewidth=0.5, color='green', marker='o', ms=0.8, label= 'psi hrts')
# ax01.plot(rTs[idxPsiTs], psiTs_s[idxPsiTs], linewidth=1.5, color='red', marker='o', ms=3, label=f'{psi1}<psi hrts<{psi2}')
# ax01.plot(rEce, psiEce_s, linewidth=0.5, color='blue', marker='*', ms=0.8, label='psi ece')
# ax01.plot(rEce[idxPsiE], psiEce_s[idxPsiE], linewidth=1.5, color='orange', marker='*', ms=3,  label=f'{psi1}<psi ece<{psi2}')
# ax01.axhline(y = psi1, c='r',ls='--',lw=.4)
# ax01.axhline(y = psi2,c='r',ls='--',lw=.4)
# ax01.set_xlim(left = leftlim-0.1, right = rightlim+0.1)
# ax01.set_ylim(bottom = psi1-0.01, top = psi2 + 0.01)
ax01.set_xlabel('R (m)')
ax01.set_ylabel('PSI')
ax01.legend()
ax01.set_title(f'JPN {shot} PSI(R) for HRTS and ECE-KK1 at t={tlim} sec')
# fig01.tight_layout()
   
