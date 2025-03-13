#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 10:05:34 2024

@author: lsenni
"""
from ppfeg import ppfs
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import myelab_V05 as mye

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
psiKk1 = vars['psiKk1']
tTs = vars['tTs']       # Chan. Te HRTS  
psiTs = vars['psiTs']    # Channel psi hrts
tEce = vars['tEce']
####################
psiTs = w.hrts.psi
idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
time_ts = tTs.t[idt]
temp_ts = tTs.v[:,idt]
psi_ts = psiTs.v[:,idt]

torFlux = w.efit.ftor
rFlux = w.efit.ftor.r  # Posizioni dei valori di flusso
vFlux = w.efit.ftor.v  # valori di flusso

rhoTs = np.zeros(tTs.v.shape)
rhoEce = np.zeros(psiKk1.shape)
 
for i,time in enumerate(tTs.t):
    psiTs_ = np.abs(psiTs.v[:,i])
    index = np.argmin(abs(torFlux.t - time))
    vFlux_ = vFlux[:,index]
    fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1]) # function to intertp the flux
    fl_int_hrts = fl_int(psiTs_, extrapolate=False) # Interpolation without extrapolating
    rhoTs[:,i] = np.sqrt(fl_int_hrts)
    
for i,time in enumerate(tEce.t):
    psi = psiKk1[:,i]
    index = np.argmin(abs(torFlux.t - time))
    vFlux_ = vFlux[:,index]
    fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1])  # function to intertp the flux
    fl_int_ece = fl_int(psi, extrapolate=False)
    rhoEce[:,i] = np.sqrt(fl_int_ece)
    # rho, _ = Flush.Flush_getFluxLabelRho(psi)
    # rhoKk1[:,i] = rho
 
# vars['time_ts'] = time_ts   # time interval t1-delta<t<t2+delta   
# vars['rho_ts'] = rho_ts

##################################################   

shot = d['shot']
rho1 = d['rho1']
rho2 = d['rho2']
tlim = vars['tlim']
tTs = vars['tTs']  
tEce = vars['tEce'] 
####################
    

# Seleziono gli indici corrispondenti all'intervallo in rho: rho1-rho2
idts = np.argmin(abs(tTs.t - tlim))

# rho mio: rhoMio
rhoTs_s = rhoTs[:,idts]  # Profilo rho hrts al tempo t=tlim   
rTs = tTs.r
selected_time_mio = tTs.t[idts]
print('selected_time_mio:',selected_time_mio)

# Rho Flush: rhoPPF
rhoPPF = vars['w'].hrts.rho.slice(t=tlim)
rhoPPF = np.abs(rhoPPF.v)

# Rho flush idts: rhoFlush
rhoPpf = vars['w'].hrts.rho
idx = np.argmin(abs(rhoPpf.t - tlim))
rhoPpf_s = rhoPpf.v[:,idx]
rhoFlush = np.abs(rhoPpf_s)
selected_time_slice = rhoPpf.t[idx]
print('selected_time_slice:',selected_time_slice)

plt.figure()
plt.plot(rTs, rhoTs_s, label='RHO mio')
plt.plot(rTs, rhoPPF, label='RHO FLUSH')
plt.xlabel('R(m)')
plt.ylabel('RHO')
plt.legend()
plt.title(f'JPN {shot} - Check plot: RHO-tor norm vs Nomalized minor radius')  

