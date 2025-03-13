#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 08:59:02 2025
@author: lsenni
Provo a calcolare a mano le coordinate di flusso anche per l'HRTS (dove è già disponibile il canale)
pechè sembrava venisse strana. 
Il calclo viene identico udando FLUSH, quindi tutto ok
NOTA: l'intervallo temporale di acquisizione del HRTS è molto più lungo
e questo causa appartenti differenze:
Il rate più veloce dell'ECE vinen 'nascosto' da questo fatto,
e quando di fa l'interpolazione non era chiaro il risultato '
"""
from ppfeg import ppfs
from scipy import signal,interpolate
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import numpy as np
import time

import myelab_V12 as mye
import my_flush

# import my_manage_file as mym
# %reset

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse
JPN = 104522 #104990 # 104522 # 104525 #104990  #96994
Tref = 6 # keV) Reference temperature: time instants with Te values above Tref are considered
timestr = time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation

d = {
    'shot' : JPN,     # 99950, 99971, 99970
    'Tref' : Tref,
    'tlim1' : 47,     # 47,     # limite inferiore selezione tempi - in secondi
    'tlim2' : 53,     # 53,     # limite superiore selezione tempi - in secondi
    'fix' : 0.01,    # intervallo aggiuntivo da considerare (in psi e rho) nella differenza tra le due
    'intPsi' : 0.1,    # intervallo in Psi e Rho su cui mediare 
    'delta' : 3,      # tempo in più di tlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.05,     # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.005,  # 0.001 Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.025,     # 0.027 #0.02
    'rho1' : 0.055,   #0.045,   # 0.0089,    
    'rho2' : 0.11,  # 0.106,     # 0.1,
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'win_len' : 15,   # window lenght: number of points for the smooth with Savitsky-golay filter
    'deg_pol' : 3,    # grado del polinomio usato per lo smooting
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : f'/home/lsenni/Python_LS/Cfr_TeTs/{timestr}_JPN_{JPN}_Plots/'  # folder to save plot
    }    

#######################################################
# Creo il dizionario con le variabili
vars = mye.dictvar(d)

print('JPN = ', d['shot'])
print('t lim1 = ', d['tlim1'])
print('t lim2 = ', d['tlim2'])

#######################################################
mye.psicalc(d,vars)   # PSI: Normalized poloidal flux coordinate
mye.rhocalc(d, vars)  # RHO: Normalized toroidal flux coordinate

#######################################################
def psicalcTs(d, vars):
    shot = d['shot']
    w = vars['w']
    ####################
    # Prove di calcolo e plots
    tTs = w.hrts.te
    zTs = w.hrts.z.v[1,0]   # antp è il canale con le coordinate della linea di vista (in metri), 
    # la seconda è l'intersezione con l'asse y. 'Position of antenna aperture centre (z,majR) '
    rTs = tTs.r
    z = np.full_like(rTs, zTs) # Retta linea di vista parallela asse r
    
    timeTs = tTs.t

    psiTscalc = np.zeros(tTs.v.shape)
    
    for i,time in enumerate(timeTs):
        ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
        # e Prendo il tempo più vicino a quello desiderato
        psi, _ = my_flush.Flush_getFlux(rTs*100, z*100) # Le coordinate vanne messe in cm!
        psiTscalc[:,i] = psi
   
    vars['psiTscalc'] = psiTscalc      
    
    return 1
#######################################################
psicalcTs(d,vars)
#######################################################
# Function to compute the RHO limits
#def rhorange(d, vars): 
shot = d['shot']
tlim1 = d['tlim1']
tlim2 = d['tlim2']
tTs = vars['tTs']  
tEce = vars['tEce'] 
psiTscalc = vars['psiTscalc']
psiTs = vars['psiTs']
psiEce = vars['psiKk1']
rhoTs = vars['rhoTs']
rhoEce = vars['rhoEce']
####################
idt = (tTs.t >= tlim1) & (tTs.t <= tlim2) # Seleziono gli indici dei tempi di interesse
rhoTs2 = rhoTs[:,idt] 
timeTs2 = tTs.t[idt]

idt = (tEce.t>= tlim1) & (tEce.t <= tlim2)
rhoEce2 = rhoEce[:,idt]
timeEce2 = tEce.t[idt]
####################

fig, ax = plt.subplots(3)
ax[0].plot(tEce.t,psiEce[40,:], label=' Full PSI ECE')
ax[0].legend()
ax[1].plot(tTs.t,psiTs.v[20,:], label='Full PSI TS')
ax[1].legend()
ax[2].plot(tTs.t,psiTscalc[20,:], label='Full PSI TS CALC')
ax[2].legend()

fig, ax = plt.subplots(2)
ax[0].plot(tEce.t,rhoEce[40,:], label='Full RHO ECE')
ax[0].legend()
ax[1].plot(tTs.t,rhoTs[20,:], label=' Full RHO TS')
ax[1].legend()

 

####################
idt = (tTs.t >= tlim1) & (tTs.t <= tlim2) # Seleziono gli indici dei tempi di interesse
psiTs2 = psiTs.v[:,idt]
rhoTs2 = rhoTs[:,idt] 
timeTs2 = tTs.t[idt]

idt = (tEce.t>= tlim1) & (tEce.t <= tlim2)
psiEce2 = psiEce[:,idt]
rhoEce2 = rhoEce[:,idt]
timeEce2 = tEce.t[idt]

print('Dimensioni nell intervallo temporale di psiTs = ',psiTs2.shape)
print('Dimensioni nell intervallo temporale di psiEce = ',psiEce2.shape)
print('Dimensioni nell intervallo temporale di rhoTs = ',rhoTs2.shape)
print('Dimensioni nell intervallo temporale di rhoEce = ',rhoEce2.shape)
 
fig, ax = plt.subplots(2)
ax[0].plot(timeEce2,psiEce2[40,:], label='PSI ECE')
ax[0].legend()
ax[1].plot(timeTs2,psiTs2[20,:], label='PSI TS')
ax[1].legend()

fig, ax = plt.subplots(2)
ax[0].plot(timeEce2,rhoEce2[40,:], label='RHO ECE')
ax[0].legend()
ax[1].plot(timeTs2,rhoTs2[20,:], label='RHO TS')
ax[1].legend()


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 