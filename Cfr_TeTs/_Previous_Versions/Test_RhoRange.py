#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 08:59:02 2025
@author: lsenni
Determine the RHO range where perform the averages:
    RHO-down and RHO-up
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
# Function to compute the RHO limits
#def rhorange(d, vars): 
shot = d['shot']
tlim = vars['tlim']
tTs = vars['tTs']  
tEce = vars['tEce'] 
rhoTs = vars['rhoTs']
rhoEce = vars['rhoEce']
####################
idts = np.argmin(abs(tTs.t - tlim))
rhoTs_s = rhoTs[:,idts] # Profilo rho hrts al tempo t=tlim   

ide = np.argmin(abs(tEce.t - tlim))
rhoEce_s = rhoEce[:,ide]     # profilo rho ece al tempo t=tlim  
rTs = tTs.r
rEce = tEce.r

rhoTs_s = rhoTs[:,idts] # Profilo rho hrts al tempo t=tlim   
rhoEce_s = rhoEce[:,ide] 
    
plt.figure()
plt.plot(rTs, rhoTs_s, linewidth=0.5, color='green', label='Rho Hrts LoS')
plt.plot(rEce, rhoEce_s, linewidth=0.5, color='blue', label='Rho Ece LoS')
plt.xlabel('Major Radius(m)')
plt.ylabel('Normalized RHO a t = XX')
plt.legend()

plt.figure()
plt.plot(tTs.t,rhoTs[20,:])
plt.xlabel('time(sec)')
plt.ylabel('Norm RHO at R= XX')
#plt.legend()
########################
# Interpolation: make TS time sampling equal to ECE

dimTs = rhoTs.shape[1]    # number of Time points TS
dimEce = rhoEce.shape[1]  # number of Time points TS

# I tempi originali di rhoTs e rhoEce
time_rhoTs = tTs.t    # Original time samplinf for TS 
time_rhoEce = tEce.t  # Original time samplinf for ECE

# New matrix with desired dimensions: same columns(time) as ECE
rhoTs_2 = np.zeros((rhoTs.shape[0], dimEce))

# Per ogni posizione radiale (riga), interpoliamo rhoTs sui tempi di rhoEce
for i in range(rhoTs.shape[0]):
    # Interpolazione dei dati di rhoTs ai tempi di rhoEce
    interpolator = interpolate.interp1d(time_rhoTs, rhoTs[i, :], kind='cubic', fill_value="extrapolate")
    rhoTs_2[i, :] = interpolator(time_rhoEce)

# Ora rhoTs_2 ha lo stesso numero di colonne di rhoEce
print("Dimensioni originali di rhoTs:", rhoTs.shape)
print("Dimensioni interpolazione di rhoTs_2:", rhoTs_2.shape)
print("Dimensioni di rhoEce:", rhoEce.shape)




plt.figure()
plt.plot(tEce.t,rhoTs_2[20,:])
plt.plot(tTs.t,rhoTs[20,:])
plt.xlabel('time(sec)')
plt.ylabel('Norm RHO at R= XX')
plt.legend()

fig01,ax0 = plt.subplots(nrows =2)
ax0.plot(tEce.t,rhoTs_2[20,:])
ax0.plot(tTs.t,rhoTs[20,:])

    
#######################################
#######################################





dimTs = rhoTs.shape[1]
dimEce = rhoEce.shape[1]
rhoTs_2 = np.zeros([rhoTs.shape[0],rhoEce.shape[1]])

# for i in rTs:
#     rhoTs_2[i,:]  = signal.resample_poly(rhoTs[i,:], up=dimEce, down=dimTs, axis=1) # posizioni ricampionate dei valori TeEce

pippo = rhoTs[30,:]
pluto = rhoEce[30,:]
fun_int = interpolate.make_interp_spline(pippo.shape[0], pluto.shape[0])


for i in (rTs):
    rhoTs_ = rhoTs[i,:]
    fl_int = interpolate.make_interp_spline(rhoTs, vFlux_/vFlux_[-1]) # function to interp the flux
    fl_int_hrts = fl_int(psiTs_, extrapolate=False) # Interpolation without extrapolating
    rhoTs[:,i] = np.sqrt(fl_int_hrts)


##
rhoEce_s = rhoEce[:,ide] 
rhoTs_2_s = rhoTs_2[:,200]
    

plt.figure()
plt.plot(tEce.t,rhoTs_2[20,:])
plt.xlabel('time(sec)')
plt.ylabel('Norm RHO at R= XX')
    
    