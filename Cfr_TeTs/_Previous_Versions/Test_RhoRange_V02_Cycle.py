#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 08:59:02 2025
@author: lsenni
Determine the RHO range where perform the averages:
    RHO-down and RHO-up
introduco d2: per indicare la mezza distanza dal centro plasma in centimetri, 
sulla quale mediare    
"""
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import time

import myelab_V12 as mye

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
tlim1 = d['tlim1']
tlim2 = d['tlim2']
eP = d['eP']
tTs = vars['tTs']  
errTs = vars['errTs']
tEce = vars['tEce'] 
psiTs = vars['psiTs']
psiEce = vars['psiKk1']
rhoTs = vars['rhoTs']
rhoEce = vars['rhoEce']

#d2 = 0.08 # mezzo intervallo intorno al centro- in centimetri
# messa sotto..poi spostare
#################### RHO
idt1 = (tTs.t >= tlim1) & (tTs.t <= tlim2) # Seleziono gli indici dei tempi di interesse
rhoTs2 = rhoTs[:,idt1] 
timeTs2 = tTs.t[idt1]
temp_ts2 = tTs.v[:,idt1]/1000 # in keV
err_ts2 = errTs.v[:,idt1]/1000  
err_ts2[err_ts2>2] = 1 # metto a 1 keV l'errore sui punti dove diverge

idt2 = (tEce.t>= tlim1) & (tEce.t <= tlim2)
rhoEce2 = rhoEce[:,idt2]
timeEce2 = tEce.t[idt2]
temp_ece2 = tEce.v[:,idt2]/1000  # in keV
err_ece2 = temp_ece2*eP

print('Dimensioni nell intervallo temporale di rhoTs = ',rhoTs2.shape)
print('Dimensioni nell intervallo temporale di rhoEce = ',rhoEce2.shape)
 
fig, ax = plt.subplots(2)
ax[0].plot(timeEce2,rhoEce2[40,:], label='RHO ECE')
ax[0].legend()
ax[0].set_title('Check Plot')
ax[1].plot(timeTs2,rhoTs2[20,:], label='RHO TS')
ax[1].legend()

########################
# Interpolation: make TS time sampling equal to ECE
# only in the selcted time interval

dimTs = rhoTs2.shape[1]    # number of Time points TS
dimEce = rhoEce2.shape[1]  # number of Time points TS

# New matrix with desired dimensions: same columns(time) as ECE
rhoTs2i = np.zeros((rhoTs2.shape[0], dimEce))
temp_ts2i = np.zeros((rhoTs2.shape[0], dimEce))
err_ts2i = np.zeros((rhoTs2.shape[0], dimEce))

# Per ogni posizione radiale (riga), interpoliamo rhoTs sui tempi di rhoEce
for i in range(rhoTs2.shape[0]):
    # Interpolazione dei dati di rhoTs ai tempi di rhoEce
    interpolator = interpolate.interp1d(timeTs2, rhoTs2[i, :], kind='cubic', fill_value="extrapolate")
    interpolator2 = interpolate.interp1d(timeTs2, temp_ts2[i, :], kind='cubic', fill_value="extrapolate")
    interpolator3 = interpolate.interp1d(timeTs2, err_ts2[i, :], kind='cubic', fill_value="extrapolate")
    rhoTs2i[i, :] = interpolator(timeEce2)        # Coord Rho del TS dopo interpolazione temporale tra t1 e t2
    temp_ts2i[i, :] = interpolator2(timeEce2)          # Temp del TS dopo interpolazione temporale tra t1 e t2
    err_ts2i[i, :] = interpolator3(timeEce2)

# Ora rhoTs_2 ha lo stesso numero di colonne di rhoEce
print("Dimensioni originali di rhoTs2:", rhoTs2.shape)
print("Dimensioni interpolazione di rhoTs2i:", rhoTs2i.shape)
print("Dimensioni di rhoEce:", rhoEce2.shape)

fig01,ax = plt.subplots(nrows =2)
ax[0].plot(timeEce2,rhoEce2[20,:], label='rhoEce2')
ax[0].set_title('Check RHOs in t1-t2 time interval')
ax[0].legend()
ax[1].plot(timeEce2,rhoTs2i[20,:], label='rhoTs2')
ax[1].legend()
#######################################
rTs = tTs.r
rEce = tEce.r

fig01, ax = plt.subplots(1,sharex = True)
ax.plot(rEce,rhoEce2[:,60],'-o', ms=2, label='rho ece')
ax.plot(rTs,rhoTs2i[:,60], '-o', ms=2, label='rho ts interpo')  
ax.set_title('Rho profiles at t=XX')
ax.legend() 

fig02, ax = plt.subplots(1,sharex = True)
ax.plot(timeEce2,temp_ts2i[20,:],'-', ms=2, label='Temp TSi')
ax.plot(timeEce2,temp_ece2[20,:], '-', ms=2, label='Temp Ece')  
ax.set_title('Time trace of Temp in different positions - DONT TRUST IT! ')
ax.legend() 

#######################################
# Inserisco ciclo per il calcolo dell'intervallo a tutti i tempi
# Rho range estimation - Consider only the ECE time scale

# for i in range(timeEce2):
#     mEce,pos_m = np.min(rhoEce2[:,i]) # min rho Ece at the ith indice
#     mTs = np.min(rhoTs2i[:,i])  # min rho Ts at the ith indice
d2 = 0.08 # mezzo intervallo intorno al centro- in centimetri
ranges = np.zeros([timeEce2.shape[0],2])  # inizializzo la matrice dei range in rho, uno per istante

temp_tsM_rho = np.zeros([rTs.shape[0],timeEce2.shape[0]])   # temp media in inttervallo rho
err_tsM_rho =np.zeros([rTs.shape[0],timeEce2.shape[0]])   # errore medio in intervallo rho

for i,_  in enumerate(timeEce2):
    mEce = np.nanmin(rhoEce2[:,i])   # min rho Ece at the i-th indice
    indE = np.nanargmin(rhoEce2[:,i]) # indice del minimo sopra calcolato
    posE = rEce[indE]                  # Position of the rho ECE minimum val
    mTs= np.nanmin(rhoTs2i[:,i])     # min rho Ts at the ith indice    
    indT = np.nanargmin(rhoTs2i[:,i])   
    posT = rTs[indT]                 # Position of the rho TS minimum val
    
    M = np.max([mEce,mTs])        # Trovo il massimo tra i due minimi
    ind = np.argmax(([mEce,mTs])) # per definire la diagnostica sulla Los della quale si trova il max:
                                  # ind = 0 il masssimo dei due minimi appartiene all'ECE...
    diag = [rhoEce2, rhoTs2i]   # Vettore con le due posizioni dei minimi 
    rDiag = [rEce,rTs]          # Costriusco un vettore con i due vettori delle posizioni
    pos = [posE,posT]
    # RHo - UP : calcolo d2 cm in meno e in più
    # rispetto alla poszione del massimo dei due minimi
    # e vedo il punto più vicino della medesima curva
    # Come RHO-UP prendo poi il massimo dei due valori trovati + 0.001
    
    diagM= diag[ind]
    r = rDiag[ind]
    
    pos1 = pos[ind]- d2
    pos2 = pos[ind] + d2
    
    idt1 = np.nanargmin(abs(r-pos1)) # indice1 più vicino sulla curva più alta
    idt2 = np.nanargmin(abs(r-pos2))
    # vettore dei valori più vicino alla posizione estrema dell'intervallo
    pippo = [v1,v2] = [rhoEce2[idt1,60],rhoEce2[idt2,60]] 
    # v2 = rhoEce2[idt2,60]
     
    siaM = np.max(pippo)
     
    rhoup = siaM + 0.0001
     #############################
     # RHO down
     
    # Trovo il punto più vicino della diag non considerata prima
    # seleziono l'altra diagnostica: minimo rta i due massimi
    indx = np.argmin(([mEce,mTs]))   # Trovo il minimo tra i due minimi
    
    diagm = diag[indx]      # seleziono la diagnostica 'piu bassa'
    id_pluto = np.nanargmin(abs(diagm[:,i]-M)) # trovo l'indice il punto con rho piu vicina
    pluto = diagm[id_pluto,i]
    rm = rDiag[indx]     # posizioni radiali della diagostica 'minore'
    pos_pluto = rm[id_pluto]
    
    # se il punto è prima del minimo aggiungo due puntio per determinare la rho down
    # se invece è dopo, ne tolgo due
    # in ambedue i casi tolgo poi un piccolo margine
    if id_pluto < indT:
        rhodown = diagm[id_pluto+2,i] - 0.001
    else:
        rhodown = diagm[id_pluto-2,60] - 0.001
    ranges[i,:] = [rhodown,rhoup]       
    
    
    for t in range(0,(rhoTs2i.shape[0])):
        mask = (rhoTs2i[:,t] >= rhodown) & (rhoTs2i[:,t] <= rhoup)    # Seleziono gli indici dei tempi di interesse 
        temp = temp_ts2i[mask,t]                    # Valori di temperaatura nelle rho selezionate
        sig = err_ts2i[mask,t]       # Errore sui singoli valori di temperatura: sigma
        ai = 1/sig**2            # Inverso del quadrato delle sigma
        somma = sum(ai)
        if somma==0:
            somma=1
        xm_ = sum(ai*temp)/somma   # Media dei valori di temperatura pesata sui singoli errori
        temp_tsM_rho[t,i] = (xm_)               # Valore delle temperature 'attese'/più probabili, nel tempo
        err_xm_ = np.sqrt(1/somma) 
        err_tsM_rho[t,i] = (err_xm_)
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 