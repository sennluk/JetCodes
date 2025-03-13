#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 08:59:02 2025
@author: lsenni
Determine the RHO range where perform the averages:
    RHO-down and RHO-up
introduco d2: per indicare la mezza distanza dal centro plasma in centimetri, 
sulla quale mediare    
V01: provo a togliere l'interpolazione e a ciclare sui tempi del TS (più lento)
prendendo per l'Ece il punto temporalmente più vicino. 
N.B.: l'Ece è molto più veloce, qiuindi ci sarà sempre un punto vicino temporalmente'
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
    'eP' : 0.03,      # Relative error assigned to the Ece data
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
# Inizio selezionando l'interv allo temporale
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

dimTs = rhoTs2.shape[1]    # number of Time points TS
dimEce = rhoEce2.shape[1]  # number of Time points TS
rTs = tTs.r
rEce = tEce.r

#######################################
# Inserisco ciclo per il calcolo dell'intervallo di RHo in tutti gli istanti di TS tra t1-t2
# Rho range estimation - Consider only the ECE time scale

d2 = 0.08 # mezzo intervallo intorno al centro- in centimetri
ranges = np.zeros([timeTs2.shape[0],2])  # inizializzo la matrice dei range in rho, uno per istante

temp_tsM_rho = np.zeros_like(timeTs2)      #inizializzo matrice Tts valori temp in t1-t2
err_tsM_rho = np.zeros_like(timeTs2)   #inizializzo matrice errTs valori Err-temp in t1-t2 

timeEce22 = np.zeros_like(timeTs2)

# valori tra t1 e t2 più vicini a quelli di TS
rhoEce22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])
temp_ece22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])
err_ece22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])

# inizializzo le mettrici delle medie, una per ogni tempo tra t1 e t2 basandosi sul campionamenti di time-TS
temp_eceM_rho = np.zeros_like(timeTs2) 
rho_eceM_rho = np.zeros_like(timeTs2)
err_eceM_rho = np.zeros_like(timeTs2) 


for i,inst  in enumerate(timeTs2):
    timeTs = inst                         # valore tempo per TS
    timeEce = np.min(abs(timeEce2-inst))  # valore tempo per Ece pià vicino
    iEce = np.argmin(abs(timeEce2-inst))  # valore indice del tempo Ece
    timeEce22[i] = timeEce2[iEce] 
    rhoEce22[:,i] = rhoEce2[:,iEce] 
    temp_ece22[:,i] = temp_ece2[:,iEce] 
    err_ece22[:,i] = err_ece2[:,iEce] 
    
    mTs= np.nanmin(rhoTs2[:,i])     # min rho Ts at the ith indice    
    indT = np.nanargmin(rhoTs2[:,i])   
    posT = rTs[indT]  
    
    mEce = np.nanmin(rhoEce2[:,iEce])   # min rho Ece at the i-th indice
    indE = np.nanargmin(rhoEce2[:,iEce]) # indice del minimo sopra calcolato
    posE = rEce[indE]                  # Position of the rho ECE minimum val
    
    M = np.max([mEce,mTs])        # Trovo il massimo tra i due minimi
    ind = np.argmax(([mEce,mTs])) # per definire la diagnostica sulla Los della quale si trova il max:
                                  # ind = 0 il masssimo dei due minimi appartiene all'ECE...
    diag = [rhoEce2, rhoTs2]      #vettore con le due posizioni dei minimi 
    rDiag = [rEce,rTs]            # Costriusco un vettore con i due vettori delle posizioni
    pos = [posE,posT]             
    # RHo - UP : calcolo d2 cm in meno e in più
    # rispetto alla poszione del massimo dei due minimi
    # e vedo il punto più vicino della medesima curva
    # Come RHO-UP prendo poi il massimo dei due valori trovati + 0.001
    
    diagM = diag[ind]  # La diagnostica che ha il massimo tra i due minimi
    radii = rDiag[ind]     # Le posiszioni della diagnostica di cui sopra
    
    pos1 = pos[ind]- d2  # Posizione a distanza d/2 dal centro plasma
    pos2 = pos[ind] + d2
    
    idt1 = np.nanargmin(abs(radii-pos1)) # indice1 più vicino sulla curva più alta
    idt2 = np.nanargmin(abs(radii-pos2))
    # vettore dei valori più vicino alla posizione estrema dell'intervallo 
    #(per la diagn che ha il max tra i due min)
    temp = [iEce,i]
    indice1 = temp[ind]
    pippo = [v1,v2] = [diagM[idt1,indice1],diagM[idt2,indice1]] 
     
    siaM = np.max(pippo) # Seleziono il pun to più alto tra i due, da scegliere come rho-up
    rhoup = siaM + 0.0001 # aggiungo un piccolo delta di sicurezza
    
    #############################
    # RHO down
    # Trovo il punto più vicino della diag non considerata prima
    # seleziono l'altra diagnostica: minimo rta i due massimi
    indx = np.argmin(([mEce,mTs]))   # Trovo il minimo tra i due minimi
    rm = rDiag[indx]     # Le posiszioni della diagnostica di cui sopra
    diagm = diag[indx]      # seleziono la diagnostica 'piu bassa'
    # trovo l'indice del punto con rho piu vicina al minimo dell'altra:
    # avendo assi temporali diversi devo secegliere indici diversi a secondo del caso    
    temp = [iEce,i]
    indice2 = temp[indx]
    id_pluto = np.nanargmin(abs(diagm[:,indice2]-M)) 
    pluto = diagm[id_pluto,i]
    
    rm = rDiag[indx]     # posizioni radiali della diagostica 'minore'
    pos_pluto = rm[id_pluto]
    
    # se il punto è prima del minimo aggiungo due puntio per determinare la rho down
    # se invece è dopo, ne tolgo due
    # in ambedue i casi tolgo poi un piccolo margine
    if id_pluto < indT:
        rhodown = diagm[id_pluto+2,indice2] - 0.001
    else:
        rhodown = diagm[id_pluto-2,indice2] - 0.001
    ranges[i,:] = [rhodown,rhoup]       
    
    

    mask = (rhoTs2[:,i] >= rhodown) & (rhoTs2[:,i] <= rhoup)    # Seleziono gli indici dei tempi di interesse 
    temp = temp_ts2[mask,i]                    # Valori di temperaatura nelle rho selezionate
    sig = err_ts2[mask,i]       # Errore sui singoli valori di temperatura: sigma
    ai = 1/sig**2            # Inverso del quadrato delle sigma
    somma = sum(ai)
    if somma==0:
        somma=1
    xm_ = sum(ai*temp)/somma   # Media dei valori di temperatura pesata sui singoli errori
    temp_tsM_rho[i] = (xm_)               # Valore delle temperature 'attese'/più probabili, nel tempo
    err_xm_ = np.sqrt(1/somma) 
    err_tsM_rho[i] = (err_xm_)

    mask1 = (rhoEce22[:,i] >= rhodown) & (rhoEce22[:,i] <= rhoup)    # Seleziono gli indici delle psi di interesse 
    temp_eceM_rho[i]  = np.mean(temp_ece22[mask1,i], axis=0)
    rho_eceM_rho[i] = np.mean(rhoEce22[mask1,i],axis=0)
    err_eceM_rho[i] = np.mean(err_ece22[mask1,i], axis=0)

   
print('Dimensioni temp_eceM_rho = ',temp_eceM_rho.shape)
print('Dimensioni temp_TsM_rho = ',temp_tsM_rho.shape)    

# Si trova la matrice 'ranges' con tutti i valoori di rhoup/down a tutti gli istanti del TS
# Si trovano qiundi le due matrici: temp_eceM_rho e temp_tsM_rho 
# che contengono i valori di Te a tutte le rho e a tutte gli istanti tra t1 e t2
# assi temporali: timeTs2 e timeEce2
# le posizioni: rTs e rEce
 
fig, ax = plt.subplots(3)
ax[0].plot(timeTs2,timeEce22, label='cfr tempi')
ax[0].legend()
ax[0].set_title('Check Plot')
ax[1].plot(timeTs2,ranges[:,1], label='rho-up')
ax[1].plot(timeTs2,ranges[:,0], label='rho-down')
ax[1].legend()
ax[2].plot(timeTs2, temp_tsM_rho, label='TS averaged')
ax[2].plot(timeEce22,temp_eceM_rho, label='ECE averaged')
ax[2].legend()

left = min(np.nanmin(temp_eceM_rho), np.nanmin(temp_tsM_rho))
right  = max(np.nanmax(temp_eceM_rho), np.nanmax(temp_tsM_rho)) 
def retta(x):
    return x

x = np.linspace(left-0.5,right+0.5,timeTs2.shape[0])   # Range in keV di dove tracciare la retta
y = retta(x)  

fig,ax=plt.subplots(1)
# ax.scatter(temp_eceM_rho, temp_tsM_rho,label='Rho Averaged ECE vs TS ') 
ax.plot(x,y,'g--', lw=.8)
ax.errorbar(temp_tsM_rho, temp_eceM_rho, xerr = err_tsM_rho, yerr = err_eceM_rho, 
              marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5, label='ECE vs TS')
ax.set_title(f'JPN {shot} - Te Ece-Michelson vs Te HRTS  for {tlim1:.2f}<t<{tlim2:.2f} (s)') #:.2f per avere 2 cifre decimali
ax.set_xlabel('Te HRTS (keV)')
ax.set_ylabel('Te Ece-Michelson (keV)')
ax.legend()
 
 
 
 
 