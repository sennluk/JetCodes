#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 11/10/2023 - lsenni
Per prove importa dati e profili
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt
from scipy.signal import resample
import pandas as pd

plt.close('all')  

# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
# FAccio l'analisi sui tempi scalti, ma visualizzo + e - due secondio
shot = 104522
tlim1 = 49  # limite inferiore selezione tempi 
tlim2 = 50 
rad = 3

w = ppfs(shot)

Tts = w.hrts.te       # canale temp elettronica HRTS
ne = w.hrts.NE        # canale temp elettronica HRTS
psi = w.hrts.psi      # Canale Normalized poloidal flux (EFIT) hrts
Tece = w.ecm1.prfl    # Can El Temp ECE Michelson

t0 = 50
idx = np.argmin(abs(Tts.t - t0))

Prof = Tts.v[:,idx]

plt.figure()
plt.plot(Tts.r,Prof)

plt.figure()
plt.plot(psi.v[:,idx],Prof)

prov = [psi.v[:,idx],Tts.v[:,idx]]
df = pd.DataFrame(prov)
df.to_csv("prova.csv")

# carlo = Tts.v
# np.savetxt('carlo.txt',())

# df = pd.DataFrame(Tts.v)
# df.to_csv("Tts.csv")
# df.to_csv(r'c:\data.pandas.txt', header=None, index=None, sep='\t', mode='a')
# np.savetxt(r'c:\data.pandas.txt',df.values, fmt='%d', delimiter='t')

# Riporto i dati del Thomson all'interno dell'intervallo temporale
# e sul raggio scelto e gli errori --> commenti nelle versioni precedenti
idt = (Tts.t >= tlim1) & (Tts.t <= tlim2)
TimeTS = Tts.t[idt]
TempTS = Tts.v[:,idt]

idr_psi = np.argmin(abs(psi.r - rad))
psiTS = psi.r[idr_psi][np.newaxis]
TempTS = Tts.v[idr_psi,:][np.newaxis,:]



# Plot 
plt.figure('1', clear=True)
plt.plot(psiTS,TempTS[:,500],label='Tts')
plt.legend()
plt.title(f'Shot n.{shot} - Te HRTS & Ece-Michelson')
plt.xlabel('psi')
plt.ylabel('Te(keV)') 



##############################################################################
##############################################################################




# Riporto i dati del Ece Michelson all'interno dell'intervallo temporale
# e sul raggio scelto
idte = (Tece.t >= tlim1) & (Tece.t <= tlim2)
TimeECE = Tece.t[idte]
TempECE = Tece.v[:,idte]
idre = np.argmin(abs(Tece.r - rad))
RadECE = Tece.r[idre][np.newaxis]
TempECE = Tece.v[idre,:][np.newaxis,:]

dim = TimeTS.size
v,t = resample(TempECE, dim, t=TimeTS, axis=1)
TempECE = v
TimeECE = t


Time = Tece.t
index = (Time >= tlim1) & (Time <= tlim2) # Definisco gli indici sui quali vengono poi
# fatte le analisi
# p = Tece.v[0,index]/1000
# q = Tts.v[0,index]/1000
# ErrY = ErrY [idx]

# Plot degli andamenti nel tempo delle temperature misurate con Ts e ECE
plt.figure('1', clear=True)
plt.plot(Tts.t,Tts.v[0,:],label='Tts')
plt.plot(Tece.t,Tece.v[0,:],label='Tece')  
plt.legend()
plt.title(f'Shot n.{shot} - Te HRTS & Ece-Michelson')
plt.xlabel('time(s)')
plt.ylabel('Te(keV)')  

# Plot della Te TS vs Te ECE + retta x=y
# plt.figure(f'{shot}_Tts_vs_Tece', clear=True)
# plt.scatter(Tece.v[0,:]/1000,Tts.v[0,:]/1000,s=5,marker='o',label='ECE vs TS')
# plt.xlim(left=2)
# plt.ylim(bottom=2)
# plt.legend()
# plt.title(f'Shot n.{shot} - Te HRTS vs Te Ece-Michelson')
# plt.xlabel('Te Ece-Michelson (keV)')
# plt.ylabel('Te HRTS (keV)')
# plt.savefig(f'{shot}_Tts_vs_Tece')
