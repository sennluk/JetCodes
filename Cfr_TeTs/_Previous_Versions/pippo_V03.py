#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 11/10/2023 - lsenni
Confornto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

Figura di confronto riportando ogni diagnostica con il proprio rate di acquisizione
e 
Figura con il rapporto facendo il reshape dell'asse temporale dell'Ece, che e' piu veloce,
su quello dell'hrts che e' piu' lento'
V03:Versione aggiustata dei plot di ~Te ECE e Te Ts vs time e uno verso l'altro
Aggiungo la divisione 1000 per mettere in keV e qualche altro particolare'
    
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt
from scipy.signal import resample

plt.close('all')  

# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
shot = 104522
tlim1 = 47  # limite inferiore selezione tempi
tlim2 = 52 
rad = 3

w = ppfs(shot)

Tts = w.hrts.te       # canale temp elettronica HRTS
Tece = w.ecm1.prfl    # Can El Temp ECE Michelson

# Riporto i dati del Thomson all'interno dell'intervallo temporale
# e sul raggio scelto e gli errori
idt = (Tts.t >= tlim1) & (Tts.t <= tlim2)
Tts.t = Tts.t[idt]
Tts.v = Tts.v[:,idt]
idr = np.argmin(abs(Tts.r - rad))
Tts.r = Tts.r[idr][np.newaxis]
Tts.v = Tts.v[idr,:][np.newaxis,:]

# Riporto i dati del Ece Michelson all'interno dell'intervallo temporale
# e sul raggio scelto
idte = (Tece.t >= tlim1) & (Tece.t <= tlim2)
Tece.t = Tece.t[idte]
Tece.v = Tece.v[:,idte]
idre = np.argmin(abs(Tece.r - rad))
Tece.r = Tece.r[idre][np.newaxis]
Tece.v = Tece.v[idre,:][np.newaxis,:]

# Plot degli andamenti nel tempo delle temperature misurate con Ts e ECE
plt.figure('1', clear=True)
plt.plot(Tts.t,Tts.v[0,:]/1000,label='Tts')
plt.plot(Tece.t,Tece.v[0,:]/1000,label='Tece')  
plt.legend()
plt.title(f'Shot n.{shot} - Te HRTS & Ece-Michelson')
plt.xlabel('time(s)')
plt.ylabel('Te(keV)')  

# Faccio il resampling dei dati Ece sulla base
# dell'asse temporale del HRTS

dim = Tts.t.size
v,t = resample(Tece.v, dim, t=Tece.t, axis=1) # Ricampiono Tece.v al numero di campioni dato da Dim
Tece.v = v
Tece.t = t

Ratio = Tts.v[0,:]/Tece.v[0,:]
re = np.linspace(1,1,dim)

# Plot dell-andamento del rapporto nel tempo + retta a Rat=1
# plt.figure('2', clear=True)
# plt.plot(Tts.t,Ratio,label='Ratio')
# plt.plot(Tts.t,re,'--')
# plt.legend()
# plt.title(shot)
# plt.xlabel('time(s)')
# plt.ylabel('Ratio Te_th / Te_ece')

def retta(x):
    return x
x = np.linspace(-1,13,dim)
y = retta(x)

# Plot della Te TS vs Te ECE + retta x=y
plt.figure(f'{shot}_Tts_vs_Tece', clear=True)
plt.scatter(Tece.v[0,:]/1000,Tts.v[0,:]/1000,s=5,marker='o',label='ECE vs TS')
plt.plot(x,y,'g--')
plt.xlim(left=2)
plt.ylim(bottom=2)
plt.legend()
plt.title(f'Shot n.{shot} - Te HRTS vs Te Ece-Michelson')
plt.xlabel('Te Ece-Michelson (keV)')
plt.ylabel('Te HRTS (keV)')
plt.savefig(f'{shot}_Tts_vs_Tece')

# plt.figure(f'{shot} Te vs time') (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
fig,(ax0,ax1) = plt.subplots(nrows=2,sharex=True)
ax0.plot(Tts.t,Tts.v[0,:]/1000,label='Tts')
ax0.plot(Tece.t,Tece.v[0,:]/1000,label='Tece')
ax0.legend()
ax0.set_title(f'Shot n.{shot} - Te HRTS & Ece-Michelson')
ax0.set_ylabel('Te (keV)')  
ax1.plot(Tts.t,Ratio,label='Ratio')
ax1.plot(Tts.t,re,'--')
ax1.legend()
ax1.set_title('Te-HRTS/ Te-ECE')
ax1.set_xlabel('time(s)')
ax1.set_ylabel('Ratio Te-hrts / Te-ece')
plt.savefig(f'{shot}_Te_vs_Time')



