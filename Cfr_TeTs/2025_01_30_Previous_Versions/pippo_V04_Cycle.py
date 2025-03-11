#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 11/10/2023 - lsenni
Confornto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

V04_Cycle: Ciclo per  Plot di sovrapposizione tra diversi spari 
           degli andamenti Te-HRTS vs Te-ECE
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt
from scipy.signal import resample

plt.close('all')  

# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
#shot = 104522
tlim1 = 47  # limite inferiore selezione tempi
tlim2 = 52 
rad = 3

shots = [104520,104521,104522,104523,104524,104525,104526]


def retta(xx):
    return xx
xx = np.linspace(-1,13)
yy = retta(xx)

fig, ax= plt.subplots(figsize = (8,6))
ax.plot(xx,yy,'g--')

for shot in shots:
        
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
    
    
    dim = Tts.t.size
    v,t = resample(Tece.v, dim, t=Tece.t, axis=1)
    Tece.v = v
    Tece.t = t
    
    x = Tece.v[0,:]/1000
    y = Tts.v[0,:]/1000
    ax.scatter(x,y,s=5,marker='o',label=f'{shot}')
    
ax.set_xlim(left=2)
ax.set_ylim(bottom=2)
ax.legend()
ax.set_title('Te HRTS vs Te Ece-Michelson')
ax.set_xlabel('Te Ece-Michelson (keV)')
ax.set_ylabel('Te HRTS (keV)')

plt.savefig('Tts_vs_Tece_MultiShot',dpi=600)   

    
#plt.savefig(f'{shot}_Tts_vs_Tece')
    #plt.close()
    
    