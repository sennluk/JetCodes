#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 11/10/2023 - lsenni
Confornto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 
    BR01:Elimino alcuni plot e
    Provo a fare le differenze
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


# Faccio il resampling dei dati Ece sulla base
# dell'asse temporale del HRTS

dim = Tts.t.size
v,t = resample(Tece.v, dim, t=Tece.t, axis=1)
Tece.v = v
Tece.t = t

Ratio = Tts.v[0,:]/Tece.v[0,:]
re = np.linspace(1,1,dim)
re2 = np.linspace(0,0,dim)

Diff = (Tts.v[0,:]-Tece.v[0,:])/1000

# Plot dell-anadmaneto della differenza nel tempo + retta a diff=0
plt.figure()
plt.plot(Tts.t,Diff,label='difference')
plt.plot(Tts.t,re2,'--')
plt.legend()
plt.title(shot)
plt.xlabel('time(s)')
plt.ylabel('Te_TS - Te_ECE')

def retta(x):
    return x
x = np.linspace(-1,13,dim)
y = retta(x)

# Plot della Te TS vs Te ECE + retta x=y
plt.figure(f'{shot}_Tts_vs_Tece', clear=True)
plt.scatter(Tece.v[0,:]/1000,Tts.v[0,:]/1000,marker='o',label='ECE vs TS')
plt.plot(x,y,'g--')
plt.xlim(left=2)
plt.ylim(bottom=2)
plt.legend()
plt.title(f'Shot n.{shot} - Te Ece-Michelson - Te HRTS')
plt.xlabel('Te Ece-Michelson (keV)')
plt.ylabel('Te HRTS (keV)')
#plt.savefig(f'{shot}_Tece-Tts')

# plt.figure(f'{shot} Te vs time') (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
fig,(ax0,ax1,ax2) = plt.subplots(nrows=3,sharex=True)
ax0.plot(Tts.t,Tts.v[0,:]/1000,label='Tts')
ax0.plot(Tece.t,Tece.v[0,:]/1000,label='Tece')
ax0.legend()
ax0.set_title(f'Shot n.{shot} - Te HRTS & Ece-Michelson')
ax0.set_ylabel('Te (keV)')  
ax1.plot(Tts.t,Ratio,label='Ratio')
ax1.plot(Tts.t,re,'--')
ax1.legend()
ax1.set_title('Te-HRTS/ Te-ECE')
#ax1.set_xlabel('time(s)')
ax1.set_ylabel('Ratio Te-hrts / Te-ece')
ax2.plot(Tts.t,Diff,label='difference')
ax2.plot(Tts.t,re2,'--')
ax2.legend()
ax2.set_title('Te_TS - Te_ECE')
ax2.set_xlabel('time(s)')
ax2.set_ylabel('Te_TS - Te_ECE')
#plt.savefig(f'{shot}_Te_vs_Time')




