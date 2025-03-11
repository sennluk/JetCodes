"""
Created on Tue Oct  3 21:45:26 2023
author: egio
created 11/10/2023 - lsenni
Confornto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

V02_Errors: Figura di confronto riportando ogni diagnostica con il proprio rate di acquisizione
            e 
            Figura con il rapporto facendo il reshape dell'asse temporale dell'Ece, che e' piu veloce,
            su quello dell'hrts che e' piu' lento'
Aggiungo le bare di errore e le perfeziono:
    
OSS/remember: se fuori dall'intervallo temporale giusto sono fuori scala
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt
from scipy.signal import resample

plt.close('all')  
# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
shot = 104522
tlim1 = 46  # limite inferiore selezione tempi
tlim2 = 54 
rad = 3

w = ppfs(shot)

Tts = w.hrts.te       # canale temp elettronica HRTS
Errts = w.hrts.dte  # Chann Errors HRTS
Tece = w.ecm1.prfl    # Can El Temp ECE Michelson

# Riporto i dati del Thomson e gli errori all'interno dell'intervallo temporale
# e sul raggio scelto e gli errori
idt = (Tts.t >= tlim1) & (Tts.t <= tlim2)
Tts.t = Tts.t[idt]
Tts.v = Tts.v[:,idt]
Errts.t = Errts.t[idt]
Errts.v = Errts.v[:,idt]
idr = np.argmin(abs(Tts.r - rad))
Tts.r = Tts.r[idr][np.newaxis]
Tts.v = Tts.v[idr,:][np.newaxis,:]
Errts.r = Errts.r[idr][np.newaxis]
Errts.v = Errts.v[idr,:][np.newaxis,:]

# Riporto i dati del Ece Michelson all'interno dell'intervallo temporale
# e sul raggio scelto
idte = (Tece.t >= tlim1) & (Tece.t <= tlim2)
Tece.t = Tece.t[idte]
Tece.v = Tece.v[:,idte]
idre = np.argmin(abs(Tece.r - rad))
Tece.r = Tece.r[idre][np.newaxis]
Tece.v = Tece.v[idre,:][np.newaxis,:]

ErrY = Errts.v[0,:]/2000 # /2 per la barra, /1000 perchè in keV
ErrY[ErrY>100] = 0

plt.figure('1', clear=True)
#plt.plot(Tts.t,Tts.v[0,:],yerr = ErrY,label='Tts')
plt.errorbar(Tts.t,Tts.v[0,:]/1000,yerr = ErrY,ecolor='g', elinewidth=.5,label='Tts')
plt.plot(Tece.t,Tece.v[0,:]/1000,label='Tece')  
plt.legend()
plt.title(shot)
plt.xlabel('time(s)')
plt.ylabel('Te Th e Ece')  

# plt.figure('2', clear=True)
# plt.errorbar(Tts.t,Tts.v[0,:],yerr = ErrY, label='Tts')
# plt.legend()
# plt.title(shot)
# plt.xlabel('time(s)')
# plt.ylabel('Te Th e Ece')  

# Faccio il resampling dei dati Ece sulla base
# dell'asse temporale del HRTS

dim = Tts.t.size
v,t = resample(Tece.v, dim, t=Tece.t, axis=1)
Tece.v = v
Tece.t = t

Ratio = Tts.v[0,:]/Tece.v[0,:]
re = np.linspace(1,1,dim)

# plt.figure('3', clear=True)
# plt.plot(Tts.t,Ratio,label='Ratio')
# plt.plot(Tts.t,re,'--')
# plt.legend()
# plt.title(shot)
# plt.xlabel('time(s)')
# plt.ylabel('Ratio Te_th / Te_ece')



