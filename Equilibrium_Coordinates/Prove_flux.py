#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 08:55:40 2024
@author: lsenni
Prove ricostruzione dati su coordinata equilibrio psi
"""

import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt

plt.close('all')  
shot = 96994
tlim = 47  # Seleziono l'istante al quale voglio fare il profilo - in secondi

w = ppfs(shot) 

HRTS = w.hrts
Tts = HRTS.te   
psi_ts = HRTS.psi
efit = w.efit
psi_rz = efit.psi  # psi(R,Z)
psin = efit.psni
psir = efit.psir
psiz = efit.psiz
RTS = HRTS.yr # coord R della linea di vista hrts 
ZTS = HRTS.z  # coord Z della linea di vista hrts

ind = np.argmin(abs(Tts.t - tlim)) # indice corrispondente al tempo selezionato
sel_time= "%.2f" %Tts.t[ind] # selected time arrotondato alla seconda cifra decimale

# Faccio la slice temporale
pr_tts = Tts.v[:,ind]  # porfilo di Te thomson
psi_tts = psi_ts.v[:,ind] # coordinata psi allo stesso istante temporale
psi_gen = psin.v[:,ind] # Valori di psi normalizzata all'istante selezionato
psi = psi_rz.v[:,ind] # Valori di psi normalizzata all'istante selezionato

plt.figure(f'Profilo di psi a t = {sel_time}')
plt.plot(psi_ts.r,psi_tts)

plt.figure()
plt.plot(psi_rz.r,psi)

plt.figure()
plt.plot(psir.v,psiz.v)

plt.figure()
plt.plot(RTS.v,ZTS.v)
plt.title('linea di vista HRTS nel piano poloidale')

# fig,(ax0,ax1) = plt.subplots(nrows=2)
# ax0.plot(psi_tts,pr_tts)
# ax0.set_title(f'Profilo di temperatura in psi, t = {sel_time}')
# ax0.set_ylabel('Te TS')  

# ax1.plot(Tts.r,pr_tts)
# ax1.set_title('Profilo di temperatura in R')
# ax1.set_ylabel('Te TS')










