#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:18:00 2024
@author: lsenni
Prove per usare le routines di Flush
OSS: le coordinate vanno messe in cm!
"""

import my_flush
import numpy as np
import ppfeg
import matplotlib.pyplot as plt


shot = 104522
w = ppfeg.ppfs(shot)       # richiamo i dati dello shot indicato

zece = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y

prfl = w.ecm1.prfl # canale del profilo
time = 48.0       # Tempo scelto
r, te = prfl.get(t=time) 

#r = np.linspace(2.5, 3.8, 1000)
z = np.full_like(r, zece)

ts, ier = my_flush.flushinit(15, shot, time)
psi, ier = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanno messe in cm!

plt.figure()
plt.plot(psi,te)

plt.figure()
plt.plot(r,te)

plt.figure()
plt.plot(r,psi)