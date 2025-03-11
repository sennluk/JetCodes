#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:54:14 2024

@author: lsenni
01: provo a ricostruire le psi nel tempo
V01_en: uso l'append per costruire l'array
"""


import my_flush
import numpy as np
import ppfeg
import matplotlib.pyplot as plt


shot = 104522
w = ppfeg.ppfs(shot)       # richiamo i dati dello shot indicato

time_ax = w.ecm1.prfl.t
prfl = w.ecm1.prfl  # Canale del profilo
zece = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y

psi = []
for i,time in enumerate(time_ax):
    r, te = prfl.get(t=time)  # Profilo a T=time
    z = np.full_like(r, zece) # Retta linea di vista parallela asse r
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi_, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
    psi.append(psi_)

# plt.figure() # psi(r) al tempo t=time
# plt.plot(r,psi)
