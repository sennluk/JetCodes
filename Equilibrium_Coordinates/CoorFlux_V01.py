#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:54:14 2024

@author: lsenni
Ricostruzione dei profili di psi nel tempo
Implemento solo per KK1 

"""
import my_flush
import numpy as np
import ppfeg
import matplotlib.pyplot as plt

shot = 104522

# def mycoord(shot):
w = ppfeg.ppfs(shot)       # richiamo i dati dello shot indicato

time_ax = w.ecm1.prfl.t
prfl = w.ecm1.prfl  # Canale del profilo
zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y
psiKk1 = np.zeros(prfl.v.shape)

for i,time in enumerate(time_ax):
    r, te = prfl.get(t=time)  # Profilo a T=time
    z = np.full_like(r, zKk1) # Retta linea di vista parallela asse r
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
    psiKk1[:,i] = psi
        
    # return psiKk1   

# pippo = mycoord(104520)     

# plt.figure() # psi(r) al tempo t=time
# plt.plot(r,psi)

# obj = enumerate(time_ax)
# print(obj)
# print(list(obj))