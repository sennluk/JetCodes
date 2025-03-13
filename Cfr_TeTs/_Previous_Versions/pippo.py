#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 21:45:26 2023
author: egio
Modified 09/10/2023 - lsenni
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt


shot = 104522

w = ppfs(shot)

# GET: andamento nel tempo ad un certoi raggio,
# Fornisce una tupla con due volonne, la prima i tempi la seconda i valori
rad=0
Tet = w.hrts.te.get(r=rad)

plt.figure('1', clear=True)
plt.plot(Tet[0],Tet[1],label='Time trend at r=')
plt.legend()
plt.title(shot)
plt.xlabel('time(s)')
plt.ylabel('Electron temperature (keV)')


# SLICE: Andamento su una slice selezionata nel tempo o nel raggio
# fornisce ancora  .v - .r - .t
Ter = w.hrts.te.slice(t=50)

plt.figure('2', clear=True)
plt.plot(Ter.r,Ter.v,label='Profile')
plt.legend()
plt.title(shot)
plt.xlabel('Radii(m)')
plt.ylabel('Electron temperature (keV)') 

Tts = w.hrts.te
rad = Tts.r
idx = rad>5

tlim1 = 46  # limite inferiore selezione tempi
tlim2 = 54 

par = Tts
def trange(par, tlim1,tlim2):
            idt = (par.t >= tlim1) & (par.t <= tlim2)
            par.t = par.t[idt]
            par.v = par.v[:,idt]
return par
Tts = par    
    

dim = Tts.t.size
Tec = w.ecm1.prfl.tresample(n=dim)
t = Tts.t
Tets = w.hrts.te.get(r=0)
Teec = w.ecm1.prfl.get(r=0)



# def Ratio(Tts, Tec):
#     return Tts/Tec
#Ratio = Tets[1]/Tec.v[30,:]

plt.figure('3', clear=True)
plt.plot(t,Tets[1],label='Tets')
plt.plot(t,Tts.v[],label='Teec')
plt.legend()
plt.title(shot)
plt.xlabel('time(s)')
plt.ylabel('Ratio') 