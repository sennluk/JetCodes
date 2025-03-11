#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 21:45:26 2023
@author: egio/ lsenni
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt

shot = 103117
perc = 0.5 # percentuale Deuterio 
R = 3.0
a = 1.0
tlim = [40.9, 53.2]
ts = 46 # time selected  

w = ppfs(shot)

Bt = 3.85
Ti = 6000
Mid = 2 # Massa ione deuterio
Mit = 3 # massa ione Trizio
Mp = 1 # D massa protoni

Meff = ((Mid*perc)+(Mit*(1-perc)))/(perc+(1-perc)) # Massa efficace
A = R/a   # Aspect ratio
eps = a/R
print('Massa efficace=',Meff, 'Aspect ratio=',A)

def qu(R, Bt, A, Ip):
    return 5*R*Bt/(A**2*Ip)  

def Bpol(Bt, A, q):
    return Bt/(A*q)         # Poloidal magnetic field

def rhostar(Meff, Ti, Bp, a):   # Raggio di larmor ioni per beta poloidale, normalizzato
    return 0.42*10**-2*np.sqrt(Ti)/(a*Bpol)

def nustar(q, R, A, Zeff, ne, Te):    # Collisionality
    return 0.8*10**-2*q*R*A**(3/2)*Zeff*ne/Te**2

def betan(rhos, nus):
    return 550*rhos**1.02 * nus**0.19

Bt = abs(w.magn.bvac.get(t=ts)[1]) # Toroidal magnetic field at t=ts (T)
Ip = abs(w.magn.ipla.get(t=ts)[1]/10**6) # Plasma current at t=ts (MA)
Bp = abs(w.efit.bpme.get(r=0,t=ts)) # Measured B poloidal at t=ts  (T)
Zeff = w.zeff.zefv.get(t=ts) # Zeff vertical
ne =  w.hrtx.nevl.get(r=0,t=ts)/10**20 # Volume averaged electron density /10**20 m^-3
Te = w.hrtx.tevl.get(r=0,t=ts)/10**3   # Volume averaged electron temperature in keV
# I valori di Ti, Te, ne sono intesi come  media di volume

_,Ti = w.hrts.te.trange(tlim).get(r=0)
t,ne = w.hrts.ne.trange(tlim).get(r=0)
tq, q95 = w.efit.q95.get(r=0)
_, Bvac = w.efit.bvac.get(r=0)
q = np.interp(t, tq, q95)
Bt = -np.interp(t, tq, Bvac)
#q = 3

q = qu(R, Bt, A, Ip)
Bp = Bpol(Bt, A, q)
rhosi = rhostar(Meff, Ti, Bp, a)
nuse = nustar(q, R, A, Zeff, ne, Te)
betn = betan(rhosi, nuse)

print('t =', ts, 'sec')
print('Measured Bt =', Bt, 'T')
print('Measured Ip =', Ip, 'MA')
print('Measured Bp =', Bp, 'T')
print('calculated Bpol (Bt/(A*q)) = ', Bpol, 'T')




fig, ax = plt.subplots(3,1,sharex=True, num='Ciao', clear=True)
fig.subplots_adjust(hspace=0)
ax[0].plot(t, betn)
w.efit.btnm.plot(ax=ax[0])
# ax[0].plot()
ax[1].plot(t, rhosi)
ax[2].plot(t, nuse)


