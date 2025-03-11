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

# Bt = 3.85
# Ti = 6000
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
    return 0.42*10**-2*(np.sqrt(Ti1))/(a*Bp)

def nustar(q, R, A, Zeff, ne, Te):    # Collisionality
    return 0.8*10**-2*q*R*A**(3/2)*Zeff*ne/Te**2

def betan(rhos, nus):
    return 550*rhos**1.02 * nus**0.19

# I valori di Ti, Te, ne sono presi come  media di volume
Bt = abs(w.magn.bvac.v) # Toroidal magnetic field on axis (T)
Ip = abs(w.magn.ipla.v)/10**6 # Plasma current at r=0 (MA)
Bp = abs(w.efit.bpme.get(r=0)[1]) # Measured B poloidal at r=0  (T)
tq = w.magn.bvac.t
Zeff = w.zeff.zefv.v[0,:] # Zeff vertical
tZeff = w.zeff.zefv.t
ne =  w.hrtx.nevl.v/10**20 # Volume averaged electron density /10**20 m^-3
t = w.hrts.ne.t
Te = (w.hrtx.tevl.v/10**3)+0.001   # Volume averaged electron temperature in keV
Ti1 = w.xcs.ti.get(r=0)[1]/4 # Ion temperature at r=0 in trange - XCS channel
tTi1 = w.xcs.ti.t
# /4 per approssimare la media di volume

# t,ne = w.hrts.ne.trange(tlim).get(r=0) # electron dens at r=0 in trange
# tq, _ = w.efit.q95.get(r=0)
_, Bvac = w.efit.bvac.get(r=0)
Btn = w.efit.btnm

# q = np.interp(t, tq, q95)
# Bt = np.interp(t, tq, Bvac)
#q = 3

q = qu(R, Bt, A, Ip)
q = np.interp(t, tq, q[0,:])
Ti1 = np.interp(t, tTi1, Ti1)
Zeff = np.interp(t, tZeff, Zeff)
Bt = np.interp(t, tq, Bt[0,:])

Bp = Bpol(Bt, A, q)
rhosi = rhostar(Meff, Ti1, Bp, a)
nuse = nustar(q, R, A, Zeff, ne, Te)
betn = betan(rhosi, nuse)


fig, ax = plt.subplots(4,1,sharex=True, num='Ciao', clear=True)
fig.subplots_adjust(hspace=0)
ax[0].plot(t,Te[0,:], label='Te')
ax[0].legend()
ax[1].plot(t, betn[0,:], label='Calc betn')
ax[1].plot(Btn.t,Btn.v[0,:], label ='EFIT Bn')
ax[1].legend()
ax[2].plot(t, rhosi, label = 'RhoS_i')
ax[2].legend()
ax[3].plot(t, nuse[0,:], label = 'NuS_e')
ax[3].legend()

# print('t =', ts, 'sec')
# print('Measured Bt =', Bt, 'T')
# print('Measured Ip =', Ip, 'MA')
# print('Measured Bp =', Bp, 'T')
# print('calculated Bpol (Bt/(A*q)) = ', Bpol, 'T')




