#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:16:33 2024
@author: lsenni
Fig 02 per poster Alto Beta Orsitto ad EPS 2024 
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import matplotlib

plt.close('all') 
 
shot = 103117   # 99950, 99971
tlim1 =  43     # limite inferiore selezione tempi - in secondi
tlim2 = 51.4 
delta = 3       # tempo in più di tlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.00      # 

w = ppfs(shot) 

pnbi = w.nbi.ptot   # Total NBI Power
idt = (pnbi.t >= tlim1) & (pnbi.t <= tlim2) 
tPnbi = pnbi.t[idt] - 40
pnbi = pnbi.v[0,idt]/10**6  # In MWatt   

pbolo = w.bolo.topi     # Tot rad power (improved) - Bolometry
idt = (pbolo.t >= tlim1) & (pbolo.t <= tlim2)
tPrad = pbolo.t[idt] - 40
prad = pbolo.v[0,idt]/10**6

betan = w.efit.btnm      # Beta normalizaed - MHD
idt = (betan.t >= tlim1) & (betan.t <= tlim2)
tBetan = betan.t[idt] - 40
betan = betan.v[0,idt]

wdia = w.efit.wdia   # Energy meas - diamag
idt = (wdia.t >= tlim1) & (wdia.t <= tlim2)
tWdia = wdia.t[idt] - 40
wdia = wdia.v[0,idt]/10**6  # MegaJoule

zeffv = w.zeff.zefv # Zeff vertical
idt = (zeffv.t >= tlim1) & (zeffv.t <= tlim2)
tZeffv = zeffv.t[idt] - 40
zeffv = zeffv.v[0,idt]

zeffh = w.zeff.zefh # Zeff Horizontal
idt = (zeffh.t >= tlim1) & (zeffh.t <= tlim2)
tZeffh = zeffh.t[idt] - 40
zeffh = zeffh.v[0,idt]

tmax = w.hrtx.tmax      # max Te max hrts (hrtx channel)
idt = (tmax.t >= tlim1) & (tmax.t <= tlim2)
tTmax = tmax.t[idt] -40
tmax = tmax.v[0,idt]/1000  # in keV

ti = w.xcs.ti  # Ion temperature 
idt = (ti.t >= tlim1) & (ti.t <= tlim2)
tTi = ti.t[idt] - 40
ti = ti.v[0,idt]/1000  # in keV

# ti = w.ks5

dens = w.hrts.ne.slice(r=rad)  # density profile hrts at r=rad
idt = (dens.t >= tlim1) & (dens.t <= tlim2)
tNe = dens.t[idt] -40
ne = dens.v[0,idt]/10**19

# qu = w.efit.qax       # Simulated q on axis
# idt = (qu.t >= tlim1) & (qu.t <= tlim2)
# tQu = qu.t[idt]
# qu = qu.v[0,idt]

nr = w.tin.rnt # KN1 neutron rate
idt = (nr.t >= tlim1) & (nr.t <= tlim2)
tNr = nr.t[idt]-40   
nr = nr.v[0,idt]/10**16

# gas = w.gash.eler # KN1 neutron rate
# idt = (gas.t >= tlim1) & (gas.t <= tlim2)
# tGas = gas.t[idt]
# gas = gas.v[0,idt]/10**21


###############################
linew = 0.5  # plot  lines dimension
fonts = 5.5    # size of the legend fonts
fs = 6       # size of the axes labels
fst = 6     # size of the ticks

# fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(nrows=9, sharex=True, num=f'JPN {shot} Trends over time')
fig, (ax1, ax2, ax3, ax5, ax6, ax8) = plt.subplots(nrows=6, sharex=True, num=f'JPN {shot} Trends over time')
plt.subplots_adjust(hspace=0)

ax1.plot(tPnbi,pnbi, lw = linew, label='$P_{NBI}$')
ax1.plot(tPrad, prad, lw = linew, label='$P_{rad}$')
ax1.set_ylabel('MW', fontsize= fs)
ax1.yaxis.set_tick_params(labelsize=fst)
ax1.set_ylim(None,25)
ax1.legend(fontsize = fonts, loc="upper right")
ax1.set_title(f'Traces JPN {shot} - $B_T$ = 2.4 T, I = 1.4 A')

ax2.plot(tBetan, betan, lw = linew, color='b', label= r'$\beta_N$')
# ax2.set_ylabel('a.u.', fontsize= fs)
ax2.yaxis.set_tick_params(labelsize=fst)
ax2.legend(fontsize = fonts, loc="upper right")

ax3.plot(tWdia, wdia, lw = linew, color='olive', label='$W_{dia}$')
ax3.set_ylabel('MJ', fontsize= fs)
ax3.yaxis.set_tick_params(labelsize=fst)
ax3.legend(fontsize = fonts, loc="upper right")

# ax4.plot(tZeffv, zeffv, lw = linew, label='$Z_{eff}$ vertical')
# ax4.plot(tZeffh, zeffh, lw = linew, label='$Z_{eff}$ horizontal')
# # ax4.set_ylabel('a.u.', fontsize= fs)
# ax4.yaxis.set_tick_params(labelsize=fst)
# ax4.legend(fontsize = fonts, loc="upper right")

ax5.plot(tTmax, tmax, lw = linew, color='r', label='$T_{e}$ max')
ax5.plot(tTi, ti, lw = linew, color='g', label='$T_{i}$') #  <$N_i$>26
ax5.set_ylabel('keV', fontsize= fs)
ax5.yaxis.set_tick_params(labelsize=fst)
ax5.legend(fontsize = fonts, loc="upper right")

ax6.plot(tNe, ne, lw = linew, color='darkblue', label=f'Density at R={rad} m')
ax6.set_ylabel('$10^{19} m^{-3}$', fontsize= fs)
ax6.yaxis.set_tick_params(labelsize=fst)
ax6.legend(fontsize = fonts, loc="upper right")

# ax7.plot(tQu, qu, lw = linew, color = 'orange',label='Simu q on axis')
# # ax7.set_ylabel('a.u.', fontsize= fs)
# ax7.yaxis.set_tick_params(labelsize=fst)
# ax7.legend(fontsize = fonts, loc="upper right")
# ax7.set_xlabel('time(s)', fontsize= fs)

ax8.plot(tNr, nr, lw = linew, color='darkblue', label='Neutron rate')
ax8.set_ylabel('$10^{16} n/s$', fontsize= fs)
ax8.yaxis.set_tick_params(labelsize=fst)
ax8.legend(fontsize = fonts, loc="upper right")
ax8.set_xlabel('time (s)')

# ax9.plot(tGas, gas, lw = linew, color='darkblue', label='Total electron flow rate')
# ax9.set_ylabel('$10^{21} Trons/s$', fontsize= fs)
# ax9.yaxis.set_tick_params(labelsize=fst)
# ax9.legend(fontsize = fonts, loc="upper right")
# ax9.xaxis.set_tick_params(labelsize=fst)

plt.savefig(f'Fig_01_JPN{shot}_Multiplot.pdf',dpi=300)

