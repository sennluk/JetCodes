#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 16:45:03 2024

@author: egio
"""
from ppfeg import jetdata, ppfs
import matplotlib.pyplot as plt
import numpy as np
from mhd_elm_analysis import load_data, COLORS
from dataclasses import fields
import equil
from tqdm import tqdm

Z0 = 0.248

def names(pr):
    return [f"{fh.name}:{fh.type}" for fh in fields(pr)]

taglio = {0:3.792,1:3.81, 2:3.791, 3:4.0, 4:3.80, 8:3.80, 9:3.80} 

shot = 102071

pr = load_data(shot, 48.6, 49.2)
w = ppfs(shot)
mhk3 = w['egio':326].mhk3
ampl = mhk3.ampl
id_redge = (pr.tek.r > 3.5) & (pr.tek.r < 3.85)
te = pr.tek.v[id_redge,:]
rr = pr.tek.r[id_redge]
id_min = te.argmin(axis=0)
r_min = rr[id_min]
eql = equil.load_equilibrium(shot, 'EFTP')
psik = np.zeros((ampl.t.size, ampl.r.size))
for i,t in enumerate(tqdm(ampl.t)):
    psik[i,:] = eql.get_psi(t, ampl.r, Z0)[0]

fig, ax = plt.subplots(4,1, sharex=True, num='figura_1', clear=True)
fig.subplots_adjust(hspace=0)
fig.align_ylabels()
ax_te, ax_ne, ax_freq, ax_ampl = ax
(pr.tek/1e3).moveto(40).plot(r=3.7, ax=ax_te)
(pr.ne_kg/1e19).moveto(40).plot(r=3.7, ax=ax_ne)

#ax_ampl.pcolormesh(ampl.t - 40, psik.T, ampl.v, vmax=200, cmap='jet', rasterized=True)
ampl.moveto(40).pcolormesh(ax=ax_ampl, vmax=200, cmap='jet', rasterized=True)
ax_ampl.plot(pr.tek.t-40, r_min, color='k')

labels = set()
for elm in pr.elms:
    for n, mode in elm.modes.items():
        for t, f in zip(mode.times, mode.frequency):
            if n not in labels:
                label = f"N={n}"
                labels.add(n)
            else:
                label = None
                
            ax_freq.plot(t-40, f/1e3, color=COLORS[n], label=label)

ax_te.set(ylim=[1,2], ylabel='Te R = 3.7 m\n(keV)')
ax_ne.set(ylim=[0,5], ylabel='ne R = 3.7 m\n($10^{19}m^{-3}$)')
ax_freq.set(ylim=[0,None], ylabel = "frequency\n(kHz)")
ax_ampl.set(ylabel='R\n(m)')
ax_freq.legend(loc=2,labelspacing=0.2, fontsize='small', borderaxespad=0)
ax_ampl.set(ylim=[3.5, 3.85])
ax[-1].set(xlim=[pr.t_start-40, pr.t_end-40], xlabel='t (s)')

for axx, ch in zip(ax, 'abcd'):
    axx.text(1.0,1.0, ' '+ch+')', transform=axx.transAxes, ha='left', va='top')


for axx in ax:
    for t in pr.t_elm:
        axx.axvline(t -40, color='k', ls='--', lw=0.5)


#fig.savefig('abstract_fig.svg',dpi=200)
#fig.savefig('abstract_fig.png',dpi=200)

