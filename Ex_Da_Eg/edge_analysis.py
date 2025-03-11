#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 17:12:10 2023

@author: egio
"""
import numpy as np
import pandas as pd
from ppfeg import jetdata, V2d, ppfs
from modenew import ModeAnalysisNew
from dataclasses import dataclass
import matplotlib.pyplot as plt
from additional_color_tables import mplrainbow
from tqdm import tqdm

def my_figures(name, nr, nc=1, **kwargs):
    fig, ax = plt.subplots(nr,nc, sharex='col', num=name, clear=True, **kwargs)
    fig.subplots_adjust(hspace=0)
    return fig, ax

def randomise(t):
    dt = t[1] - t[0]
    dts = (np.random.rand(t.size) - 0.5)*dt*0.8
    return t + dts

def rand_plot(v, r=None, ax=None, **kwargs):
    t, vv = v.get(r=r)
    # add random number to t in roder to plot it better
    t = randomise(t)
    ax.plot(t, vv, **kwargs)
    
def delupper(*ax):
    for axx in ax:
        axx.yaxis.get_major_locator().set_params(prune='upper')

def checkne(nev, redge):
    t, ne = nev.get(r=redge)
    idt = np.argmin(np.abs(t-49.5))
    return ne[idt] < 0.54e20

def peak(ve, redge, r0):
    t, ve_edge = ve.get(r=redge)
    _, ve0 = ve.get(r=r0)
    t = randomise(t)
    ve_peak = ve0/ve_edge
    return t, ve_peak
        

REDGE = 3.75
df = pd.read_csv('pulse_data.csv')
df = df[df.plasma_current > 2.3]

fig, ax = my_figures('Ne', 5, figsize=(7,8))
fig.align_ylabels()
ax_ne, ax_te, ax_pe, ax_nepeak,ax_nbi = ax
for line in tqdm(df.itertuples(), total=len(df)):
    shot = line.pulse_number
    ipla = line.plasma_current
    if shot == 104472:
        continue
    ne = jetdata(shot, 'prts','ne', userid='egio')
    te = jetdata(shot, 'prts','te', userid='egio')
    pe = jetdata(shot, 'prts','pe', userid='egio')
    if not ne.ok():
        continue

    nbi = jetdata(shot, 'nbi', 'ptot')
        
    t, ne_peak = peak(ne, REDGE, 3.0)
    #t, pe_peak = peak(pe, REDGE, 3.0)
    
    if ipla > 2.7:
        color='r'
    else:
        color = 'b'
    
    alpha = 0.3
    
    rand_plot(ne/1e20, ax=ax_ne, r=REDGE, ls='', marker='o', color=color, mec = 'none', 
              mfc=color, ms=4.0, alpha=alpha)
    rand_plot(te/1e3, ax=ax_te, r=REDGE, ls='', marker='o', color=color, mec = 'none', 
              mfc=color, ms=4.0, alpha=alpha)
    rand_plot(pe/1e4, ax=ax_pe, r=REDGE, ls='', marker='o', color=color, mec = 'none', 
              mfc=color, ms=4.0, alpha=alpha)
    
    #ne.plot(ax=ax_ne, r=REDGE, ls='', marker='.', color=color, mfc=color, alpha=0.2)
    #te.plot(ax=ax_te, r=REDGE, ls='', marker='.', color=color, mfc=color, alpha=0.2)
    #pe.plot(ax=ax_pe, r=REDGE, ls='', marker='.', color=color, mfc=color, alpha=0.2)
    ax_nepeak.plot(t, ne_peak, ls='', marker='o', color=color, mfc=color, mec = 'none', 
                   alpha=alpha, ms=4.0)
    (nbi/1e6).plot(ax=ax_nbi, r=0, ls='-', color=color, alpha=0.2)

ax_ne.plot(0,0, ls='', marker='o', color=color, mec = 'none', 
              mfc='b', ms=4.0, alpha=2*alpha, label = 'Ip = 2.5 MA')
ax_ne.plot(0,0, ls='', marker='o', color=color, mec = 'none', 
              mfc='r', ms=4.0, alpha=2*alpha, label = r'Ip $\geq$ 3.0 MA')
ax_ne.legend(loc=2)
ax_te.set(ylim=[0,1], ylabel='Te @ R = 3.75\n(keV)')
ax_ne.set(ylim=[0,1.6], ylabel='ne @ R = 3.75\n($10^{20}m^{-3}$)')
ax_pe.set(ylim=[0,2], ylabel='pe @ R = 3.75\n($10^4$ Pa)')
ax_nepeak.set(ylim=[0,1.5], ylabel='Ne0/ne_edge')
delupper(ax_te, ax_pe, ax_nepeak)

ax_nbi.set(ylim=[0,35], ylabel='NBI (MW)')    
ax[-1].set(xlim=[48.5,51], xlabel='t (s)')    