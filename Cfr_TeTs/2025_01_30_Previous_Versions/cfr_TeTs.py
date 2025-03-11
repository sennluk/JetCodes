#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09/10/2023 by lsenni
"""

import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt

%reset -f        
plt.close('all')  

def Ratio(Te, Ts):
    return Te/Ts

shot = 104522
w = ppfs(shot)
ts,Ts = w.hrts.te.get(r=0)    # TS Electron temperature
te,Te = w.ecm1.prfl.get(r=0)       # ECE Mich(KK1) El. temperature
ter,Ter = w.KK3.prfl.get(r=0)  #  Ece Rad (KK£) el temp

# Ts = Ts1.v[1,]   # Prendo la prima colonna di valori



plt.figure('Te Comparison')
plt.plot(ts,Ts,label='Te - Thomson')
plt.plot(te,Te,label='Te - ECE Michelson')
# plt.plot(Ts2.t,Ts2.v[1,:],label='Te2 - HRTS')
# plt.xlabel('time(s)')
# plt.ylabel('Electron temperature (keV)')
# plt.title(['Shot (1) n.',shot1, 'Shot (1) n.',shot2])
plt.legend()
plt.show()

plt.figure('Te Ratio')
plt.plot(ts,Ratio,label='Te - Thomson')