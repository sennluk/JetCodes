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

def RATIO(Te, Ts):
    return 550*rhos**1.02 * nus**0.19

shot1 = 104522
w = ppfs(shot1)
Ts1 = w.hrts.te

shot2 = 104521
w = ppfs(shot2)
Ts2 = w.hrts.te



plt.figure('Te Comparison')
plt.plot(Ts1.t,Ts1.v[1,:],label='Te1 - HRTS')
plt.plot(Ts2.t,Ts2.v[1,:],label='Te2 - HRTS')
plt.xlabel('time(s)')
plt.ylabel('Electron temperature (keV)')
plt.title([shot1, shot2])
plt.legend()