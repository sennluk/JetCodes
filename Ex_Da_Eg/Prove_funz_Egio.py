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

shot = 104522
w = ppfs(shot)

ts,Ts = w.hrts.te.get(r=0)    # TS Electron temperature
te,Te = w.ecm1.prfl.get(r=0)       # ECE Mich(KK1) El. temperature
# ter,Ter = w.KK3.prfl.get(r=0)  #  Ece Rad (KK£) el temp



plt.figure('Pippo')
plt.plot(ts,Ts,label='1')
plt.plot(te,Te,label='2')
plt.legend()
plt.show()

