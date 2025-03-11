"""
Created on 28/01/2025 
@author: lsenni. Prova di calcolo dei residui
"""
import numpy as np
from ppfeg import ppfs
import matplotlib.pyplot as plt

plt.close('all')

shot = 104522
w = ppfs(shot)

sig4 = w.hrts.sig4  # Segnali misurati dal canale 4 dell'HRTS
fit4 = w.hrts.fit4  # Valori Fit del canale 4 HRTS

di = sig4-fit4
# di = np.zeros_like(sig4)
# for i in sig4.t:
#     di[:,i] = sig4[:,i] - fit4[:,i]

# Plot segnali misurati vs posizione
plt.figure('1')
plt.plot(sig4.r,sig4.v)
plt.title('Measured HRTS Signals Ch4 vs positions')
plt.xlabel('R (m)')
plt.ylabel('HRTS Ch4 Signal')

# Plot differenze vs posizione
plt.figure('2')
plt.plot(sig4.r, di.v)
plt.title('Sig - Fit ch4 vs positions')
plt.xlabel('R (m)')
plt.ylabel('Differences')    

ind = 250
time = round(sig4.t[ind],2)
# Plot andamento segnale e fit a dato tempo
plt.figure('3')
plt.plot(sig4.r,sig4.v[:,ind], label='sig4')
plt.plot(fit4.r, fit4.v[:,ind], label = 'fit4')
plt.legend()
plt.title(f'Signal and Fit at t = {time} vs positions')
plt.xlabel('R (m)')
plt.ylabel('Differences')  

idx = 30
pos = round(sig4.r[idx],2)
# Plot andamento segnale nel tempo a data posizione
plt.figure('4')
plt.plot(sig4.t,di.v[idx,:])
plt.title(f'Difference at R = {pos} vs time')
plt.xlabel('time (s)')
plt.ylabel('Differences')  