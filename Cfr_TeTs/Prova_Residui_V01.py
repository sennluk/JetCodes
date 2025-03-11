"""
Created on 28/01/2025 
@author: lsenni. Prova di calcolo dei residui
"""
import numpy as np
from ppfeg import ppfs
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.close('all')

shot = 104522
w = ppfs(shot)

sig1 = w.hrts.sig1  # Segnali misurati dal canale 4 dell'HRTS
fit1 = w.hrts.fit1  # Valori Fit del canale 4 HRTS
sig2 = w.hrts.sig2  # Segnali misurati dal canale 4 dell'HRTS
fit2 = w.hrts.fit2  # Valori Fit del canale 4 HRTS
sig3 = w.hrts.sig3  # Segnali misurati dal canale 4 dell'HRTS
fit3 = w.hrts.fit3  # Valori Fit del canale 4 HRTS
sig4 = w.hrts.sig4  # Segnali misurati dal canale 4 dell'HRTS
fit4 = w.hrts.fit4  # Valori Fit del canale 4 HRTS

csqr = w.hrts.csqr  # Fit ChiSquare 

di1 = sig1 - fit1
di2 = sig2 - fit2
di3 = sig3 - fit3
di4 = sig4 - fit4

inst = 50  # Istante temporale in secondi
pos = 3    # Posizione in metri per il calcolo a singola posizione

idx_time = np.argmin(abs(sig3.t - inst)) # Index of the selected time instant
idx_pos = np.argmin(abs(sig3.r - pos))   # Index of the selected position
t_inst = sig3.t[idx_time]                # Exact time instant
pos = sig3.r[idx_pos]                    # Exact position 
 
print('Selected time =',t_inst)
print(f'Selected position = {pos:.2f}')  # 2 cifre decimali anche senza il round
print(f'Selected time = {t_inst:.2f}')   # 2 cifre decimali anche senza il round

#%%   # Sezione dei plots

# # Plot segnali misurati vs posizione
# plt.figure('1 - Ch3 - Sig vs Pos')
# plt.plot(sig3.r,sig3.v)
# plt.title('Measured HRTS Signals Ch3 vs positions')
# plt.xlabel('R (m)')
# plt.ylabel('HRTS Ch3 Signal')

# plt.figure('1 - Ch4 - Sig vs Pos')
# plt.plot(sig4.r,sig4.v)
# plt.title('Measured HRTS Signals Ch4 vs positions')
# plt.xlabel('R (m)')
# plt.ylabel('HRTS Ch4 Signal')

# Plot Residuals vs position 
plt.figure('2 -Res vs Pos')
plt.plot(sig1.r, di1.v, color = 'c')
plt.plot(sig2.r, di2.v, color = 'g')
plt.plot(sig3.r, di3.v, color = 'b')
plt.plot(sig4.r, di4.v, color = 'r')
plt.title('Residuals Ch1,2,3,4 vs position')
plt.xlabel('R (m)')
plt.ylabel('Residuals')  
cyan_line =  mlines.Line2D([], [], color='c', label='Ch1')
green_line = mlines.Line2D([], [], color='g', label='Ch2') 
blue_line = mlines.Line2D([], [], color='b', label='Ch3')
red_line = mlines.Line2D([], [], color='r', label='Ch4')
plt.legend(handles=[cyan_line, green_line, blue_line, red_line])

# Plot andamento segnale e fit a dato tempo
plt.figure('3 - Sig e Fit vs Pos')
plt.plot(sig1.r,sig1.v[:,idx_time], 'c-', label='sig1')
plt.plot(fit1.r, fit1.v[:,idx_time], 'c--', label = 'fit1')
plt.plot(sig2.r,sig2.v[:,idx_time], 'g-', label='sig2')
plt.plot(fit2.r, fit2.v[:,idx_time],'g--', label = 'fit2')
plt.plot(sig3.r,sig3.v[:,idx_time], 'b-', label='sig3')
plt.plot(fit3.r, fit3.v[:,idx_time], 'b--', label = 'fit3')
plt.plot(sig4.r,sig4.v[:,idx_time], 'r-', label='sig4')
plt.plot(fit4.r, fit4.v[:,idx_time], 'r--', label = 'fit4')
plt.legend()
plt.title(f'Signal and Fit at t = {t_inst:.2f} s vs positions')
plt.xlabel('R (m)')
plt.ylabel('A.U.')  

# Plot andamento segnale nel tempo a data posizione
plt.figure('4 - Res vs time')
plt.plot(sig1.t,di1.v[idx_pos,:], label='Ch1')
plt.plot(sig2.t,di2.v[idx_pos,:], label='Ch2')
plt.plot(sig3.t,di3.v[idx_pos,:], label='Ch3')
plt.plot(sig4.t,di4.v[idx_pos,:], label='Ch4')
plt.title(f'Residuals at R = {pos:.2f} m vs time')
plt.xlabel('time (s)')
plt.ylabel('Residuals')
plt.legend()  

# Plot andamento residuals vs Fit 
plt.figure('5 - Residuals vs Fit at t')
plt.scatter(fit1.v[:,idx_time], di1.v[:,idx_time], marker = '1', label='Ch1')
plt.scatter(fit2.v[:,idx_time], di2.v[:,idx_time], marker = '2', label='Ch2')
plt.scatter(fit3.v[:,idx_time], di3.v[:,idx_time], marker = '3', label='Ch3')
plt.scatter(fit4.v[:,idx_time], di4.v[:,idx_time], marker = '4', label='Ch4')
plt.axhline(y=0, ls = '--', lw = 0.8)
plt.title(f'Residuals vs Fit at t = {t_inst:.2f} s')
plt.xlabel('Fitted value')
plt.ylabel('Residuals')
plt.legend()  

# Plot Chisquare Trend at t
plt.figure('6 - Chi^2 at t')
plt.plot(csqr.r, csqr.v[:,idx_time], label='Chi^2')
plt.axhline(y = 1, ls='--', lw = 0.8)
plt.title(f'Chi^2 at t = {t_inst:.2f} s')
plt.xlabel('Position (m)')
plt.ylabel('Chi^2')
plt.legend()  

# Plot Chisquare Trend at R
plt.figure('7 - Chi^2 at R')
plt.plot(csqr.t, csqr.v[idx_pos,:], label='Chi^2')
plt.axhline(y = 1, ls='--', lw = 0.8)
plt.title(f'Chi^2 at R = {pos:.2f} m')
plt.xlabel('time (s)')
plt.ylabel('Chi^2')
plt.legend()  








