"""
Created 11/10/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

Partendo da cfr_te_Fish_No1_V04
richiamo singole funsioni da fiule libreria: myelab.py
Funzioni:
    rprof(shot, rad, tlim1, tlim2, delta): restituisce l'andamenti temporale 
                       di Te per Ece e Ts per la coordinata radiale 'rad'
                       nell'intervallo temporale tra tlim1 e tlim2
    psicalc(shot,tim): restituisce andamento di PSI(R) e Te(PSI) 
                       per lo shoot 'shot'a dato tempo tim. Elimino i 
punti dove psi_ts è minore di zero, perchè corrispondono a punti oltre la linea di vista
Sovreppondo con line aspessa l'intervallo di psi sul quale si calcola la media'
                       
                       
"""
import numpy as np
import ppfeg
import my_flush
from ppfeg import ppfs, jetdata, V2d
import matplotlib.pyplot as plt
from scipy.signal import resample
import myelab as mye

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse

shot = 104522 # 99950, 99971
tlim1 = 47.5    # limite inferiore selezione tempi - in secondi
tlim2 = 53 
delta = 1     # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.0     # raggio al quale viene fatta l'analisi (in metri)
psi1 = 0.01   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
psi2 = 0.02   # 
eP = 0.02  # Relative error assigned to the Ece data 
tlim = (tlim1+tlim2)/2

w = ppfs(shot)
print('JPN = ', shot)
print('t lim1 = ', tlim1)
print('t lim2 = ', tlim2)
#######################################################
# Plot of the EFIT position of the magnetic over time
mye.magax(shot, w, tlim1, tlim2)

# Plot Te time trend at R=Rad for Ece-KK1 and HRTS + Errorbars
mye.rprof(shot, w, rad, tlim1, tlim2, delta, eP)
# return 1 plot

# PSI profiles at t=tlim and psi over time at r=rad
# and interval of psi1-psi2 evidenced
mye.psicalc(shot, w, tlim, psi1, psi2)
# return 2 plots

# Mean values comoputation for Te's in the psi1-psi2 interval
mye.meancalc(shot, w, tlim1, tlim2, delta, psi1, psi2)
# return 1 plot


# print('PSI 1 = ', psi1)
# print('PSI 2 = ', psi2)

#############################################





#############################################


###############################################################################

################################### Resampling and direct comparison - HRTS vs ECE 


# dim = timeTs.size
# v,t = resample(tempEceM, dim, t=timeEce, axis=1)
# tempEceM = v
# timeEce = t

# v,t = resample(errEce, dim, axis=1)  #t=timeEce,
# errEce = v
# timeEce = t

# def retta(x):
#     return x
# x = np.linspace(-1,13,dim)   # Range in keV di dove tracciare la retta
# y = retta(x)  

# plt.figure(f'{shot}_Tts_vs_Tece: Mean values over {psi1}<PSI<{psi2}', clear=True)
# plt.errorbar(tempEceM,xm,lw = linew, xerr = errEceM, yerr = errXm, ecolor='g', elinewidth=.3, label='ECE vs TS')
# plt.plot(x,y,'g--', lw=.8)
# plt.xlim(left=4)
# plt.ylim(bottom=4)
# plt.legend()
# plt.title(f'Shot n.{shot} - Te HRTS vs Te Ece-Michelson (Mean values over psi interval)')
# plt.xlabel('Te Ece-Michelson (keV)')
# plt.ylabel('Te HRTS (keV)')


