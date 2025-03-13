"""
Created 11/10/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

Partendo da cfr_te_Fish_No1_V04
richiamo singole funzioni da file libreria: myelab.py
Funzioni:
    rprof(shot, rad, tlim1, tlim2, delta): restituisce l'andamenti temporale 
                       di Te per Ece e Ts per la coordinata radiale 'rad'
                       nell'intervallo temporale tra tlim1 e tlim2
    psicalc(shot,tim): restituisce andamento di PSI(R) e Te(PSI) 
                       per lo shoot 'shot'a dato tempo tim. Elimino i 
punti dove psi_ts è minore di zero, perchè corrispondono a punti oltre la linea di vista
Sovrapponedo con linea spessa l'intervallo di psi sul quale si calcola la media'

Divido le funzione che fanno i calcoli da quelle che fanno i plot                       
                       
"""
import numpy as np
from ppfeg import ppfs
import matplotlib.pyplot as plt
import myelab as mye

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse

shot = 104522 # 99950, 99971
tlim1 = 47.5    # limite inferiore selezione tempi - in secondi
tlim2 = 53 
delta = 1     # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.0     # raggio al quale viene fatta l'analisi (in metri)
psi1 = 0.01   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
psi2 = 0.025   # 
rho1 = 0.1
rho2 = 0.2
eP = 0.02  # Relative error assigned to the Ece data 
tlim = (tlim1+tlim2)/2

w = ppfs(shot)
print('JPN = ', shot)
print('t lim1 = ', tlim1)
print('t lim2 = ', tlim2)
#######################################################
# Plot of the EFIT position of the magnetic over time
mye.magax(shot, w, tlim1, tlim2)  # return 1 plot

# Plot Te time trend at R=Rad for Ece-KK1 and HRTS + Errorbars
mye.tprof(shot, w, rad, tlim1, tlim2, delta, eP) # return 1 plot

# PSI profiles at t=tlim and psi over time at r=rad
# and interval of psi1-psi2 evidenced
psiKk1 = mye.psicalc(shot, w)  
mye.psifig(shot, w, psiKk1, tlim, psi1, psi2)  # return 2 plots

# Perform the rho calculation for HRTS and ECE
rhoTs, rhoEce = mye.rhocalc(shot, w, psiKk1, tlim)
mye.rhofig(shot, w, tlim, rho1, rho2, rhoTs, rhoEce)

# Mean values comoputation for Te's in the psi1-psi2 interval
temp_tsM, err_tsM, xm, err_xm, temp_eceM, err_eceM = mye.meancalc(shot, w, tlim1, tlim2, delta, psi1, psi2, eP, psiKk1, rhoEce, rhoTs) # return 1 plot
mye.fig_psi_mean(shot, w, tlim1, tlim2, delta, psi1, psi2, eP, temp_tsM, err_tsM, xm, err_xm, temp_eceM, err_eceM) # return 2 plots

# Resample the faster diagnotic and plot the direct comparison(mean in psi) with errorbars
mye.fig_cfr_psi(shot, w, tlim1, tlim2, delta, psi1, psi2, xm, err_xm, temp_eceM, err_eceM)