"""
Created 11/10/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 
V04: Utilizzo i dizionari per la gestione degli input/putput
d: dizionario dei parametri impostati all'inizio
vars: dizonario della variabili calcolate              
Ogni versione del Main si accoppia con una libreria con lo stesso numero:
l'ultima con la versione della lbreria senza numero versione'
il grafico di confronto è fatto limitatamente ll'intervallo t1-t2
Per fare questo ricalcolo le medie                       
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import myelab_V04 as mye

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse

# Creo dizionario dei parametri
d = {
    'shot' : 104520, # 99950, 99971
    'tlim1' : 48,    # limite inferiore selezione tempi - in secondi
    'tlim2' : 52.5, 
    'delta' : 2,     # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.0,     # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.01,   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.025,   
    'rho1' : 0.1,    
    'rho2' : 0.2,
    'eP' : 0.02 }    # Relative error assigned to the Ece data

shot = d['shot']     
tlim1 = d['tlim1']   
tlim2 = d['tlim2'] 
w = ppfs(shot) 

vars = {
    'tlim' : (tlim1+tlim2)/2,        # tempo medio al quale vengono calcolati i profili per plot di controllo
    'w' : w,                         # canale dati del ppfeg con tutti i dati relativi allo sparo in oggetto
    'tTs' : w.hrts.te,               # Chann Te HRTS  
    'errTs' : w.hrts.dte,            # Chann Errors HRTS.dte
    'tEce' : w.ecm1.prfl,            # Chann Te Ece-Kk1
    'errEce' : (w.ecm1.prfl)*d['eP'],  # Error associated to the Kk1 meas: ab 2%
    'psiTs' : w.hrts.psi,            # Chann PSI HRTS 
    'psiKk1' : None,
    'rhoTs' : None,
    'rhoEce' : None, 
    'temp_tsM' : None, 
    'err_tsM': None, 
    'xm': None, 
    'err_xm' : None, 
    'temp_eceM' : None, 
    'err_eceM' : None, 
    'xm12': None, 
    'err_xm12' : None, 
    'temp_eceM12' : None, 
    'err_eceM12' : None, 
    'temp_tsM_rho' : None,
    'err_tsM_rho' : None,
    'temp_eceM_rho' : None,
    'err_eceM_rho' : None}   

# vars['tlim'] = (tlim1+tlim2)/2  

print('JPN = ', shot)
print('t lim1 = ', tlim1)
print('t lim2 = ', tlim2)
#######################################################
# Plot of the EFIT position of the magnetic over time
mye.magax(d,vars)  # return 1 plot

# Plot Te time trend at R = Rad for Ece-KK1 and HRTS + Errorbars
mye.tprof(d,vars) # return 1 plot

# PSI profiles at t=tlim and psi over time at r=rad
# and interval of psi1-psi2 evidenced
mye.psicalc(d,vars)  
mye.psifig(d,vars)  # return 2 plots

# Perform the rho calculation for HRTS and ECE
mye.rhocalc(d, vars)
mye.rhofig(d, vars)

# Mean values computation for Te's in the psi1-psi2 interval in the (t1-Delta) - (t2 + delta) time interval
mye.meancalc(d, vars) # return 1 plot
mye.fig_psi_mean(d, vars) # return 2 plots

# Mean values computation for Te's in the psi1-psi2 interval in the t1/t2 time interval
mye.meancalc_12(d, vars) 
# Resample the faster diagnostic and plot the direct comparison(mean in psi) with errorbars
mye.fig_cfr_psi(d, vars)

# Mean values computation for Te's in the rho1-rho2 interval in the t1/t2 time interval
mye.mean_calc_rho(d, vars)
# Resample the faster diagnostic and plot the direct comparison(mean in rho) with errorbars
mye.fig_cfr_rho(d, vars)