"""
Created 11/10/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 
V05:
introduco il save per i plot: selezionando 1 crea la cartella
e ci salva dentro i plot                      
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import myelab_V05 as mye
import my_manage_file as mym

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse

# Creo dizionario dei parametri
d = {
    'shot' : 104522,  # 99950, 99971
    'tlim1' : 47.5,     # limite inferiore selezione tempi - in secondi
    'tlim2' : 52.5, 
    'delta' : 3,      # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.0,      # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.008,   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.015,   
    'rho1' : 0.01,    
    'rho2' : 0.15,
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : '/home/lsenni/Python_LS/Cfr_TeTs/prova/'  # folder to save plot
    }    

if d['savefigs'] == 1: 
    mym.create_folder(d)
    # mym.save_print_output_to_file(d,vars)
else:
    print('You are not saving the plots')
        
shot = d['shot']     
tlim1 = d['tlim1']   
tlim2 = d['tlim2'] 
w = ppfs(shot) 

vars = {
    'tlim' : (tlim1+tlim2)/2,          # tempo medio al quale vengono calcolati i profili per plot di controllo
    'w' : w,                           # canale dati del ppfeg con tutti i dati relativi allo sparo in oggetto
    'tTs' : w.hrts.te,                 # Chann Te HRTS  
    'errTs' : w.hrts.dte,              # Chann Errors HRTS.dte
    'tEce' : w.ecm1.prfl,              # Chann Te Ece-Kk1
    'errEce' : (w.ecm1.prfl)*d['eP'],  # Error associated to the Kk1 meas: ab 2%
    'psiTs' : w.hrts.psi,              # Chann PSI HRTS 
    'psiKk1' : None,
    'time_ts' : None,
    'time_ece' : None,
    'rho_ts' : None,
    'rho_ece' : None, 
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

# import myelab_V05 as mye
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

if d['savefigs'] == 1:
    mye.restore_stdout()