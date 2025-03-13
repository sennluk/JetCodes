"""
Created 11/10/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 
V09: Aggiungo la funzione cfrea dizionario variabili                     
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import myelab_V09_no_auto as mye
# import my_manage_file as mym
# %reset

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse

d = {
    'shot' : 99950,  # 99950, 99971, 99970
    'tlim1' : 40, # 47,     # limite inferiore selezione tempi - in secondi
    'tlim2' : 60, # 53, 
    'delta' : 3,      # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.07,      # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.01,   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.08,   # 0.02
    'rho1' : 0.01,# 0.0089,    
    'rho2' : 0.15,# 0.1,
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'win_len' : 15, # window lenght: number of points for the smooth with Savitsky-golay filter
    'deg_pol' : 3, # grado del polinomio usato per lo smooting
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : '/home/lsenni/Python_LS/Cfr_TeTs/Dummy_Folder_pippo/'  # folder to save plot
    }    

#######################################################
# Creo il dizionario con le variabili
vars = mye.dictvar(d)

print('JPN = ', d['shot'])
print('t lim1 = ', d['tlim1'])
print('t lim2 = ', d['tlim2'])

#######################################################
# Multiplot: shot charachteristics
mye.multiplot(d,vars)

# Plot of the EFIT position of the magnetic over time
mye.magax(d,vars)  # return 1 plot

mye.psicalc(d,vars)  
mye.rhocalc(d, vars)

mye.def_range(d,vars) # calcola psi1,2 e rho1,2 al tempo tlim

# Plot Te time trend at R = Rad for Ece-KK1 and HRTS + Errorbars
mye.tprof(d,vars) # return 1 plot

# PSI profiles at t=tlim and psi over time at r=rad
# and interval of psi1-psi2 evidenced
mye.psifig(d,vars)  # return 2 plots

# import myelab_V05 as mye
# Perform the rho calculation for HRTS and ECE
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

# if d['savefigs'] == 1:
#     mye.restore_stdout()
