"""
Created 21/06/2024 - lsenni
Partendo dal cfr_Te_Fun_v10
Voglio calcolare e confrontare i valori di Rho e PSI
per un dato sparo         
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import myelab_V10 as mye
# import my_manage_file as mym
# %reset

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse
JPN = 96994
Tref = 7 # keV) Reference temperature: time instant whiht Te values above Tref are considered

d = {
    'shot' : JPN,     # 99950, 99971, 99970
    'Tref' : Tref,
    'tlim1' : 43,     # 47,     # limite inferiore selezione tempi - in secondi - inizializzazione
    'tlim2' : 57,     # 53,   # inizializzazione
    'delta' : 3,      # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.07,     # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.001,   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.9,     # 0.02
    'rho1' : 0.001,   # 0.0089,    
    'rho2' : 0.9,     # 0.1,
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'win_len' : 15,   # window lenght: number of points for the smooth with Savitsky-golay filter
    'deg_pol' : 3,    # grado del polinomio usato per lo smooting
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : '/home/lsenni/Python_LS/Cfr_TeTs/Dummy_Folder/'  # folder to save plot
    }    

#######################################################
# Creo il dizionario con le variabili
vars = mye.dictvar(d)

print('JPN = ', d['shot'])
print('t lim1 = ', d['tlim1'])
print('t lim2 = ', d['tlim2'])
#######################################################
# Define tlim1 e tlim2
mye.tdef(d,vars)
print('JPN = ', d['shot'])
print('Tref = ', d['Tref'])
print('t lim1 = ', d['tlim1'])
print('t lim2 = ', d['tlim2'])

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
mye.rhofig2(d, vars)
