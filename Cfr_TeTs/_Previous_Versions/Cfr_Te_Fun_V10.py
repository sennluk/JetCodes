"""
Created 20/06/2024 - lsenni - partendo dalle versioni precedenti
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

Definizione automatica dei  limiti di PSI e RHO:
    i valori sono inizializzati nel dizionario d
    poi messi nel dizionarion vars

V09: Aggiungo la funzione crea dizionario variabili         
V10: 
- definisco il JPN all'inizio
- inserisco la possibilità, se selezionata, 
di definire una temperatura di riferimento Tref al disopra della quale si fanno le analisi. 
- La definizione degli intervalli di psi e rho avviene in modo automatico:
    
- Aggiungo multiplot
- Aggiungo nomi figure        
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import myelab_V10 as mye
import time

# import my_manage_file as mym
# %reset

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse
JPN = 104522 # 104990  #96994
Tref = 6 # keV) Reference temperature: time instants with Te values above Tref are considered
timestr = time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation

d = {
    'shot' : JPN,     # 99950, 99971, 99970
    'Tref' : Tref,
    'tlim1' : 47,     # 47,     # limite inferiore selezione tempi - in secondi - inizializzazione
    'tlim2' : 53,     # 53,     # inizializzazione
    'interv' : 0.1,    # intervallo sul quale si andrà a mediare, in m, determinando i valori di psi e rho
    'delta' : 3,      # tempo in più di tlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.05,     # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.0001,  # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.9,     # 0.02
    'rho1' : 0.001,   # 0.0089,    
    'rho2' : 0.9,     # 0.1,
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'win_len' : 15,   # window lenght: number of points for the smooth with Savitsky-golay filter
    'deg_pol' : 3,    # grado del polinomio usato per lo smooting
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : f'/home/lsenni/Python_LS/Cfr_TeTs/{timestr}_JPN_{JPN}_Plots/'  # folder to save plot
    }    

#######################################################
# Creo il dizionario con le variabili
vars = mye.dictvar(d)

print('JPN = ', d['shot'])
print('t lim1 = ', d['tlim1'])
print('t lim2 = ', d['tlim2'])

#######################################################
# Define tlim1 e tlim2: if mye.tdef is selecrted the range is automaticcaly selected based on the temperature value Tref

# mye.tdef(d,vars)
# print('Tref = ', d['Tref'])
# print('t lim1 = ', d['tlim1'])
# print('t lim2 = ', d['tlim2'])

# Multiplot: shot charachteristics
mye.multiplot(d,vars)

# Plot of the EFIT position of the magnetic axis over time
mye.magax(d,vars)  # return 1 plot

# Calculates the PSI and the RHO coordinates of the lines of sights of the two diagnostics
mye.psicalc(d,vars)   # PSI: Normalized poloidal flux coordinate
mye.rhocalc(d, vars)  # RHO: Normalized toroidal flux coordinate

# calcola psi1,2 e rho1,2 prendendo comne riferimento la situazione al tempo t=tlim
# definisce ilminimo e prende un intervallo pari a 2*interv in m
mye.def_range(d,vars) 

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
