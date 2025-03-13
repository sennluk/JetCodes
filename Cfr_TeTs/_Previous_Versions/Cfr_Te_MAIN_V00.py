"""
Created 28/01/2025 - lsenni - Based on the Cfr_Te_Fun_V11.py 
(Cfr_Te_Fun_V11.py is the 'previuos MAIN used together with myelab_V12.py,
 containing all the functions called in the MAIN.
The Cfr_Te_Main_V00.py recall various function contained in Cfr_Te_LIB_V00.py 
to copute the comparison of the Electron temperatures measured by ECE and HRTS at JET.
To avoid spikes and similar probles a time interval is uatomatically selected --> see the library
To optimise the comparison the normalized toroidal flux coordinate (RHO) is used.
To optimize the robustness of the comparinson, the range of the RHO values for the average
is computed at every instant, thus taking into account the movement of the plsama center on the poliodal plane
with respet to the lines of sights of the diagnostics. (See libraries for details)
Parameters are scvad in the 'd' dictionary
Variables are stored in the 'vars' dict. 

Averages on the PSI window is no longer implemented,
only RHO

To be added:
    - Choose the ref temperature
    - choose manually the ti and tf: time window for the analysis
    - Coose manually the RHO winndow for the average
- inserisco la possibilità, se selezionata, 
di definire una temperatura di riferimento Tref al disopra della quale si fanno le analisi. 
- E' possibile scegliere se inserire manualmente o calcolare in automaticoi range di psi e rho utilizzati per le medie,
  commentando la funzione def_range (per il calcolo in automatico), si utilizzazo i valori di inzializzazione 
N.B.: the range in meters for the average in RHO is 0.08, as martellato in LIB library
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import Cfr_Te_LIB_V00 as mye
import Cfr_Te_PLOTS_V00 as graph
import time

# import my_manage_file as mym
# %reset

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse
JPN = 104525 #104990 # 104522 # 104525 #104990  #96994
Tref = 6 # keV) Reference temperature: time instants with Te values above Tref are considered
timestr = time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation

d = {
    'shot' : JPN,     # 99950, 99971, 99970
    'Tref' : Tref,
    'window_size' : 5,  # Numero punti per media mobile calcolo di ti e tf
    # 'tlim1' : 47,     # 47,     # limite inferiore selezione tempi - in secondi
    # 'tlim2' : 53,     # 53,     # limite superiore selezione tempi - in secondi
    'np' : 2,       # numeri di punti di acquisizione che vanno aggiunti per determinare rho1 e rho2
    'fix' : 0.01,    # intervallo aggiuntivo da considerare (in psi e rho) nella differenza tra le due
    'intPsi' : 0.1,    # intervallo in Psi e Rho su cui mediare 
    'delta' : 3,      # tempo in più di tlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.05,     # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.003,  # 0.001 Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.027,     # 0.027 #0.02
    # 'rho1' : 0.045,   #0.045,   # 0.0089,    
    # 'rho2' : 0.106,  # 0.106,     # 0.1,
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'win_len' : 15,   # window lenght: number of points for the smooth with Savitsky-golay filter
    'deg_pol' : 3,    # grado del polinomio usato per lo smooting
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : f'/home/lsenni/Python_LS/Cfr_TeTs/{timestr}_JPN_{JPN}_Plots/'  # folder to save plot
    }    

#######################################################
# Creo il dizionario con le variabili
vars = mye.dictvar(d)
#######################################################
# Compute the time window ti-tf for the analysis
mye.tdef(d,vars)

# Multiplot: general shot charachteristics
# graph.multiplot(d,vars)

# Plot of the EFIT position of the magnetic axis over time
# graph.magax(d,vars)  # return 1 plot

# Calculates the PSI and the RHO coordinates of the lines of sights of the two diagnostics
mye.psicalc(d,vars)   # PSI: Normalized poloidal flux coordinate
mye.rhocalc(d, vars)  # RHO: Normalized toroidal flux coordinate

# calcola psi1,2 e rho1,2 prendendo come riferimento la situazione al tempo t=tlim
# vedi dettagli nella libreria. Alternativamente si possono inserire manuamente i valori
mye.def_range_av(d,vars) 
# mye.man_range(d,vars)

# Plot Te time trend at R = Rad for Ece-KK1 and HRTS + Errorbars
# graph.tprof(d,vars) # return 1 plot

# Perform the plots for evaluating the RHO calculation for HRTS and ECE
# graph.rhofigs(d, vars)

# Plot the direct comparison (mean in RHO) with errorbars - based on the slowest diagnostic (HRTS)
graph.fig_cfr_rho(d, vars)
