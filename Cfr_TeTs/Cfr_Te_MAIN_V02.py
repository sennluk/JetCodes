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
V01: partendo dalla V=00 (e librerie legate) funzionante, vado eliminando cose superflue

"""
import matplotlib.pyplot as plt
import Cfr_Te_LIB_V02 as mye
import Cfr_Te_PLOTS_V02 as graph
import time

# %reset

plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot
JPN = 104522 #104990 # 104522 # 104525 #104990  #96994
# DTE3 shot list: 104990,991,994,995,999   
# 104520,521,522,523,524,526

timestr = time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation

d = {
    'shot' : JPN,        # 99950, 99971, 99970
    'Tref' : 1,          # (keV) Reference temperature: time instants with Te values above Tref are considered
    'window_size' : 10,  # Number of points for the moving average used for ti and tf computation
    'min_increase': 0.0, # Minimum increase in the time window to be considered for the ti 
    'd2' : 0.08,         # Half interval (in cm) around the plasma center used to average values of HRTS and ECE 
    'np' : 1,            # Number of acquisition data points to add in computing the RHO ranges
    'fix' : 0.01,        # Delta to be safer in the RHO and PSI ranges determination
    'intPsi' : 0.1,      # PSI range for the average  
    'delta' : 3,         # Time (in sec) to be added to the ti-tf interval for multiplosts figure 
    'rad' : 3.05,        # Major radius coordinate (in m) for the single position comparison 
    'psi1' : 0.003,      # Smallest PSI value for the PSI range
    'psi2' : 0.027,      # Highest PSI value for the PSI range
    'eP' : 0.03,         # Relative error assigned to the Ece data
    'win_len' : 15,      # Window lenght: number of points for the smooth with Savitsky-golay filter
    'deg_pol' : 3,       # Poly degree used for the smoothing 
    'savefigs' : 1,      # 1--> save plots, 0 don't save.
    'mypath' : f'/home/lsenni/Python_LS/Cfr_TeTs/{timestr}_JPN_{JPN}_Plots/'  # folder to save plots
    }    

#######################################################
## Creo il dizionario con le variabili
vars = mye.dictvar(d)
#######################################################
## Compute the time window ti-tf for the analysis
mye.tdef(d,vars)
print(' ###   tdef ok!!!!  ####')
## Calculates the PSI and the RHO coordinates of the lines of sights of the two diagnostics
mye.psicalc(d,vars)   # PSI: Normalized poloidal flux coordinate
mye.rhocalc(d, vars)  # RHO: Normalized toroidal flux coordinate
print('### PSI e RHO calc ok !!!!! ###')
## calcola psi1,2 e rho1,2 prendendo come riferimento la situazione al tempo t=tlim
## vedi dettagli nella libreria. Alternativamente si possono inserire manuamente i valori
mye.def_range_av(d,vars) 
# mye.man_range(d,vars)
print('### Def range ok !!!!! ###')

graph.multiplot(d,vars)       ## Multiplot: general shot charachteristics
graph.tprof(d,vars)        ## Plot Te time trend at R = Rad for Ece-KK1 and HRTS + Errorbars
# graph.magax(d,vars)        ## Plot of the EFIT position of the magnetic axis over time
# graph.psifig(d,vars)       ## Perform the plots for evaluating the PSI (manual) calculation for HRTS and ECE at t=tlim
# graph.rhofig(d, vars)      ## Perform the plots for evaluating the RHO (automatic) calculation for HRTS and ECE at t=tlim
# graph.fig_cfr_rho(d, vars) ## Plot the direct comparison (mean in RHO) with errorbars - based on the slowest diagnostic (HRTS)








