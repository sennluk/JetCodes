"""
Created 29/01/2025 - lsenni - 
The 'CYCLE version is equal to the MAIN version except that introduce the cycle over 
the shots of interests, creating a Database - (Python dictionary) 
with 'd', and 'vers' correspondig to every shot
Save the DB in .json
"""
import matplotlib.pyplot as plt
import Cfr_Te_LIB_V01 as mye
import Cfr_Te_PLOTS_V01 as graph
import time
import json

# %reset

plt.close('all')  

# Selection of the shots for the Database
shots = [104559, 104560, 104574]# 104520, 104549, 104521, 104522, 104523, 104524, 104525, 104526, 
          # 104547, 104548 ,104549, 104550, 104551, 104553, 104554,104555,104558,104559,104560,
          # 104574, 104575, 104990, 104991,104994] 
# 104995: PSI not present
# 
DB = {} # Empty general database 

saveDB = 1 # 0: don't save, 1: save the DB
timestr = time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation
for shot in shots:
    JPN = shot
    d = {
        'shot' : JPN, # 99950, 99971, 99970
        'Tref' : 1.5,   # (keV) Reference temperature: time instants with Te values above Tref are considered
        'window_size' : 3,   # Numero punti per media mobile calcolo di ti e tf
        'min_increase': 0.0, # (keV) minimum increase in the tiome window to be considered for the ti
        'd2' : 0.08,         # Mezzo intervallo, in cm, sul quale vengono effettuate le medie (rho1, e rho 2)
        'np' : 2,        # numeri di punti di acquisizione che vanno aggiunti per determinare rho1 e rho2
        'fix' : 0.01,    # intervallo aggiuntivo da considerare (in psi e rho) nella differenza tra le due
        'intPsi' : 0.1,    # intervallo in Psi e Rho su cui mediare 
        'delta' : 3,       # tempo in più di tlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
        'rad' : 3.05,      # raggio al quale viene fatta l'analisi (in metri)
        'psi1' : 0.003,    # 0.001 Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
        'psi2' : 0.027,    # 0.027 #0.02
        'eP' : 0.02,       # Relative error assigned to the Ece data
        'win_len' : 15,    # window lenght: number of points for the smooth with Savitsky-golay filter
        'deg_pol' : 3,     # grado del polinomio usato per lo smooting
        'savefigs' : 1,    # 1--> save plots, 0 don't save.
        'mypath' : f'/home/lsenni/Python_LS/Cfr_TeTs/{timestr}_Plots/JPN_{JPN}_Plots/'  # folder to save plot
        }       
    #######################################################
    ## Creo il dizionario con le variabili
    vars = mye.dictvar(d)
    #######################################################
    ## Compute the time window ti-tf for the analysis
    mye.tdef(d,vars)
    
    ## Multiplot: general shot charachteristics
    # graph.multiplot(d,vars)
    
    ## Plot of the EFIT position of the magnetic axis over time
    # graph.magax(d,vars)  # return 1 plot
    
    ## Calculates the PSI and the RHO coordinates of the lines of sights of the two diagnostics
    mye.psicalc(d,vars)   # PSI: Normalized poloidal flux coordinate
    mye.rhocalc(d, vars)  # RHO: Normalized toroidal flux coordinate
    
    ## calcola psi1,2 e rho1,2 prendendo come riferimento la situazione al tempo t=tlim
    ## vedi dettagli nella libreria. Alternativamente si possono inserire manuamente i valori
    mye.def_range_av(d,vars) 
    ## mye.man_range(d,vars)
    
    # Plot Te time trend at R = Rad for Ece-KK1 and HRTS + Errorbars
    # graph.tprof(d,vars) # return 1 plot
    
    ## Perform the plots for evaluating the RHO calculation for HRTS and ECE
    # graph.rhofigs(d, vars)
    
    ## Plot the direct comparison (mean in RHO) with errorbars - based on the slowest diagnostic (HRTS)
    # graph.fig_cfr_rho(d, vars)
    
    # Put data of the shot in the Database DB
    DB[shot] = {'par':d, 'var':vars} 

## Save the Dictionary in .json. If exist do not overwrite
n = len(shots)

import pickle

if saveDB == 1:
    try:
        with open(f"DB_{n}_Shots_V01.pkl", "wb") as f:
            pickle.dump(DB, f)
    except FileExistsError:
        print(" 'DB.json' already exist!")


# if saveDB == 1:    
#     try:
#         with open(f"DB_{n}_Shots.json", "w") as f:    # 'x': non sovrascrive, 'w' sovrascrive
#             json.dump(DB, f, indent=4)     # Scrivi nel file
#     except FileExistsError:
#         print(" 'DB.json' already exist!")

## Per riaprire il dizonario
# with open("DB.json", "r") as f:
#     DB = json.load(f)

# import h5py

# # Salvataggio in formato HDF5
# with h5py.File('DB.h5', 'w') as f:
#     for key, value in DB.items():
#         f.create_dataset(key, data=value)

# # # Caricamento del file HDF5
# # with h5py.File('DB.h5', 'r') as f:
# #     loaded_db = {key: f[key][()] for key in f.keys()}


