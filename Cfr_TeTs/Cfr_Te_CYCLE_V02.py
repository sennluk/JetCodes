"""
Created 29/01/2025 - lsenni - 
The 'CYCLE version is equal to the MAIN version except that introduce the cycle over 
the shots of interests, creating a Database - (Python dictionary) 
with 'd', and 'vers' correspondig to every shot
Save the DB in .json
"""
import numpy as np
import matplotlib.pyplot as plt
import Cfr_Te_LIB_V02 as mye
import Cfr_Te_PLOTS_V02 as graph
import time
import json
import pickle
# %reset

plt.close('all')  

# Selection of the shots for the Database
shots =  [104520, 104521, 104522, 104523, 104524, 104525, 104526, 
            104547, 104548 ,104549, 104550, 104551, 104553, 104554,104555,104558,104559,104560,
            104574, 104575, 104990, 104991,104994] 
# 104995: PSI not present
# 
DB = {} # Empty general database 

saveDB = 1 # 0: don't save, 1: save the DB
timestr = time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation

fig,ax=plt.subplots(1, num = 'ECE vs HRTS')
for shot in shots:
    # plt.close('all')  
    JPN = shot
    d = {
        'shot' : JPN,        # 99950, 99971, 99970
        'Tref' : 1,          # (keV) Reference temperature: time instants with Te values above Tref are considered
        'window_size' : 10,  # Number of points for the moving average used for ti and tf computation
        'min_increase': 0.0, # Minimum increase in the time window to be considered for the ti 
        'd2' : 0.08,         # Half interval (in cm) around the plasma center used to average values of HRTS and ECE 
        'np' : 1,            # Number of acquisition data points to add in computing the RHO ranges
        'fix' : 0.01,        # Delta to be safer in the RHO and PSI ranges determination
        'intPsi' : 0.1,      # PSI ramnge for the average  
        'delta' : 3,         # Time (in sec) to be added to the ti-tf interval for multiplosts figure 
        'rad' : 3.05,        # Major radius coordinate (in m) for the single position comparison 
        'psi1' : 0.003,      # Smallest PSI value for the PSI range
        'psi2' : 0.027,      # Highest PSI value for the PSI range
        'eP' : 0.02,         # Relative error assigned to the Ece data
        'win_len' : 15,      # Window lenght: number of points for the smooth with Savitsky-golay filter
        'deg_pol' : 3,       # Poly degree used for the smooting 
        'savefigs' : 0,      # 1--> save plots, 0 don't save.
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
    ## Multiplot: general shot charachteristics
    # graph.multiplot(d,vars)

    ## Plot of the EFIT position of the magnetic axis over time
    # graph.magax(d,vars)  # return 1 plot

    # Plot Te time trend at R = Rad for Ece-KK1 and HRTS + Errorbars
    # graph.tprof(d,vars) # return 1 plot

    ## Perform the plots for evaluating the PSI calculation for HRTS and ECE at t=tlim
    ##  PSI values are manually inserted
    # graph.psifig(d,vars)

    ## Perform the plots for evaluating the RHO calculation for HRTS and ECE at t=tlim
    ##  RHO values are automatically evaluated and contained in 'ranges' array
    # graph.rhofig(d, vars)

    ## Plot the direct comparison (mean in RHO) with errorbars - based on the slowest diagnostic (HRTS)
    # graph.fig_cfr_rho(d, vars)
    
    ######################
    # shot = d['shot']
    # tlim1 = vars['ti']
    # tlim2 = vars['tf']
    # temp_tsM_rho = vars['temp_tsM_rho']
    # err_tsM_rho = vars['err_tsM_rho']
    # temp_eceM_rho = vars['temp_eceM_rho']
    # err_eceM_rho = vars['err_eceM_rho']
    # timeTs2 = vars['timeTs2']
    # timeEce22 = vars['timeEce22']
    # ranges = vars ['ranges']
    # ###
    
             
    # left = min(np.nanmin(temp_eceM_rho), np.nanmin(temp_tsM_rho))
    # right  = max(np.nanmax(temp_eceM_rho), np.nanmax(temp_tsM_rho)) 
    # def retta(x):
    #     return x
    
    # x = np.linspace(left-0.5,right+0.5,timeTs2.shape[0])   # Range in keV di dove tracciare la retta
    # y = retta(x)  
    
    # # fig,ax=plt.subplots(1, num = 'ECE vs HRTS')
    # # ax.scatter(temp_eceM_rho, temp_tsM_rho,label='Rho Averaged ECE vs TS ') 
    # ax.plot(x,y,'g--', lw=.8)
    # ax.scatter(temp_tsM_rho, temp_eceM_rho, marker='o', s=3, color = 'blue', facecolors='none', edgecolors='r', linewidth=.2, label='ECE vs TS')
    # # ax.errorbar(temp_tsM_rho, temp_eceM_rho, xerr = err_tsM_rho, yerr = err_eceM_rho, 
    # #               marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5, label='ECE vs TS')
    # ax.set_title(f'JPN {shot} - Te Ece-Michelson vs Te HRTS  for {tlim1:.2f}<t<{tlim2:.2f} (s)') #:.2f per avere 2 cifre decimali
    # ax.set_xlabel('Te HRTS (keV)')
    # ax.set_ylabel('Te Ece-Michelson (keV)')
    # ax.legend()
    
    
    
    # Put data of the shot in the Database DB
    DB[shot] = {'par':d, 'var':vars} 

## Save the Dictionary in .json. If exist do not overwrite
n = len(shots)

if saveDB == 1:
    try:
        with open(f"DB_{n}_Shots_V02.pkl", "wb") as f:
            pickle.dump(DB, f)
    except FileExistsError:
        print(" 'DB.pkl' already exist!")

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


