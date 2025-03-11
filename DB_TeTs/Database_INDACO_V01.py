#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 14:12:23 2023
@author: lsenni
Creazione Database per INDACO
NOTA: Utilizzo il modulo di Edmondo Giovannozzi--> va caricato prima dello script direttamentee da shell

Tutti i dati vengono raccolti in un Dizionario, con ugual dimensioni 
corrispondenti al numero di punti temporali del TS, la diagnostica più lenta.
5/12/2023: Prima prova con ciclo per creare il dizionario che contiene tutti i dati,
le chiavi sono il numero degli shots, i valori sono i dataframe con tutti i dati

OSS: va risolto il problema della mancanza di dati--> il ppfeg si blocca!
Nota: Da apr 2022, shot 100424 --> cambio linea di vista HRTS (+ No Lidar?)
"""

import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt
from scipy.signal import resample
import pandas as pd

# Seleziono gli impulsi di interesse

# Pre-DTE3
pulses_DD_hybrid_all = [96951,96953,96954,96955,96957,96958,96765,96767,96768,96769,96770,96771,96773,
                        96776,96777,96779,96717,96718,96719,96720,96721,96722,
                        96723,96724,96725,96726,96663,96664,96665,96668,96670,96671,96674,96676,96677,
                        96678,96679,96680,96681]
pulses_DD_hybrid = [97451,97452,97455,97457,97458,97460,97589,97591,97594,97604,97605,97683,97684,
                    97685,97686,97732,97734,97735,97736,97737,97740,97741,97742,97779,97780,97781,
                    97783,97784,97785,97786,97787,97788,97790,97791,97796,97797,97798,97799,97844,
                    97847,97852,97896,97898,97976,97977]
pulses_DT_11 = [99815,99817,99818] # Hybrid like? Y.Kazakov experiments
pulses_TT_hybrid = [98913,99162,99163,99164,99273,99274]
pulses_DT_hybrid = [99448,99449,99450,99452,99455,99541,99542,99543,99544,99594,99595,99596,99760,
                    99761,99866,99867,99868,99869,99908,99910,99912,99914,99949,99950,99951,99953] #no LIDAR 99887,99527
# pulses_all_hybrid = np.array(pulses_DD_hybrid_all + pulses_DD_hybrid_bis + pulses_TT_hybrid+pulses_DT_hybrid 
#  + pulses_DT_hybrid)
# pulses_old = np.array(pulses_DD_hybrid_all + pulses_DD_hybrid + pulses_DT_11 + pulses_TT_hybrid + pulses_DT_hybrid)

# Pulses DTE3:

# Pulses per Fishbones
pulses_fi = [96994, 99869, 99801, 99802, 99870, 99950, 99970, 92415, 95986] #, 100793]
    
# Decide on saving file options
save = 2 # 0: non salva, 1:salva , 2: salva .npy
filename = 'Database'
##################
# shot = 104522
pulses = pulses_fi #[99815] #,9917,99818]
ddas = ['hrts','ecm1','kk3','nbi','icrh']  # canali richiamati nel seguito

db = {}
for shot in pulses:
    
    w = ppfs(shot)
    for ddas in w:
        rad = 3 # Position selected
        #########################
        # Slection of the data to be saved
        HRTS = w.hrts.te       # Te by HRTS ( Time, radial positions)
        Errts = w.hrts.dte    # Chann Errors HRTS
        Dens = w.hrts.ne       # Densisty Ne from HRTS
        ErrDens = w.hrts.dne   # Errors on density
        PSI = w.hrts.psi      # EFIT psi HRTS
        KK1 = w.ecm1.prfl     # Te by ECE Michelson (Time, radial positions)
        KK3 = w.kk3.tprf  # Te by ECE Radiometer(Time, radial positions)
        NBI = w.nbi.ptot  # Total power of the NBI beam
        ICRH = w.icrh.ptot # Total power of ICRH
        
    #########################
    # Selection of the radial position of the data 
    idr = np.argmin(abs(HRTS.r - rad))   # Seleziono gli indici dei raggi di interesse (il più vicino)
    Rts = HRTS.r[idr][np.newaxis]        # Selected position
    Tts = HRTS.v[idr,:][np.newaxis,:]    # TS Te in the selected position 
    time_ts = HRTS.t
    RErrts = Errts.r[idr][np.newaxis]     # Posistion of the error
    Errts = Errts.v[idr,:][np.newaxis,:]  # TS Errors on Te in the selected position 
    Ne = Dens.v[idr,:][np.newaxis,:]
    ErrNe = ErrDens.v[idr,:][np.newaxis,:]
    Rpsi = PSI.r[idr][np.newaxis]
    PSI = PSI.v[idr,:][np.newaxis,:]
    
    ErrY = Errts[0,:]/2000   # Errors on Tts divìded by 2 for the +- bar, and /1000 per keV
    ErrY[ErrY>100] = 0       # check on error values --> >1000 wrong value 
    
    ErrN = ErrNe[0,:]/2    # Errors on Ne divìded by 2 for the +- bar
    ErrN[ErrN>1e20] = 0    # Check on error values --> >10^20 wrong value 
    
    idre = np.argmin(abs(KK1.r - rad))
    Rkk1 = KK1.r[idre][np.newaxis]
    Tkk1 = KK1.v[idre,:][np.newaxis,:]
    
    # X mode, op-mode, II harm x mode, 3rd harm x-mode
    
    idre = np.argmin(abs(KK3.r - rad))
    Rkk3 = KK3.r[idre][np.newaxis]
    Tkk3 = KK3.v[idre,:][np.newaxis,:]
    
    idre = np.argmin(abs(NBI.r - rad))
    RPnbi = NBI.r[idre][np.newaxis]
    Pnbi = NBI.v[idre,:][np.newaxis,:]
    
    idre = np.argmin(abs(ICRH.r - rad))
    RPicrh = ICRH.r[idre][np.newaxis]
    Picrh = ICRH.v[idre,:][np.newaxis,:]
    
    ##########################
    # Resample time axixs based on the slowest acq rate
    dim = time_ts.size
    
    Time = KK1.t
    v,t = resample(Tkk1, dim, t=Time, axis=1)
    Tkk1 = v
    
    Time = KK3.t
    v,t = resample(Tkk3, dim, t=Time, axis=1)
    Tkk3 = v
    
    Time = NBI.t
    v,t = resample(Pnbi, dim, t=Time, axis=1)
    Pnbi = v
    
    Time = ICRH.t
    v,t = resample(Picrh, dim, t=Time, axis=1)
    Picrh = v
    
    ##########################
    # Dataframe creation
    
    # Define the columns names
    mycol = ['time','psi','Tts','ErrTS','Ne','ErrNe','Tkk1','Tkk3','Pnbi','Picrh']
    df = pd.DataFrame(columns=mycol)
    # Assign values to the columns
    df.time = time_ts          # Time in sec
    df.psi = PSI.T             # Psi
    df.Tts = Tts.T[:,0]        # TeTS in KeV
    df.ErrTS = ErrY.T          # Error (KeV) estimate for TS measurements as elaborated above
    df.Ne = Ne.T[:,0]          # Density in m-3
    df.ErrNe = ErrN.T          # Error on density
    df.Tkk1 = Tkk1.T[:,0]/1000 #Te KK1 in keV   
    df.Tkk3 = Tkk3.T[:,0]
    df.Pnbi = Pnbi.T[:,0]
    df.Picrh = Picrh.T[:,0]

    # Dictionary voice of the current pulse/shot
    db[shot] = df
    print('Processing pulse '+str(shot))
    
TotShots = len(pulses)
TotCol = len(df.columns)
print('Total number of processed Shots = '+str(TotShots))
print('Total number of columns = '+str(TotCol))
Legend = df.columns
print(Legend)

#####################################################
# Save Dictionary to a .csv file
if save == 1: 
    import csv
            # Open a csv file for writing
    with open("DataBase_op.csv", "w", newline="") as fp:
        # Create a writer object
        writer = csv.DictWriter(fp, fieldnames=db.keys())
        # Write the header row
        writer.writeheader()
        # Write the data rows
        writer.writerow(db)
        print('Done writing dict to a csv file')

#####################################################
# Save Dictionary using NumPy
if save == 2: 
    np.save('database_op.npy', db) 
    
    # Load
# read_dictionary = np.load('my_file.npy',allow_pickle='TRUE').item()
# print(read_dictionary['hello']) # displays "world"