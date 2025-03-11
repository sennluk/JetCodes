#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 14:12:23 2023
@author: lsenni
Creazione Database per INDACO
NOTA: Utilizzo il modulo di Edmondo Giovannozzi--> va caricato prima dello script direttamentee dal shell

Tutti i dati vengono raccolti in un DataFrame, con ugual dimensioni 
corrispondenti al numero di punti temporali del TS, la diagnostica più lenta.
5/12/2023: Prima prova con ciclo per creare il dizionario che contiene tutti i dati,
le chiavi sono il numero degli shots, i valori sono i dataframe con tutti i dati
  
12/12: Le diverse finestre temporali di acquisizione delle diverse diagnostiche causano alcuni problemi.
--> Modifico provando ad introdurre un intervallo temporale fisso di acquisizione: t1-t2 per tutte le diagnostiche
Scelto in base all'esistenza dei dati delle singole diagnostiche.

Da fare:
    Medie su più posizioni radiali?
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
                    99761,99866,99867,99868,99869,99908,99910,99912,99914,99949,99950,99951,99953] 
#no LIDAR 99887,99527
# pulses_all_hybrid = np.array(pulses_DD_hybrid_all + pulses_DD_hybrid_bis + pulses_TT_hybrid+pulses_DT_hybrid + pulses_DT_hybrid)
# pulses_old = np.array(pulses_DD_hybrid_all + pulses_DD_hybrid + pulses_DT_11 + pulses_TT_hybrid + pulses_DT_hybrid)

# Pulses DTE (starting at shot n.?):
    
P00 = [95679, 94700] # Kazakov fast ions    
P01 = [99801, 100793] # Kiptily
P02 = [99869, 99643] # Fishbones
P03 = [99950, 96994] # MHD 49 - 54 sec
P04 = [99971, 99970] # Record
P05 = [104522, 104523, 104524, 104525, 104526] # Non thermal -Tritium rich
P06 = [104547,104548,104549,104550,104551] # ICRH modulation - up to 12 KeV
P07 = [104553,104554,104555,104556,104557,104558,104559,104560] # Neon Seeded
P08 = [104574, 104575, 104580] # Hybrid High field
P09 = [104990,104991, 10499, 4104495, 104991,104994, 104502, 104511] # And more: see text

##################
# pulses = [99981]
# pulses = [99815, 99817, 99818]
pulses = [99815]
# pulses = np.array(P05 + P06)
db = {}
for shot in pulses:
    
    w = ppfs(shot)
    # t1 = 40
    # t2 = 60
    rad = 3 # Position selected
    #########################
    # Slection of the data to be saved
    HRTS = w.hrts.te       # Te by HRTS (Time, radial positions)
    Errts = w.hrts.dte    # Chann Errors HRTS
    Dens = w.hrts.ne       # Density Ne from HRTS
    ErrDens = w.hrts.dne   # Errors on density
    PSI = w.hrts.psi      # EFIT psi HRTS
    KK1 = w.ecm1.prfl     # Te by ECE Michelson (Time, radial positions)
    KK3 = w.kk3.tprf  # Te by ECE Radiometer(Time, radial positions)
    NBI = w.nbi.ptot  # Total power of the NBI beam
    ICRH = w.icrh.ptot # Total power of ICRH
    
    # To choose the time interval for each shot
    t1 = np.max([np.min(HRTS.t),np.min(KK1.t), np.min(KK3.t), np.min(NBI.t), np.min(ICRH.t)])
    t2 = np.min([np.max(HRTS.t),np.max(KK1.t), np.max(KK3.t), np.max(NBI.t), np.max(ICRH.t)])
    
    #########################
    # Selection of time interval and the radial position of the data 
    
    idt = (HRTS.t >= t1) & (HRTS.t <= t2) # Seleziono intervallo temporale
    TimeTS = HRTS.t[idt]
    TempTS = HRTS.v[:,idt]
    Errts.t = Errts.t[idt]
    Errts.v = Errts.v[:,idt]
    DensTS = Dens.v[:,idt]
    ErrDensTS = ErrDens.v[:,idt]
    PSIefit = PSI.v[:,idt]
    
    idr = np.argmin(abs(HRTS.r - rad))   # Indices of the selected postions
    Rts = HRTS.r[idr][np.newaxis]        # Selected position
    Tts = TempTS[idr,:][np.newaxis,:]    # TS Te in the selected position 
    time_ts = TimeTS                     # Time axis from HRTS diagn
    RErrts = Errts.r[idr][np.newaxis]     # Posistion of the error
    Errts = Errts.v[idr,:][np.newaxis,:]  # TS Errors on Te in the selected position 
    Ne = DensTS[idr,:][np.newaxis,:]
    ErrNe = ErrDensTS[idr,:][np.newaxis,:]
    Rpsi = PSI.r[idr][np.newaxis]
    PSI = PSIefit[idr,:][np.newaxis,:]
    
    ErrY = Errts[0,:]/2000   # Errors on Tts divìded by 2 for the +- bar, and /1000 per keV
    ErrY[ErrY>100] = 0       # check on error values --> >1000 wrong value 
    
    ErrN = ErrNe[0,:]/2    # Errors on Ne divìded by 2 for the +- bar
    ErrN[ErrN>1e20] = 0    # Check on error values --> >10^20 wrong value 
    
    idte = (KK1.t >= t1) & (KK1.t <= t2)
    TimeKK1 = KK1.t[idte]
    TempKK1 = KK1.v[:,idte]
    idre = np.argmin(abs(KK1.r - rad))
    Rkk1 = KK1.r[idre][np.newaxis]
    Tkk1 = TempKK1[idre,:][np.newaxis,:]
    
    # X mode, op-mode, II harm x mode, 3rd harm x-mode
    
    idte = (KK3.t >= t1) & (KK3.t <= t2)
    TimeKK3 = KK3.t[idte]
    TempKK3 = KK3.v[:,idte]
    idre = np.argmin(abs(KK3.r - rad))
    Rkk3 = KK3.r[idre][np.newaxis]
    Tkk3 = TempKK3[idre,:][np.newaxis,:]
    
    idte = (NBI.t >= t1) & (NBI.t <= t2)
    TimeNBI = NBI.t[idte]
    Pnbi = NBI.v[:,idte]
    idre = np.argmin(abs(NBI.r - rad))
    RPnbi = NBI.r[idre][np.newaxis]
    Pnbi = Pnbi[idre,:][np.newaxis,:]
    
    idte = (ICRH.t >= t1) & (ICRH.t <= t2)
    TimeICRH = ICRH.t[idte]
    Picrh = ICRH.v[:,idte]
    idre = np.argmin(abs(ICRH.r - rad))
    RPicrh = ICRH.r[idre][np.newaxis]
    Picrh = Picrh[idre,:][np.newaxis,:]
    
    ##########################
    # Resample time axixs based on the slowest acq rate
    dim = time_ts.size
    
    Time = TimeKK1
    v,t = resample(Tkk1, dim, t=Time, axis=1)
    Tkk1 = v
    TimeKK1 = t
    
    Time = TimeKK3
    v,t = resample(Tkk3, dim, t=Time, axis=1)
    Tkk3 = v
    timeKK3 = t
    
    Time = TimeNBI
    v,t = resample(Pnbi, dim, t=Time, axis=1)
    Pnbi = v
    
    Time = TimeICRH
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
    df.Tts = Tts.T[:,0]/1000        # TeTS in KeV
    df.ErrTS = ErrY.T          # Error (KeV) estimate for TS measurements as elaborated above
    df.Ne = Ne.T[:,0]          # Density in m-3
    df.ErrNe = ErrN.T          # Error on density
    df.Tkk1 = Tkk1.T[:,0]/1000 # Te KK1 in keV   
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
print('legend:', Legend)

# test plots
# fig,(ax0,ax1) = plt.subplots(nrows=2,sharex=True)
# ax0.plot(df.time,df.Tts)
# ax1.plot(TimeKK1,df.Tkk1)

# plt.figure()
# plt.scatter(db[99815].Tkk1,db[99815].Tkk3)

# shotn = 99817
# plt.figure()
# plt.scatter(db[shotn].Tts,db[shotn].Tkk1)
#####################################################
# Save Dictionary using NumPy

# np.save('DB_12_12.npy', db) 

# Load
# read_dictionary = np.load('my_file.npy',allow_pickle='TRUE').item()
# print(read_dictionary['hello']) # displays "world"

#####################################################
# Save Dictionary to a .csv file
# import csv
#        # Open a csv file for writing
# with open("DataBase.csv", "w", newline="") as fp:
#     # Create a writer object
#     writer = csv.DictWriter(fp, fieldnames=db.keys())

#     # Write the header row
#     writer.writeheader()

#     # Write the data rows
#     writer.writerow(db)
#     print('Done writing dict to a csv file')
