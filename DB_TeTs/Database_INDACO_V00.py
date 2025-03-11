#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 14:12:23 2023
@author: lsenni
Prova creazione Database per INDACO
NOTA: Utilizzo il modulo di Edmondo Giovannozzi--> va caricato prima dello script direttamentee dal shell
Data List:
    

"""

import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt
from scipy.signal import resample
import pandas as pd

# Seleziono gli impulsi di interesse

pulses_DT_11 = [99815,99817,99818] # Hybrid like? Y.Kazakov experiments
pulses_DD_hybrid_all = [96951,96953,96954,96955,96957,96958,96765,96767,96768,96769,96770,96771,96773,
                        96776,96777,96779,96717,96718,96719,96720,96721,96722,
                        96723,96724,96725,96726,96663,96664,96665,96668,96670,96671,96674,96676,96677,
                        96678,96679,96680,96681]
pulses_DD_hybrid_bis = [97451,97452,97455,97457,97458,97460,97589,97591,97594,97604,97605,97683,97684,
                        97685,97686,97732,97734,97735,97736,97737,
		97740,97741,97742,97779,97780,97781,97783,97784,97785,97786,97787,97788,97790,97791,97796,97797,
        97798,97799,97844,97847,97852,
		97896,97898,97976,97977]
pulses_TT_hybrid = [98913,99162,99163,99164,99273,99274]
pulses_DT_hybrid = [99448,99449,99450,99452,99455,99541,99542,99543,99544,99594,99595,99596,99760,99761,
                    99866,99867,
                    99868,99869,99908,99910,99912,99914,99949,99950,99951,99953] #no LIDAR 99887,99527

pulses_all_hybrid = np.array(pulses_DD_hybrid_all+pulses_DD_hybrid_bis+pulses_TT_hybrid+pulses_DT_hybrid)

##################
shot = 104522
w = ppfs(shot)
rad = 3 # Position selected
#########################
# Slection of the data to be saved
HRTS = w.hrts.te       # Te by HRTS ( Time, radial positions)
Errts = w.hrts.dte    # Chann Errors HRTS
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
Rpsi = PSI.r[idr][np.newaxis]
PSI = PSI.v[idr,:][np.newaxis,:]

ErrY = Errts[0,:]/2000 # Errors divìded by 2 for the +- bar, and /1000 per keV
ErrY[ErrY>100] = 0       # check on error values --> >1000 wrong value 

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
v,t = resample(Tkk1, dim, t=time_ts, axis=1)
Tkk1 = v
# TimeECE = t
v,t = resample(Tkk3, dim, t=time_ts, axis=1)
Tkk3 = v

v,t = resample(Pnbi, dim, t=time_ts, axis=1)
Pnbi = v

v,t = resample(Picrh, dim, t=time_ts, axis=1)
Picrh = v

##########################
# Dataframe creation

# Define the columns names
mycol = ['time','psi','Tts','ErrTS','Tkk1','Tkk3','Pnbi','Picrh']
df = pd.DataFrame(columns=mycol)
# Assign values to the columns
df.time = time_ts
df.psi = PSI.T
df.Tts = Tts.T[:,0]
df.ErrTS = ErrY.T  # Error estimate for TS measurements as elaborated above
df.Tkk1 = Tkk1.T[:,0]
df.Tkk3 = Tkk3.T[:,0]
df.Pnbi = Pnbi.T[:,0]
df.Picrh = Picrh.T[:,0]

