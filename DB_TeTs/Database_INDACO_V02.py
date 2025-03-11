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

# pulses_fi = [96994, 99869, 99801, 99802, 99870, 99950, 99970, 92415, 95986, 100793] 
pulses_prova = [95986, 100793] # prova: in uno esistono tutti i dati, nel 100793 manca il Michelson (ecm1)

# Decide on saving file options# Decide on saving file options
save = 0 # 0: non salva, 1:salva .csv , 2: salva .npy
filename = 'Database'

##################
# shot = 104522
pulses = pulses_prova #[99815] #,9917,99818]
ddas = ['hrts','ecm1','kk3','nbi','icrh']  # canali richiamati nel seguito


db = {}
for shot in pulses:
    
    w = ppfs(shot)
    rad = 3 # Position selected
    #########################
    # Slection of the data to be saved
    HRTS = w.hrts.te       # Te by HRTS ( Time, radial positions)
    Errts = w.hrts.dte    # Chann Errors HRTS
    Dens = w.hrts.ne       # Densisty Ne from HRTS
    ErrDens = w.hrts.dne   # Errors on density
    PSI = w.hrts.psi      # EFIT psi HRTS
    # KK1 = w.ecm1.prfl     # Te by ECE Michelson (Time, radial positions)
    KK3 = w.kk3.tprf    # Te by ECE Radiometer(Time, radial positions)
    NBI = w.nbi.ptot    # Total power of the NBI beam
    ICRH = w.icrh.ptot  # Total power of ICRH
  