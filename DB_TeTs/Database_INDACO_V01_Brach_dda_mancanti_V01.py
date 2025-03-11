#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 14:12:23 2023
@author: lsenni
Creazione Database per INDACO
NOTA: Utilizzo il modulo di Edmondo Giovannozzi--> va caricato prima dello script direttamente da shell

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
import mylib as my

# pulses_fi = [96994, 99869, 99801, 99802, 99870, 99950, 99970, 92415, 95986, 100793] 
pulses_prova = [95986, 100793] # prova: in uno esistono tutti i dati, nel 100793 manca il Michelson (ecm1)

# Decide on saving file options# Decide on saving file options
save = 0 # 0: non salva, 1:salva .csv , 2: salva .npy
filename = 'Database'

##################
# shot = 104522
pulses = pulses_prova #[99815] #,9917,99818]
chans = ['hrts','ecm1','kk3','nbi','icrh']  # canali richiamati nel seguito

db = {}
for shot in pulses:
    w, HRTS, Errts, Dens, ErrDens, PSI, KK1, KK3, NBI, ICRH = my.mydata(shot)
    