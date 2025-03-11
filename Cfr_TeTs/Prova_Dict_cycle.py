#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:02:46 2025

@author: lsenni
"""

shots = [104522, 104525, 104575]

d = {
    'window_size' : 5,  # Numero punti per media mobile calcolo di ti e tf
    # 'tlim1' : 47,     # 47,     # limite inferiore selezione tempi - in secondi
    # 'tlim2' : 53,     # 53,     # limite superiore selezione tempi - in secondi
    'd2' : 8,       # Mezzo intervallo, in cm, sul quale vengono effettuate le medie (rho1, e rho 2)
    'np' : 2,       # numeri di punti di acquisizione che vanno aggiunti per determinare rho1 e rho2
    'fix' : 0.01,  
    }
 
vars = {
     'ti' : None, # Tempo iniziale finestra analisi
     'tf' : None, # Tempo finale foinestra analisi
     'tlim' : None,          # tempo medio al quale vengono calcolati i profili per plot di controllo
     'Ti' : None,  # Temp HRTS/max dell'istante iniziale
     'Tf' : None,  # Temp HRTS/max dell'istante finale
     'rho1' : None,    # quelli calcolati con porocedura automatica
     }

DB = {}

for shot in shots:
    d = d
    vars = vars
    DB[shot] = {'d':d, 'vars':vars}

