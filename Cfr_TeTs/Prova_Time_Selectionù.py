#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 09:52:06 2025

@author: lsenni
"""
import numpy as np
import my_flush
from ppfeg import ppfs
from scipy import interpolate
import my_manage_file as mym
##################################################
# Definisco il dizionario delle variabili

# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot
JPN = 104990 #104990 # 104522 # 104525 #104990  #96994
# DTE3 shot list: 104990,991,994,995,999   
# 104520,521,522,523,524,526

timestr = 5# time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation

d = {
    'shot' : JPN,     # 99950, 99971, 99970
    'Tref' : 1,  # keV) Reference temperature: time instants with Te values above Tref are considered
    'window_size' : 10,  # Numero punti per media mobile calcolo di ti e tf
    'min_increase': 0.0, # minimum increase in the tiome window to be considered for the ti 
    'd2' : 0.08,       # Mezzo intervallo, in cm, sul quale vengono effettuate le medie (rho1, e rho 2)
    'np' : 2,       # numeri di punti di acquisizione che vanno aggiunti per determinare rho1 e rho2
    'fix' : 0.01,    # intervallo aggiuntivo da considerare (in psi e rho) nella differenza tra le due
    'intPsi' : 0.1,    # intervallo in Psi e Rho su cui mediare 
    'delta' : 3,      # tempo in più di tlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
    'rad' : 3.05,     # raggio al quale viene fatta l'analisi (in metri)
    'psi1' : 0.003,  # 0.001 Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
    'psi2' : 0.027,     # 0.027 #0.02
    'eP' : 0.02,      # Relative error assigned to the Ece data
    'win_len' : 15,   # window lenght: number of points for the smooth with Savitsky-golay filter
    'deg_pol' : 3,    # grado del polinomio usato per lo smooting
    'savefigs' : 0,   # 1--> save plots, 0 don't save.
    'mypath' : f'/home/lsenni/Python_LS/Cfr_TeTs/{timestr}_JPN_{JPN}_Plots/'  # folder to save plot
    }    
#######

shot = d['shot']
# timestr = time.strftime("%Y%m%d-%H%M%S")

w = ppfs(shot)   

if d['savefigs'] == 1: 
    mym.create_folder(d)
    mym.save_param(d)
    # mym.save_print_output_to_file(d,vars)
else:
    print('You are not saving the plots')
 
vars = {
    'w' : w,                           # canale dati del ppfeg con tutti i dati relativi allo sparo in oggetto
    'ti' : None, # Tempo iniziale finestra analisi
    'tf' : None, # Tempo finale foinestra analisi
    'tlim' : None,          # tempo medio al quale vengono calcolati i profili per plot di controllo
    'Ti' : None,  # Temp HRTS/max dell'istante iniziale
    'Tf' : None,  # Temp HRTS/max dell'istante finale
    'psi1' : d['psi1'],    # inizializzo i valori di psi1,2 e rho1,2
    'psi2' : d['psi2'],    # questi possono essere usati se non si vuole usare
    'rho1' : None,    # quelli calcolati con porocedura automatica
    'rho2' : None,
    'rhodown' : None,     # rho
    'rhoup' : None,
    'tTs' : w.hrts.te,                 # Chann Te HRTS  
    'errTs' : w.hrts.dte,              # Chann Errors HRTS.dte
    'tEce' : w.ecm1.prfl,              # Chann Te Ece-Kk1
    'errEce' : (w.ecm1.prfl)*d['eP'],  # Error associated to the Kk1 meas: ab 2%
    'psiTs' : w.hrts.psi,              # Chann PSI HRTS 
    'psiKk1' : None,
    'psiTscalc' : None,
    'rhoKk1' : None,
    'time_ts' : None,
    'time_ece' : None,
    'rhoTs' : None,
    'rhoEce' : None, 
    'temp_tsM' : None, 
    'err_tsM': None,
    'timeTs2' : None,  # TS time instants in the ti-tf window
    'timeEce22' : None, # ECE time instants closest to the TS ones in the ti-tf window     '
    'ranges' : None,  # RHO ranges as calculated by the automatic procedure - TS interpolated in t1<t<t2
    'xm': None, 
    'err_xm' : None, 
    'temp_eceM' : None, 
    'err_eceM' : None, 
    'xm12': None, 
    'err_xm12' : None, 
    'temp_eceM12' : None, 
    'err_eceM12' : None, 
    'temp_tsM_rho' : None,
    'err_tsM_rho' : None,
    'temp_eceM_rho' : None,
    'err_eceM_rho' : None} 


################

# shot = d['shot']
w = vars['w']
Tref = d['Tref']
window_size = d['window_size']
min_increase = d['min_increase']
#####################################
tmax = w.hrtx.tmax/1000     # max Te max hrts (hrtx channel) in keV
data = tmax.v[0,:]
# prendo come tlim l'istante in cui la Tmax è massima
valid_indices = np.where(tmax.t <= 60)[0]
tlim = tmax.t[np.nanargmax(tmax.v[0,valid_indices])] 
ind_tlim = np.nanargmax(tmax.v[0,valid_indices])
TMAX = data[ind_tlim]

def moving_average_np(data, window_size):
    """Calcola la media mobile semplice usando NumPy."""
    return np.convolve(data, np.ones(window_size) / window_size, mode='valid')

def find_rising_point_np(data, window_size, threshold, min_increase):
    """Trova il primo punto in cui i dati iniziano a salire."""
    # Calcola la media mobile
    smoothed_data = moving_average_np(data, window_size)
    
    # Trova il primo punto in cui il valore aumenta
    for i in range(1, len(smoothed_data)):
        if smoothed_data[i] > smoothed_data[i - 1]+ min_increase and smoothed_data[i] > threshold:
            return i + window_size - 1  # Indice relativo all'array originale
    
    return None  # Nessun punto di crescita trovato

def find_return_point_np(data, start_index, target_value):
    """Trova il punto in cui i dati tornano al valore target dopo l'indice dato."""
    for i in range(start_index + 1, len(data)):
        if data[i] <= target_value and i > i+20 and i<len(data)-50 : #+20:  # Modifica qui se vuoi confrontare con <= o ==
            return i
        data_subset = data[ind_tlim + 1:]  # Prendi solo la parte dopo start_index
        closest_index = np.argmin(np.abs(data_subset - target_value))  # Trova il valore più vicino
        closest_index += start_index + 1  # Riporta l'indice nel riferimento originale
   
    return closest_index  # Restituisce l'indice del valore più vicino
    # Se non troviamo un punto di ritorno, cerchiamo il primo valore diverso da zero
    # for i in range(start_index + 1, len(data)):
    #     if data[i] != 0:
    #         return i  # Restituisce il primo valore diverso da zero
    # return None  # Nessun punto di ritorno trovato

# def find_last_valid_point(data, start_index, threshold):
#     """Trova l'ultimo valore di 'data' (dopo 'start_index') che sia diverso da zero
#        e maggiore di 'threshold', scorrendo da destra verso sinistra."""
    
#     for i in range(len(data) - 1, start_index, -1):  # Scansiona all'indietro
#         if data[i] > threshold and data[i] != 0:
#             return i  # Restituisce l'indice dell'ultimo valore valido
    
#     return None  # Nessun valore trovato
#####################
rising_index = find_rising_point_np(data, window_size, Tref, min_increase)

if rising_index is not None:
    rising_value = data[rising_index]
    print(f"I dati iniziano a salire all'indice {rising_index}, valore: {rising_value}")
    
    # Trova il punto di ritorno
    return_index = find_return_point_np(data, ind_tlim, rising_value-0.2)
    Tfin = data[return_index]
    if return_index is not None:
        print(f"L'inice finale è: {return_index}, corrispondente ad una valore: {Tfin}")
    else:
        print(f"I dati non tornano più al valore {rising_value}.")
    #     # RV = rising_value+3
    #     return_index = closest_index# find_return_point_np(data, ind_tlim, RV)  # Cerca il primo valore ≠ 0 dopo tlim

    #     if return_index is not None:
    #         print(f"Ho trovato un valore diverso da zero all'indice {return_index}")
    #     else:
    #         print("Nessun valore diverso da zero trovato dopo tlim.")

        # print(f"I dati non tornano più al valore {rising_value}.")
        # pluto = find_return_point_np(data, ind_tlim, 0)
        # return_index = pluto-1
        
        # target = np.min(abs(data[0,ind_tlim]-rising_value))
        # return_index = find_return_point_np(data, ind_tlim, target) 
        # return_index = np.where(tmax.t[rising_index] + 15)
        # return_index = find_return_point_np(data, rising_index, rising_value+2)
        # print(f"I dati non tornano più al valore {rising_value}.")
        # print("Prendo un valore a + 15 secondi")
else:
    print("Non ci sono punti in cui i dati iniziano a salire.")

#####################
Ti = rising_value          # Temp istante inizale
Tf = data[return_index]    # Temp istante finale
   
ti = tmax.t[rising_index]  # Istante iniziale
tf= tmax.t[return_index]   # Istante finale

tlim1 = vars['ti']
tlim2 = vars['tf']

vars['ti'] = ti
vars['tf'] = tf
vars['Ti'] = Ti
vars['Tf'] = Tf
vars['tlim'] = tlim

print('JPN = ', d['shot'])
print("Automated tlim1 =", ti)
print("Automated tlim2 =", tf)
print('t lim1 = ', tlim1)
print('t lim2 = ', tlim2)
print('Instant of the max Tmax=', tlim)