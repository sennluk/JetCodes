#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 12:29:13 2024
@author: lsenni- 25/01/2024 Aiuti e consigli da E.Gio.
Prova creazione libreria per database, 
dove metto un po' di funzioni da richiamare per snellire il tutto:
- min_time: calcola il minimo tempo, se esiste, altrimenti lo mette a zero
- max_time: calcola il massimo tempo, se esiste, altrimenti lo mette a 80     
- get_at_radius: estrae i dati nella posizione e nell'intervallo temporale scelti
                 controllando che il dato sia ok (esista)
- mypulses: selezione degli impulsi di interesse                
- mydata: richiama le DDA e acquisisce i dati se presenti, altrimenti mette dei valori non rilevanti
- myslice: taglia l'intervallo temporale di interesse eseleziona la coordinata radiale di interesse
- myres: Resample dei dati sulla base del Thomson
- mysave: salva in .csv o in  .npy    

"""
import numpy as np
import ppfeg 
from ppfeg import ppfs, jetdata, V2d
from scipy.signal import resample
from dataclasses import dataclass

###########################
# Seleziono gli impulsi di interesse

def mypulses():
   
    pulses = [95986,100793] #pulses_fi [95986]
    # pulses = np.array(P05 + P06)
    return pulses

##################
######################
#####################
# Richiama le DDA e acquisisce i dati se presenti, altrimenti mette dei valori non rilevanti
# Al fine di non bloccare l'esecuzione bisogna costruire un asse dei temi verosimile--> lo baso su quello dell'hrts
@dataclass
class Dda4db:
    hrts: ppfeg.V2d
    errTs: ppfeg.V2d
    dens: ppfeg.V2d 
    errDens: ppfeg.V2d
    psi: ppfeg.V2d
    kk1: ppfeg.V2d
    kk3: ppfeg.V2d
    nbi: ppfeg.V2d
    icrh: ppfeg.V2d

def mydata(shot):
    print('Processing pulse', shot)
    hrts = jetdata(shot, 'hrts', 'te')
    errTs = jetdata(shot, 'hrts', 'dte')
    dens = jetdata(shot, 'hrts', 'ne')
    errDens = jetdata(shot, 'hrts', 'dne')
    psi = jetdata(shot, 'hrts', 'psi')
    kk1 = jetdata(shot, 'ecm1', 'prfl')
    N = hrts.t.size
    if not kk1.ok():                     # Se il dato esiste/è ok
        print('kk1 channel not ok')
        kk1 = V2d(r=np.zeros_like(hrts.r), t = np.linspace(40,80,N) , v=np.zeros_like(hrts.v))  # metto tutti zero mantenendo le stesse dimensione dell'hrts
         # t = np.arange(40,80,(80-40)/(N-1)) opp. t = np.linspace(40,80,N) 
    kk3 = jetdata(shot, 'kk3', 'tprf')
    nbi = jetdata(shot, 'nbi', 'ptot')
    icrh = jetdata(shot, 'icrh','ptot')
    if not icrh.ok():                     # Se il dato esiste/è ok
        print('icrh channel not ok')
        icrh = V2d(r=np.zeros_like(hrts.r), t = np.arange(40,80,(80-40)/N), v=np.zeros_like(hrts.v))
            
    return Dda4db(hrts=hrts, errTs=errTs, dens=dens, errDens=errDens, psi=psi, kk1=kk1, kk3=kk3, nbi=nbi, icrh=icrh)

######################          
# Selection of time interval and the radial position of the data 
# --> utilizza le funzioni min_time, max_time, e get_at_radius
# can: canale, DDA
@dataclass
class Data4db:
    timeTs: np.ndarray
    tempTs: np.ndarray
    errTs: np.ndarray
    dens: np.ndarray
    errDens: np.ndarray
    psiTs: np.ndarray
    timeKk1: np.ndarray
    tempKk1: np.ndarray
    timeKk3: np.ndarray
    tempKk3: np.ndarray
    timeNbi: np.ndarray
    pNbi: np.ndarray
    timeIcrh: np.ndarray
    pIcrh: np.ndarray    
    
def min_time(can):
    return can.t.min() if can.ok() else 0

def max_time(can):
    return can.t.max() if can.ok() else 80

def get_at_radius(can, t1, t2, rad):
    idte = (can.t >= t1) & (can.t <= t2)
    can_t = can.t[idte]
    can_v = can.v[:,idte]           # can.v data
    idre = np.argmin(abs(can.r - rad)) 
    can_r = can.r[idre][np.newaxis]
    can_rt = can_v[idre,:][np.newaxis,:] # dato calcolato nella posizione e nelll'intervallo di tempo selezionato
    return can_t, can_rt  

def myslice(rad, p:Dda4db): #hrts, errTs, dens, errDens, psi, kk1, kk3, nbi, icrh):
    # To choose the time interval for each shot
    t1 = np.max([min_time(dd) for dd in [p.hrts, p.kk1, p.kk3, p.nbi, p.icrh]])
    t2 = np.min([max_time(dd) for dd in [p.hrts, p.kk1, p.kk3, p.nbi, p.icrh]])
    if t2<t1:
        print('Wrong time axis selection')
    timeTs, tempTs = get_at_radius(p.hrts, t1, t2, rad)
    _, errTs = get_at_radius(p.errTs, t1, t2, rad)
    _, dens = get_at_radius(p.dens, t1, t2, rad)
    _, errDens = get_at_radius(p.errDens, t1, t2, rad)
    _, psiTs = get_at_radius(p.psi, t1, t2, rad)
    
    errTs = errTs[0,:]/2000      # Errors on Tts divìded by 2 for the +- bar, and /1000 per keV
    errTs[errTs > 100] = 0       # Check on error values --> >1000 wrong value 
    errDens = errDens[0,:]/2     # Errors on Ne divìded by 2 for the +- bar
    errDens[errDens > 1e20] = 0  # Check on error values --> >10^20 wrong value 
    
    timeKk1, tempKk1 = get_at_radius(p.kk1, t1, t2, rad)
    timeKk3, tempKk3 = get_at_radius(p.kk3, t1, t2, rad)
    timeNbi, pNbi = get_at_radius(p.nbi, t1, t2, rad)
    timeIcrh, pIcrh = get_at_radius(p.icrh, t1, t2, rad)
    
    return Data4db(timeTs=timeTs, tempTs=tempTs, errTs=errTs, dens=dens, errDens=errDens, psiTs=psiTs, timeKk1=timeKk1, tempKk1=tempKk1, timeKk3=timeKk3, tempKk3=tempKk3, timeNbi=timeNbi, pNbi=pNbi, timeIcrh=timeIcrh, pIcrh=pIcrh)
        
###########################
##########################

 # Resample time axis based on the slowest acq rate (HRTS in this case)
def myres(p:Data4db): #timeTs, timeKk1, tempKk1, timeKk3, tempKk3, timeNbi, pNbi, timeIcrh, pIcrh):
     dim = p.timeTs.size
     
     time = p.timeKk1
     v,t = resample(p.tempKk1, dim,  t=time, axis=1) # , domain='time'
     timeKk1 = t
     tempKk1 = v
         
     time = p.timeKk3
     v,t = resample(p.tempKk3, dim, t=time, axis=1)
     timeKk3 = t
     tempKk3 = v
     
     time = p.timeNbi
     v,t = resample(p.pNbi, dim, t=time, axis=1)
     pNbi = v
     
     time = p.timeIcrh
     v,t = resample(p.pIcrh, dim, t=time, axis=1)
     pIcrh = v
     
     return timeKk1, tempKk1, timeKk3, tempKk3, pNbi, pIcrh 
 
###########################
##########################
# Save db with different extesion
def mysave(save,filename,db):
    # Save Dictionary to a .csv file
    if save == 1:
        import csv
                # Open a csv file for writing
        with open(f'{filename}.csv', "w", newline="") as fp:
            # Create a writer object
            writer = csv.DictWriter(fp, fieldnames=db.keys())
            
            # Write the header row
            writer.writeheader()
            
            # Write the data rows
            writer.writerow(db)
            print('Done writing dict to a csv file')
            
    #####################################################
    # Save Dictionary using NumPy
    if save ==2:
        np.save(f'{filename}.npy', db) 
    return()

# Load
# read_dictionary = np.load('my_file.npy',allow_pickle='TRUE').item()
# print(read_dictionary['hello']) # displays "world" 
