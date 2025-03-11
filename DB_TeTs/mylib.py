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
- myslice: taglia l'intervallo temporale di interesse e seleziona la coordinata radiale di interesse
- myres: Resample dei dati sulla base del Thomson
- mysave: salva in .csv o in  .npy    

"""
import numpy as np
from ppfeg import ppfs, jetdata, V2d
from scipy.signal import resample

###########################
# Seleziono gli impulsi di interesse

def mypulses():
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
    pulses_fi = [96994, 99869, 99801, 99802, 99870, 99950, 99970, 92415, 95986, 100793] # Identified to study ishbones related to discr
    
    pulses = [95986,100793] #pulses_fi [95986]
    # pulses = np.array(P05 + P06)
    return pulses

##################
######################
#####################
# Richiama le DDA e acquisisce i dati se presenti, altrimenti mette dei valori non rilevanti
# Al fine di non bloccare l'esecuzione bisogna costruire un asse dei temi verosimile--> lo baso su quello dell'hrts
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
        kk1 = V2d(r=np.zeros_like(hrts.r), t = np.arange(40,80,(80-40)/N), v=np.zeros_like(hrts.v))  # metto tutti zero mantenendo le stesse dimensione dell'hrts
    kk3 = jetdata(shot, 'kk3', 'tprf')
    nbi = jetdata(shot, 'nbi', 'ptot')
    icrh = jetdata(shot, 'icrh','ptot')
    if not icrh.ok():                     # Se il dato esiste/è ok
        print('icrh channel not ok')
        icrh = V2d(r=np.zeros_like(hrts.r), t = np.arange(40,80,(80-40)/N), v=np.zeros_like(hrts.v))
            
    return hrts, errTs, dens, errDens, psi, kk1, kk3, nbi, icrh

######################          
# Selection of time interval and the radial position of the data 
# --> utilizza le funzioni min_time, max_time, e get_at_radius
# can: canale, DDA

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

def myslice(rad, hrts, errTs, dens, errDens, psi, kk1, kk3, nbi, icrh):
    # To choose the time interval for each shot
    t1 = np.max([min_time(dd) for dd in [hrts, kk1, kk3, nbi, icrh]])
    t2 = np.min([max_time(dd) for dd in [hrts, kk1, kk3, nbi, icrh]])
    if t2<t1:
        print('Wrong time axis selection')
    timeTs, tempTs = get_at_radius(hrts, t1, t2, rad)
    _, errTs = get_at_radius(errTs, t1, t2, rad)
    _, dens = get_at_radius(dens, t1, t2, rad)
    _, errDens = get_at_radius(errDens, t1, t2, rad)
    _, psiTs = get_at_radius(psi, t1, t2, rad)
    
    errTs = errTs[0,:]/2000      # Errors on Tts divìded by 2 for the +- bar, and /1000 per keV
    errTs[errTs > 100] = 0       # Check on error values --> >1000 wrong value 
    errDens = errDens[0,:]/2     # Errors on Ne divìded by 2 for the +- bar
    errDens[errDens > 1e20] = 0  # Check on error values --> >10^20 wrong value 
    
    timeKk1, tempKk1 = get_at_radius(kk1, t1, t2, rad)
    timeKk3, tempKk3 = get_at_radius(kk3, t1, t2, rad)
    timeNbi, pNbi = get_at_radius(nbi, t1, t2, rad)
    timeIcrh, pIcrh = get_at_radius(icrh, t1, t2, rad)
    
    return timeTs, tempTs, errTs, dens, errDens, psiTs, timeKk1, tempKk1, timeKk3, tempKk3, timeNbi, pNbi, timeIcrh, pIcrh
        
###########################
##########################

 # Resample time axis based on the slowest acq rate (HRTS in this case)
def myres(timeTs, timeKk1, tempKk1, timeKk3, tempKk3, timeNbi, pNbi, timeIcrh, pIcrh):
     dim = timeTs.size
     
     time = timeKk1
     v,t = resample(tempKk1, dim,  t=time, axis=1) # , domain='time'
     timeKk1 = t
     tempKk1 = v
         
     time = timeKk3
     v,t = resample(tempKk3, dim, t=time, axis=1)
     timeKk3 = t
     tempKk3 = v
     
     time = timeNbi
     v,t = resample(pNbi, dim, t=time, axis=1)
     pNbi = v
     
     time = timeIcrh
     v,t = resample(pIcrh, dim, t=time, axis=1)
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
