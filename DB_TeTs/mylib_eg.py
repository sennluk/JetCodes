#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 12:29:13 2024
@author: lsenni- 25/01/2024 Aiuti e consigli da E.Gio.
Prova creazione libreria per database, 
dove metto un po' di funzioni da richiamare per snellire il tutto:
- min_time: calcola il minimo tempo, se esiste, altrimenti lo mette a zero
- get_at_radius: estrae i dati nella posizione e nell'intervallo temporale scelti
- max_time: calcola il massimo tempo, se esiste, altrimenti lo mette a 80     
- mydata: richiama le DDA e acquisisce i dati se presenti, altrimenti mette un None
- myslice: taglia l'intervallo temporale di interesse eseleziona la coordinata radiale di interesse
- myres: Resample dei dati sulla base del Thomson
- mysave: salva in .csv o in  .npy    

"""
import numpy as np
from ppfeg import ppfs, jetdata, V2d
from scipy.signal import resample

######################
######################

def min_time(v):
    return v.t.min() if v.ok() else 0

def max_time(v):
    return v.t.max() if v.ok() else 80

def get_at_radius(v, t1, t2, rad):
    if v.ok():
        idte = (v.t >= t1) & (v.t <= t2)
        TimeV = v.t[idte]
        TempV = v.v[:,idte]
        idre = np.argmin(abs(v.r - rad))
        RV = v.r[idre][np.newaxis]
        TV = TempV[idre,:][np.newaxis,:]
    else:
        TimeV = np.full((1,1),np.nan)
        TempV = np.full((1,1),np.nan)
        RV = np.full((1,1), np.nan)
        TV = np.full((1,1),np.nan)
    return TimeV, TempV, RV, TV


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
    
    pulses = pulses_fi
    # pulses = np.array(P05 + P06)
    return(pulses)

##################
######################
#####################
def mydata(shot):
    w = ppfs(shot)
    print('Processing pulse', shot)
    dda = dir(w)
    HRTS = jetdata(shot, 'HRTS','te')
    Errts = jetdata(shot, 'HRTS','dte')
    Dens = jetdata(shot, 'hrts','ne')
    ErrDens = jetdata(shot, 'hrts','dne')
    PSI = jetdata(shot, 'hrts','psi')
    KK1 = jetdata(shot, 'ecm1','prfl')
    KK3 = jetdata(shot, 'kk3','tprf')
    NBI = jetdata(shot, 'nbi','ptot')
    ICRH = jetdata(shot, 'icrh','ptot')
    if not ICRH.ok():
        ICRH = V2d(r=np.zeros(1), t=np.array([40.0, 80.0]), v=np.array([0.0,0.0]))
            
    return(w, HRTS, Errts, Dens, ErrDens, PSI, KK1, KK3, NBI, ICRH)
######################
######################          
# Selection of time interval and the radial position of the data 
   
def myslice(rad, HRTS, Errts, Dens, ErrDens, PSI, KK1, KK3, NBI, ICRH):
    # To choose the time interval for each shot
    t1 = np.max([min_time(v) for v in [HRTS, KK1, KK3, NBI, ICRH]])
    t2 = np.min([max_time(v) for v in [HRTS, KK1, KK3, NBI, ICRH]])
    
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
    ErrY[ErrY > 100] = 0       # check on error values --> >1000 wrong value 
    
    ErrN = ErrNe[0,:]/2    # Errors on Ne divìded by 2 for the +- bar
    ErrN[ErrN > 1e20] = 0    # Check on error values --> >10^20 wrong value 
    
    # if KK1.ok():
        # idte = (KK1.t >= t1) & (KK1.t <= t2)
        # TimeKK1 = KK1.t[idte]
        # TempKK1 = KK1.v[:,idte]
        # idre = np.argmin(abs(KK1.r - rad))
        # Rkk1 = KK1.r[idre][np.newaxis]
        # Tkk1 = TempKK1[idre,:][np.newaxis,:]
    # else:
        # TimeKK1 = np.nan
        # TempKK1 = np.nan
        # Rkk1 = np.nan
        # Tkk1 = np.nan

    TimeKK1, TempKK1, Rkk1, Tkk1 = get_at_radius(KK1, t1, t2, rad)

    #KK1 = KK1.trange([t1, t2])
    #TimeKK1, Tkk1 = KK1.get(r=rad)
    
    print(Rkk1, Tkk1)
    
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
    return(TimeTS, Tts, ErrY, PSI, Ne, ErrN, TimeKK1, Tkk1, TimeKK3, Tkk3, TimeNBI, Pnbi, TimeICRH, Picrh)
        
###########################
##########################

 # Resample time axixs based on the slowest acq rate
def myres(TimeTS, Tkk1, TimeKK1, Tkk3, TimeKK3, Pnbi, TimeNBI, Picrh, TimeICRH):
     dim = TimeTS.size
     
     Time = TimeKK1
     v,t = resample(Tkk1, dim, t=Time, axis=1)
     Tkk1 = v
     TimeKK1 = t
     
     Time = TimeKK3
     v,t = resample(Tkk3, dim, t=Time, axis=1)
     Tkk3 = v
     TimeKK3 = t
     
     Time = TimeNBI
     v,t = resample(Pnbi, dim, t=Time, axis=1)
     Pnbi = v
     
     Time = TimeICRH
     v,t = resample(Picrh, dim, t=Time, axis=1)
     Picrh = v
     return(TimeKK1, Tkk1, TimeKK3, Tkk3, Pnbi, Picrh) 
 
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
