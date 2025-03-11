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

V03: utilizzo le functions che metto nella mylib.py, contenuta all'interno della
stessa cartella:
    mydata: richiama le DDA e acquisisce i dati se presenti, altrimenti mette un None
    myslice: taglia l'intervallo temporale di interesse eseleziona la coordinata radiale di interesse
    myres: Resample dei dati sulla base del Thomson
    mysave: salva in .csv o in  .npy    
"""
import pandas as pd
import mylib as my

# saving file options
save = 2 # 0: non salva, 1:salva .csv , 2: salva .npy
filename = 'DB_Fish'
rad = 3 # Position selected (m)

##################
# shot = 104522
# pulses = pulses_prova #[99815] #,9917,99818]
pulses = my.mypulses()   # Shot selection from mylib.py

chans = ['hrts','ecm1','kk3','nbi','icrh']  # canali richiamati nel seguito

db = {}
for shot in pulses:
    w, HRTS, Errts, Dens, ErrDens, PSI, KK1, KK3, NBI, ICRH = my.mydata(shot)
    TimeTS, Tts, ErrY, PSI, Ne, ErrN, TimeKK1, Tkk1, TimeKK3, Tkk3, TimeNBI, Pnbi, TimeICRH, Picrh = my.myslice(rad, HRTS, Errts, Dens, ErrDens, PSI, KK1, KK3, NBI, ICRH)
    TimeKK1, Tkk1, TimeKK3, Tkk3, Pnbi, Picrh = my.myres(TimeTS, Tkk1, TimeKK1, Tkk3, TimeKK3, Pnbi, TimeNBI, Picrh, TimeICRH)
       
    ##########################
    # Dataframe creation
    
    # Define the columns names
    mycol = ['time','psi','Tts','ErrTS','Ne','ErrNe','Tkk1','Tkk3','Pnbi','Picrh']
    df = pd.DataFrame(columns=mycol)
    # Assign values to the columns
    df.time = TimeTS           # Time in sec
    df.psi = PSI.T             # Psi
    df.Tts = Tts.T[:,0]/1000   # TeTS in KeV
    df.ErrTS = ErrY.T          # Error (KeV) estimate for TS measurements as elaborated above
    df.Ne = Ne.T[:,0]          # Density in m-3
    df.ErrNe = ErrN.T          # Error on density
    df.Tkk1 = Tkk1.T[:,0]/1000 # Te KK1 in keV   
    df.Tkk3 = Tkk3.T[:,0]
    df.Pnbi = Pnbi.T[:,0]
    df.Picrh = Picrh.T[:,0]

    # Dictionary voice of the current pulse/shot
    db[shot] = df
    
  
TotShots = len(pulses)
TotCol = len(df.columns)
print('Total number of processed Shots = '+str(TotShots))
print('Total number of columns = '+str(TotCol))
Legend = df.columns
print('legend:', Legend)


#####################################################
# Save Dictionary to a .csv file
my.mysave(save,filename,db)