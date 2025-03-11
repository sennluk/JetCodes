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
V04: dopo l'incontro con egio cambio il nome alle variabili, e cambio il file mylib:
    necessario qualche accorgimento per la ricostruzione dei dati quando mancanti, in particolare per l'asse dei tempi
    e per non bloccare le successive elaborazioni
    
"""
import pandas as pd
import mylib as my

# saving file options
save = 0 # 0: non salva, 1:salva .csv , 2: salva .npy
filename = 'DB_Fish'
rad = 3 # Position selected (m)

##################
# shot = 104522
# pulses = pulses_prova #[99815] #,9917,99818]
pulses = my.mypulses()   # Shot selection from mylib.py

chans = ['hrts','ecm1','kk3','nbi','icrh']  # canali richiamati nel seguito

db = {}
for shot in pulses:
    hrts, errTs, dens, errDens, psi, kk1, kk3, nbi, icrh = my.mydata(shot)
    timeTs, tempTs, errTs, dens, errDens, psiTs, timeKk1, tempKk1, timeKk3, tempKk3, timeNbi, pNbi, timeIcrh, pIcrh = my.myslice(rad, hrts, errTs, dens, errDens, psi, kk1, kk3, nbi, icrh)
    timeKk1, tempKk1, timeKk3, tempKk3, pNbi, pIcrh = my.myres(timeTs, timeKk1, tempKk1, timeKk3, tempKk3, timeNbi, pNbi, timeIcrh, pIcrh)
       
    ##########################
    # Dataframe creation
    
    # Define the columns names
    mycol = ['time','psi','Tts','ErrTS','Ne','ErrNe','Tkk1','Tkk3','Pnbi','Picrh']
    df = pd.DataFrame(columns=mycol)
    # Assign values to the columns
    df.time = timeTs                # Time in sec
    df.psi = psiTs.T                # Psi
    df.Tts = tempTs.T[:,0]/1000     # TeTS in KeV
    df.ErrTs = errTs.T              # Error (KeV) estimate for TS measurements as elaborated above
    df.Ne = dens.T[:,0]             # Density in m-3
    df.ErrNe = errDens.T            # Error on density
    df.Tkk1 = tempKk1.T[:,0]/1000   # Te KK1 in keV   
    df.Tkk3 = tempKk3.T[:,0]
    df.Pnbi = pNbi.T[:,0]
    df.Picrh = pIcrh.T[:,0]

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