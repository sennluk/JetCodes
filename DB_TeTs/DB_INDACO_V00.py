#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 14:12:23 2023
@author: lsenni
Creazione Database per INDACO:
    Tutti i dati vengono raccolti in un Dizionario, con ugual dimensioni 
    corrispondenti al numero di punti temporali del TS, la diagnostica più lenta.
    A partire dalla versione  Database_Indaco_V05, provo a ricostruire il Db
    utilizzando la coordinata di equilibrio psi. Per tutte le osservazioni 
    vedere i coommenti alle versioni precedenti
    
NOTA: Utilizzo il modulo di Edmondo Giovannozzi--> va caricato prima dello script direttamente da shell
Utilizzo il file di libreria contenente tutte le funzioni lib_db.py

DB_Indaco_V00: 
Inserisco altre colonne e gli equilibri oper le diagnostiche che non hanno i canali gia fatti
"""
import pandas as pd
import lib_db as lib

# saving file options
save = 0 # 0: non salva, 1:salva .csv , 2: salva .npy
filename = 'Jet_DataBase'

rad = 3 # Radial position selected (m)
psi_1 = 0.07
psi_2 = 0.12 # Selected position in eq coordinate

##################
# shot = 104522
# pulses = pulses_prova #[99815] #,9917,99818]
pulses = lib.mypulses()   # Shot selection from mylib.py

chans = ['hrts','ecm1','kk3','nbi','icrh']  # canali richiamati nel seguito

db = {}
for shot in pulses:
    dda4db = lib.mydata(shot)
    data4db = lib.myslice(rad, dda4db)
    timeKk1, tempKk1, timeKk3, tempKk3, pNbi, pIcrh = lib.myres(data4db)
       
    ##########################
    # Dataframe creation
    
    # Define the columns names
    col_names = ['time','psiTs','Tts','ErrTS','Ne','ErrNe','Tkk1','Tkk3','Pnbi','Picrh']
    df = pd.DataFrame(columns=col_names)
    # Assign values to the columns
    df.time = data4db.timeTs                # Time in sec
    df.psiTs = data4db.psiTs.T                # Psi
    df.Tts = data4db.tempTs.T[:,0]/1000     # TeTS in KeV
    df.ErrTs = data4db.errTs.T              # Error (KeV) estimate for TS measurements as elaborated above
    df.Ne = data4db.dens.T[:,0]             # Density in m-3
    df.ErrNe = data4db.errDens.T            # Error on density
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
lib.mysave(save,filename,db)