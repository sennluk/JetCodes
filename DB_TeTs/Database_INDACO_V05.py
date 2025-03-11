#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 14:12:23 2023
@author: lsenni
Creazione Database per INDACO
NOTA: Utilizzo il modulo di Edmondo Giovannozzi--> va caricato prima dello script direttamente da shell

Tutti i dati vengono raccolti in un Dizionario, con ugual dimensioni 
corrispondenti al numero di punti temporali del TS, la diagnostica più lenta.

V05:  utilizzo le classsi, ma Nota Bene: i nomi con la maiuscola inizale sono le classi, 
se lo stesso nome ha la minuscola iniziale è un oggetto/metodo. 
inserisco le coordinate di equilibrio --> funzione my_flush di egio
"""
import pandas as pd
import mylib_B_Classi as my

# saving file options
save = 0 # 0: non salva, 1:salva .csv , 2: salva .npy
filename = 'DB_Fish'
rad = 3 # Radial position selected (m)
psi = 1.2 # Selected position in eq coordinate

##################
# shot = 104522
# pulses = pulses_prova #[99815] #,9917,99818]
pulses = my.mypulses()   # Shot selection from mylib.py

chans = ['hrts','ecm1','kk3','nbi','icrh']  # canali richiamati nel seguito

db = {}
for shot in pulses:
    dda4db = my.mydata(shot)
    data4db = my.myslice(rad, dda4db)
    timeKk1, tempKk1, timeKk3, tempKk3, pNbi, pIcrh = my.myres(data4db)
       
    ##########################
    # Dataframe creation
    
    # Define the columns names
    mycol = ['time','psiTs','Tts','ErrTS','Ne','ErrNe','Tkk1','Tkk3','Pnbi','Picrh']
    df = pd.DataFrame(columns=mycol)
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
my.mysave(save,filename,db)