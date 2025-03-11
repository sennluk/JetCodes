#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:16:33 2024
@author: lsenni
Fig 02 per poster Alto Beta Orsitto ad EPS 2024 
Vengono riportati i dati del Beta Normalizzato vs la potenza NBI
nell'istante di massimo BetaN, per tre valori distinti di campo magnetico,
al variare di Pnbi
_EFTP: Confronto tra betaN calcolato con EFIT e copn EFTP
"""

from ppfeg import ppfs
import matplotlib.pyplot as plt
import numpy as np
import csv
# import mylib as my

plt.close('all') 

pulses17 =  [102305,102306, 102308, 102312,102433,102442] # Shots at Bt = 1.7 T: 102307, 102309, 102311, 102450,102451 STRAN3
pulses20 =  [102705, 102710,105318, 105320, 105321, 105322, 105323]  # Shots at Bt = 2.0 T - 105317,105319,
pulses24 =  [103114,103115,103116,103117,103122,103123] # Shots at Bt = 2.4 T


# shot = 103117  
tlim1 = 43    # limite inferiore selezione tempi - in secondi
tlim2 = 52 
delta = 3      # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.00      # 

cfr = 2 # n.cifre decimal approx valori calcolati

# shots = pulses14
bts = (17,20,24)

defit = {'tab%s' % bt:[] for bt in bts}
deftp = {'tab%s' % bt:[] for bt in bts}

for t,bt in enumerate(bts):
    shots = eval(('pulses'+str(bts[t])))
    print(shots)
    num = len(shots)
    print(num)
    
    defit['tab%s' % bt] = np.ones((num,4))
    defit[('pulses'+str(bts[t]))] = shots
    deftp['tab%s' % bt] = np.ones((num,4))
    deftp[('pulses'+str(bts[t]))] = shots
    
    print(t,bt)

    for i,shot in enumerate(shots):
        w = ppfs(shot) 
        
        pnbi_ = w.nbi.ptot   # Total NBI Power
        idt = (pnbi_.t >= tlim1) & (pnbi_.t <= tlim2) 
        tPnbi = pnbi_.t[idt]
        pnbi = pnbi_.v[0,idt]/10**6   
        
        # betan_ = w.efit.btnm      # Beta normalizaed - MHD
        # idt = (betan_.t >= tlim1) & (betan_.t <= tlim2)
        # tBetan = betan_.t[idt]
        # betan = betan_.v[0,idt]
        
        # betan2_ = w.eftp.btnm      # Beta normalizaed - MHD
        # idt = (betan2_.t >= tlim1) & (betan2_.t <= tlim2)
        # tBetan2 = betan2_.t[idt]
        # betan2 = betan2_.v[0,idt]
############# altri canali: btnd        
        betan_ = w.efit.btnd      # Beta normalizaed - MHD
        idt = (betan_.t >= tlim1) & (betan_.t <= tlim2)
        tBetan = betan_.t[idt]
        betan = betan_.v[0,idt]
        
        betan2_ = w.eftp.btnd      # Beta normalizaed - MHD
        idt = (betan2_.t >= tlim1) & (betan2_.t <= tlim2)
        tBetan2 = betan2_.t[idt]
        betan2 = betan2_.v[0,idt]
#########################################        
        maxb = max(betan)
        t_maxb = tBetan[np.where(betan == maxb)]
        
        maxb2 = max(betan2)
        t_maxb2 = tBetan2[np.where(betan2 == maxb2)]
        
        ind_time_Pnbi = np.argmin(abs(tPnbi-t_maxb))
        t_pnbi = tPnbi[ind_time_Pnbi]
        val_pnbi = pnbi[ind_time_Pnbi]
        
        ind_time_Pnbi2 = np.argmin(abs(tPnbi-t_maxb2))
        t_pnbi2 = tPnbi[ind_time_Pnbi2]
        val_pnbi2 = pnbi[ind_time_Pnbi2]
        
        print(i, shot)
        # Creo il dizionario con i dati
        defit['tab%s' % bt][i,:] = ["%.2f" %t_maxb, "%.2f" %maxb, "%.2f" %t_pnbi, "%.2f" %val_pnbi]
        deftp['tab%s' % bt][i,:] = ["%.2f" %t_maxb2, "%.2f" %maxb2, "%.2f" %t_pnbi2, "%.2f" %val_pnbi2]
        
        plt.figure(f'JPN{shot}')
        plt.plot(tBetan, betan, label = 'BetaN-Efit')
        plt.plot(tBetan2, betan2, label = 'BetaN-Eftp')
        plt.xlabel('time (sec)')
        plt.ylabel('Beta N')
        plt.legend()
        plt.savefig(f'V2_JPN{shot}_BetaN_EFIT_EFTP.pdf',dpi=300)
    
# pippo = tPnbi[np.where(tPnbi == time_maxb)] # altro metodo per fare la stessa cosa
# ##############################
linew = 0.5  # plot  lines dimension
fonts = 5.5    # size of the legend fonts
fs = 8       # size of the axes labels
fst = 6     # size of the ticks

fig, ax = plt.subplots(nrows=1, sharex=True, num='BetaN vs Pnbi')

ax.scatter(defit['tab17'][:,3], defit['tab17'][:,1], lw = linew, marker = 'v', color = 'blue', label='Bt = 1.7 T - EFIT')
ax.scatter(defit['tab20'][:,3], defit['tab20'][:,1], lw = linew, marker = '^', color = 'springgreen', label='Bt = 2.0 T - EFIT')
ax.scatter(defit['tab24'][:,3], defit['tab24'][:,1], lw = linew, marker = '>', color = 'm', label='Bt = 2.4 T - EFIT')

ax.scatter(deftp['tab17'][:,3], deftp['tab17'][:,1], lw = linew, marker = '.', color = 'blue', label='Bt = 1.7 T - EFTP')
ax.scatter(deftp['tab20'][:,3], deftp['tab20'][:,1], lw = linew, marker = ',', color = 'springgreen', label='Bt = 2.0 T - EFTP')
ax.scatter(deftp['tab24'][:,3], deftp['tab24'][:,1], lw = linew, marker = 'o', color = 'm', label='Bt = 2.4 T - EFTP')

ax.set_xlabel('Pnbi (MW)')
ax.set_ylabel('BetaN', fontsize= fs)
ax.yaxis.set_tick_params(labelsize=fst)
ax.legend(fontsize = fonts) #, loc="upper right"
ax.set_title(' BetaN vs Pnbi for $B_T$ = 1.7 - 2.0 - 2.4  T')

plt.savefig('V2_Fig_02_Beta_EFIT_E_EFTP_vs_Pnbi.pdf',dpi=300)

save = 0 # 0: non salva, 1:salva .csv , 2: salva .npy
filename = 'BetaN_e BetaN2_vs_Pnbi_V2'
# my.mysave(1,filename,d)


if save == 1:
            # Open a csv file for writing
    with open(f'{filename}_EFIT.csv', "w", newline="") as fp:
        # Create a writer object
        writer = csv.DictWriter(fp, fieldnames=defit.keys())
        # Write the header row
        writer.writeheader()
        # Write the data rows
        writer.writerow(defit)
        print('Done writing dict to a cs')
        
        with open(f'{filename}_EFTP.csv', "w", newline="") as fp:
            # Create a writer object
            writer = csv.DictWriter(fp, fieldnames=deftp.keys()) 
            # Write the header row
            writer.writeheader()
            # Write the data rows
            writer.writerow(deftp)
            print('Done writing dict to a cs')

