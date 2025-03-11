#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:17:34 2024

@author: lsenni
"""


import numpy as np
import ppfeg
import my_flush
from ppfeg import ppfs
import matplotlib.pyplot as plt
from scipy import signal,interpolate


plt.close('all')  
# save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e la posizione radiale di interesse

shot = 104522 # 99950, 99971
tlim1 = 47.5    # limite inferiore selezione tempi - in secondi
tlim2 = 53 
delta = 1     # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.0     # raggio al quale viene fatta l'analisi (in metri)
psi1 = 0.01   # Intervallo in psi su cui mediare: 0.01-0.1 per vecchia configurazione
psi2 = 0.02   # 
eP = 0.02  # Relative error assigned to the Ece data 
tlim = (tlim1+tlim2)/2

w = ppfs(shot)
############################################
# Calcolo rho poloidale normalizzato - a titolo di esempio metto la funzione qua sotto
# che poi utilizzerò nel ciclo
polflu = w.efit.ftor  # poloidal flux - Segno negativo per il verso del campo
# Funzione che peremette il calcolo della funzione da interpolare quando si fornirà
# in input la funzione sulla baase della quale fare l'interpolazione - al tempo di indice 50
fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,50]/w.efit.ftor.v[-1,50])
# Flusso interpolato sui dati della PSI HRTS
fl_int_hrts = fl_int(w.hrts.psi.v[:,50]) 
rho = np.sqrt(fl_int_hrts)

# for i,time in enumerate(time_ece):
#      psi = psiKk1[i]
#      time = time_ece[i]
#      fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i])
#      fl_int_hrts = fl_int(psiKk1[i])
#      rho_ece[i] = np.sqrt(fl_int_hrts)

plt.figure()
w.efit.ftor.plot()
plt.xlabel('time (sec)')
plt.ylabel('Rho (Weber)')
plt.title(f'JPN {shot} Pol flux over time at different PSI norm values')

timer = 50

# profilo rho a dato tempo: Che sarebbe lo stesso di:w.efit.ftor.plot(t=timer)
rr,ft = w.efit.ftor.get(t=timer)
plt.figure()
plt.plot(rr,ft)
plt.xlabel('PSI norm')
plt.ylabel('Rho (Weber)')
plt.title(f'JPN {shot} Pol flux profile at t= {timer}' )

# Profilo di Rho normalizzato a dato tempo
plt.figure()    
plt.plot(rr,ft/ft[-1])
plt.xlabel('PSI norm')
plt.ylabel('Pol Flux normalized')
plt.title(f'JPN {shot} Normalized poloidal flux profile at t= {timer}' )


# for i,time in enumerate(time_ece):
#     psi = psiKk1[i]
#     time = time_ece[i]
#     fl_int = w.efit.ftor.slice(t=time)
#     rho_ece[i] = rhos[]
###############################################################################
# Prove di calcolo e plots

tTs = w.hrts.te       # Chan Te HRTS  
psiTs = w.hrts.psi    # Channel psi hrts

tEce = w.ecm1.prfl 
zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y
time_ece = tEce.t
time_ts = tTs.t

psiKk1 = np.zeros(tEce.v.shape)
rho_ece = np.zeros(tEce.v.shape)
rho_ts = np.zeros(tTs.v.shape)

for i,time in enumerate(time_ece):
    r, te = tEce.get(t=time)  # Profilo a T=time
    z = np.full_like(r, zKk1) # Retta linea di vista parallela asse r
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
    psiKk1[:,i] = psi
        
    fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i]) # function to intertpolate the flux
    fl_int_ece = fl_int(psi)
    rho_ece[:,i] = np.sqrt(fl_int_ece)

for i,time in enumerate(time_ts):
    psi_th = psiTs.v[:,i]
    fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i]) # function to intertpolate the flux
    fl_int_hrts = fl_int(psi_th)
    rho_ts[:,i] = np.sqrt(fl_int_hrts)

# Profilo temperatura in rho al tempo di indice 350    
plt.figure()
plt.scatter(rho_ece[:,350],tEce.v[:,350])    
    
plt.figure()
plt.scatter(psiTs.v[:,100],tTs.v[:,100])  

###############################################################################
###############################################################################
# Calcolo del profilo a dato tempo 

rEce = r
idx= np.argmin(abs(time_ece - tlim))
psi_ece = psiKk1[:,idx]

idx = np.argmin(abs(time_ts-tlim))
psi_ts = psiTs.v[:,idx]
idd = psi_ts>=0
psi_ts = psi_ts[psi_ts>=0]

tempEce = tEce.slice(t=tlim)
tempTs = tTs.slice(t=tlim)
temp_ts = tempTs.v[idd]
rad_ts = psiTs.r[idd]
# Seleziono gli indicii corrispondenti all'intervallo in psi:psi1-psi2
idxPsiE = ((psi_ece >= psi1) & (psi_ece <= psi2))  # Indici dell'array di PsiEce
idxPsiTs = ((psi_ts >= psi1) & (psi_ts <= psi2))

#############################################

tTs = w.hrts.te       # Chan Te HRTS  
psiTs = w.hrts.psi    # Channel psi hrts
tEce = w.ecm1.prfl 
zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y
time_ece = tEce.t
time_ts = tTs.t
psiKk1 = np.zeros(tEce.v.shape)

idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = tTs.t[idt]
tempTs_ = tTs.v[:,idt]/1000 # in keV
psiTs_ = psiTs.v[:,idt]
errTs_ = w.hrts.dte.v[:,idt]/1000  
# errTs_ = errTs_[errTs_<2]  # prendo solo i pounti in cui errTs è minore di 2 keV 
errTs_[errTs_>2] = 1 # metto a 1 keV l'errore sui punti dove diverge
            
# Ciclo sulle Selezione intervallo psi e media tra psi1 e psi2: 
# Per il HRTS

tempTsM = np.zeros(tempTs_.shape[1])  
psiTsM_ = np.zeros(psiTs_.shape[1])
errTsM = np.zeros(tempTs_.shape[1])
psiTsM = []
xm = []
errXm = []

for i in range(0,(psiTs_.shape[1])):
    mask = (psiTs_[:,i] >= psi1) & (psiTs_[:,i] <= psi2)    # Seleziono gli indici dei tempi di interesse 
    tempTsM[i]  = np.mean(tempTs_[mask,i], axis=0)   # Media aritmetica sui valori di psi selzionati
    psiTsM_[i] = np.mean(psiTs_[mask,i],axis=0)      # media aritmetica della posizione psi
    errTsM[i] = np.mean(errTs_[mask,i], axis=0)
    temp = tempTs_[mask,i]                    # Valroi di temperaatura nelle psi selezionate
    sig = errTs_[mask,i]       # Errore sui singoli valori di temperatura: sigma
    ai = 1/sig**2            # Inverso del quadrato delle sigma
    somma = sum(ai)
    if somma==0:
        somma=1
    xm_ = sum(ai*temp)/somma   # Media dei valori di temperatura pesata sui singoli errori
    xm.append(xm_)               # Valore delle temperature 'attese'/più probabili, nel tempo
    errXm_ = np.sqrt(1/somma) 
    errXm.append(errXm_)
        
### Medie sull'intervallo di psi scelto      
# Calcolo della psi ece per tutti i tempi - psiKk1    
for i,time in enumerate(time_ece):
    r, te = tEce.get(t=time)  # Profilo a T=time
    z = np.full_like(r, zKk1) # Retta linea di vista parallela asse r
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
    psiKk1[:,i] = psi

idt = (time_ece >= tlim1-delta) & (time_ece <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeEce = time_ece[idt]
tempEce = w.ecm1.prfl.v[:,idt]/1000  # in keV
psiEce = psiKk1[:,idt]
errEce = tempEce*0.02

tempEceM = np.zeros(tempEce.shape[1])
psiEceM = np.zeros(psiEce.shape[1])
errEceM = np.zeros(tempEce.shape[1])
# Colcolo delle medie nell'intervallo di psi tra psi1 e psi2 nell'intervallo di tempo scelto
for i in range(0,(psiEce.shape[1])):
    mask1 = (psiEce[:,i] >= psi1) & (psiEce[:,i] <= psi2)    # Seleziono gli indici delle psi di interesse 
    tempEceM[i]  = np.mean(tempEce[mask1,i], axis=0)
    psiEceM[i] = np.mean(psiEce[mask1,i],axis=0)
    errEceM[i] = np.mean(errEce[mask1,i], axis=0)

# Check plot of the PSI-HRTS and PSI-ECE mean values considered for the evaluation of the temperature
plt.figure('Check plot of the averaged PSI values')
# plt.scatter(timeTs,psiTsM)
plt.plot(timeTs,psiTsM_, label='Mean PSI - HRTS')
plt.plot(timeEce,psiEceM, label='Mean PSI - ECE KK1')
plt.xlabel('Time (sec)')
plt.ylabel('PSI')
plt.title('Selected PSI mean values over time')
# plt.ylim(psi1,psi2)    
plt.legend()

####################### Plot in psi
fig03,(ax03) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend')
linew = 0.7
ax03.lw = 0.5
ax03.errorbar(timeTs,tempTsM,  lw = linew, color = 'b', 
             yerr= errTsM, ecolor='g', elinewidth= 0.2, label='HRTS mean')    # arithmetic mean HRTS
ax03.errorbar(timeTs,xm,lw = linew, color='g',
              yerr = errXm, ecolor='g', elinewidth=.3, label='HRTS w-mean') # xm is the weighted mean for HRTS
# ax03.errorbar(timeEce,tempEceM, lw = linew, color='orange', 
             # yerr= errEceM, ecolor='r', elinewidth= 0.2, label='ECE mean')      # arithmetic mean ECE Michelson
ax03.plot(timeEce,tempEceM,lw = linew, color='orange', label='ECE mean' )             
ax03.fill_between(timeEce,tempEceM+errEceM/2,tempEceM-errEceM/2,color = 'darkorange',alpha=0.3)
ax03.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax03.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax03.legend(fontsize=8)
ax03.set_title(f'{shot} - Mean Profile for {psi1}<PSI<{psi2}')
ax03.set_ylabel('Te (keV)') 
ax03.set_xlabel('time (s)')
########################## Plot in Psi seconda versione

####################### Plot in psi
fig05,(ax05) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend')

errXm = np.array(errXm)
linew = 0.7
ax05.lw = 0.5
ax05.plot(timeTs,xm,lw=linew,color='darkolivegreen',label='HRTS w-mean')
ax05.fill_between(timeTs,xm+errXm/2,xm-errXm, color='g',alpha=0.3) # xm is the weighted mean for HRTS
ax05.plot(timeEce,tempEceM,lw = linew, color='orange', label='ECE mean' )             
ax05.fill_between(timeEce,tempEceM+errEceM,tempEceM-errEceM/2,color = 'darkorange',alpha=0.3)
ax05.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax05.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax05.legend(fontsize=8)
ax05.set_title(f'{shot} - Mean Profile for {psi1}<PSI<{psi2}')
ax05.set_ylabel('Te (keV)') 
ax05.set_xlabel('time (s)')

#############################################

dim = timeTs.size
tempEceM_  = signal.resample_poly(tempEceM, up=150, down=402)
errEceM_ = signal.resample_poly(errEceM, up=150, down=402)
# v,t = resample(tempEceM, dim, t=timeEce, axis=0)
# v = signal.resample(tempEceM, dim)
# tempEceM_ = v
# timeEce_ = t

# v,t = signal.resample(errEceM, dim, t=timeEce, axis=0)  #t=timeEce,
# errEce_ = v


def retta(x):
    return x
x = np.linspace(-1,13,dim)   # Range in keV di dove tracciare la retta
y = retta(x)  

# plt.figure(f'{shot}_Tts_vs_Tece: Mean values over {psi1}<PSI<{psi2}', clear=True)
# # plt.scatter(tempEceM_,xm)
# plt.errorbar(tempEceM_,xm, xerr = errEceM_, yerr = errXm, ecolor='g', elinewidth=.3, label='ECE vs TS')
# plt.plot(x,y,'g--', lw=.8)
# plt.xlim(left=4)
# plt.ylim(bottom=4)
# plt.legend()
# plt.title(f'Shot n.{shot} - Te HRTS vs Te Ece-Michelson (Mean values over psi interval)')
# plt.xlabel('Te Ece-Michelson (keV)')
# plt.ylabel('Te HRTS (keV)')


q = xm
p = tempEceM_
#############################################################################
# Plot della Te TS vs Te ECE + retta x=y
plt.figure(f'{shot}_Tts_vs_Tece', clear=True)
# plt.scatter(Tece.v[0,:]/1000,Tts.v[0,:]/1000,s=5,marker='o',label='ECE vs TS')
plt.errorbar(p,q,xerr = errEceM_, yerr = errXm,marker='o', markersize=3,ecolor='g',linestyle='none', 
             elinewidth=.5,label='ECE vs TS')
#plt.scatter(p,q,marker='o', markersize=3,ecolor='g',linestyle='none', label='ECE vs TS')      
plt.plot(x,y,'g--', lw=.8)
plt.xlim(left=4)
plt.ylim(bottom=4)
plt.legend()
plt.title(f'Shot n.{shot} - Te HRTS vs Te Ece-Michelson')
plt.xlabel('Te Ece-Michelson (keV)')
plt.ylabel('Te HRTS (keV)')
# if save == 1: 
#     plt.savefig(f'{shot}_Tts_vs_Tece_Err',dpi=600)










