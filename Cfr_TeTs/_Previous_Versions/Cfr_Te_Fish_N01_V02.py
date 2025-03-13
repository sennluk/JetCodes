"""
Created 11/10/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

Partendo da Pippo06 verso il confronto degli sapri per analisi correlazione fishbones
Provo inserendo le funzioni di egio, anzichè fare i passaggi a mano

OBJ: focalizzaimo su DTE3
HRTS e ECM1 con errori (definito sper HRTS, per ECE stimiamo un 2%,
                        massimo errore riscontrato)
Coordinate di equilibrio - Psi
Media su un intervallo di Psi
V01: provo ciclando sulle colonne per il calcolo della media di psi nel valori desiderati
V02: uso 'slice' invece di 'get' rta le funzioni di Egio
"""
import numpy as np
import ppfeg
import my_flush
from ppfeg import ppfs, jetdata, V2d
import matplotlib.pyplot as plt
from scipy.signal import resample

plt.close('all')  
save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
# Faccio l'analisi sui tempi scalti, ma visualizzo + e - due secondio
shot = 104522 #99950, 99971
tlim1 = 48  # limite inferiore selezione tempi - in secondi
tlim2 = 53 
delta = 2   # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.3     # raggio al quale viene fatta l'analisi (in metri)
psi1 = 0.07   # Intervallo in psi su cui mediare
psi2 = 0.10

w = ppfs(shot)    # Acquisizione del canale relativo allo shot selezionato
# Acquisisco anche i dati sul JPF dei modi MHD
# f1 = jetdata(shot, userid='JETJPF',jpf='C1H-G101')

# DDAs = ['hrts','ecm1','C1H-G101','C1H-G102'] # Inutile
#########################
#########################
tTs = w.hrts.te       # Chan Te HRTS
errTs = w.hrts.dte    # Chann Errors HRTS
psiTs = w.hrts.psi    # Channel psi hrts

tEce = w.ecm1.prfl    # Chan El Temp ECE Michelson
###############################################################################
# Plot time trend at a specific position R

sTs = tTs.slice(r=rad)    # time slice for TS-HRTS ò the nearest R to rad
sErrTs = errTs.slice(r=rad)
idt = (sTs.t >= tlim1-delta) & (sTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = sTs.t[idt]
tempTs = sTs.v[0,idt]/1000   # in keV
errTs = sErrTs.v[0,idt]
errTs = errTs/1000 # /1000 per keV
errTs[errTs>100] = 0       # Controllo che l'errore non sia troppo elevato--> errato
posTs = sTs.r
print('TS position=',posTs)

sEce = tEce.slice(r=rad)    # time slice for TS-HRTS ò the nearest R to rad
idt = (sEce.t >= tlim1-delta) & (sEce.t <= tlim2+delta)
timeEce = sEce.t[idt]
tempEce = sEce.v[0,idt]/1000
errEce = tempEce*0.02   # Assegno un errore generico del 2% sulla misura
posEce = sEce.r
print('Ece position = ', posEce)

linew = 0.7
# fig,(ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5, sharex=True, num=f'{shot} - Time trends')
fig,(ax0) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -R- Time trend')

ax0.lw = 0.5
ax0.errorbar(timeTs,tempTs,lw = linew, yerr = errTs,ecolor='g', elinewidth=.3,label='Tts')
ax0.errorbar(timeEce,tempEce,lw = linew, yerr = errEce,ecolor='r', elinewidth=.3,label='Tece')
ax0.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax0.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax0.legend(fontsize=8)
ax0.set_title(f'{shot} - Te Profile at Rts={posTs} m and Rece={posEce}')
ax0.set_ylabel('Te (keV)')  
ax0.set_xlabel('Time (s)')

# ax1.plot(timeF1,f1, lw = linew, label='G1H-G101')
# ax1.set_ylabel('C1H-G101')
###############################################################################
# Plot in PSI
# Acq dati, selezione punto radiale e intervallo temporale

idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = tTs.t[idt]
tempTs = tTs.v[:,idt]/1000 # in keV
psiTs = psiTs.v[:,idt]
errTs = w.hrts.dte.v[:,idt]/1000  
# errTs = errTs[idt]

# Ciclo sulle Selezione intervallo psi e media rta psi1 e psi2:
    
# Per il HRTS

# sig = errTs    # sigma sull'HRTS
# ai = 1/sig**2   # coeff a
# xm = sum(ai*sig)/sum(ai) # valor più probabile 
# errXm = np.sqrt(1/sum(ai))  # Errore da associare al valore più probabile
    
tempTsM = np.zeros(tempTs.shape[1])
psiTsM = np.zeros(psiTs.shape[1])
xm = []
# xm = np.zeros_like(tempTs)
# errXm = np.zeros_like(errTs)

for i in range(0,(psiTs.shape[1])):
    mask = (psiTs[:,i] >= psi1) & (psiTs[:,i] <= psi2)    # Seleziono gli indici dei tempi di interesse 
    tempTsM[i]  = np.mean(tempTs[mask,i], axis=0)
    psiTsM[i] = np.mean(psiTs[mask,i],axis=0)
    temp = tempTs[mask,i]
    # temps.append(temp)
    sig = errTs[mask,i]
    # # sig[i] = errTs[mask][:,i]
    ai = 1/sig**2 
    xm_ = sum(ai*sig)/sum(ai)
    xm.append(xm_)
    # errXm[i] = np.sqrt(1/sum(ai)) 
    
# Check plot of the PSI mean values considered for the evaluation of the temperature
plt.figure()
# plt.scatter(timeTs,psiTsM)
plt.plot(timeTs,psiTsM)
plt.xlabel('Time (sec)')
plt.ylabel('PSI')
plt.title('selected PSI mean values over time')
plt.ylim(psi1,psi2)

######################## PSI KK1 calc

# w = ppfeg.ppfs(shot)       # richiamo i dati dello shot indicato

time_ax = w.ecm1.prfl.t
prfl = w.ecm1.prfl  # Canale del profilo
zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y
psiKk1 = np.zeros(prfl.v.shape)

for i,time in enumerate(time_ax):
    r, te = prfl.get(t=time)  # Profilo a T=time
    z = np.full_like(r, zKk1) # Retta linea di vista parallela asse r
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
    psiKk1[:,i] = psi

idt = (time_ax >= tlim1-delta) & (time_ax <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeEce = time_ax[idt]
tempEce = w.ecm1.prfl.v[:,idt]/1000  # in keV
psiEce = psiKk1[:,idt]
errEce = tempEce*0.02
### Medie sull'intervallo dipsi scelto
tempEceM = np.zeros(tempEce.shape[1])
psiEceM = np.zeros(psiEce.shape[1])

for i in range(0,(psiEce.shape[1])):
    mask1 = (psiEce[:,i] >= psi1) & (psiEce[:,i] <= psi2)    # Seleziono gli indici delle psi di interesse 
    tempEceM[i]  = np.mean(tempEce[mask1,i], axis=0)
    psiEceM[i] = np.mean(psiEce[mask1,i],axis=0)
    
####################### Plot in psi
fig,(ax0) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend')

ax0.lw = 0.5
ax0.plot(timeTs,tempTsM, lw = linew, label='HRTS')
ax0.plot(timeEce,tempEceM, lw = linew, label='ECE')
# ax0.plot(timeEce,tempEce, lw = 1,label='ECE Michelson')
ax0.legend(fontsize=8)
ax0.set_title(f'{shot} - Mean Profile for {psi1}<PSI<{psi2}')
ax0.set_ylabel('Te (keV)') 
ax0.set_xlabel('Time (s)')


###################################
# psiTs_pro = psiTs[mask]
# tempTs = tempTs[mask]

# nrow = mask[:,1].sum()
# ncols = mask.shape[1]
# pippo = psiTs_pro.reshape(nrow,ncols)

# Errts.r = Errts.r[idr][np.newaxis]
# Errts.v = Errts.v[idr,:][np.newaxis,:]



# te,Te = w.ecm1.prfl.get(r=rad)  # time  and ECE Mich(KK1) el. temp.
# idt = (te >= tlim1-delta) & (te <= tlim2+delta)
# TimeEM = te[idt]
# TempEM = Te[idt]/1000

# idt = (f1.t >= tlim1-delta) & (f1.t <= tlim2+delta)
# timeF1 = f1.t[idt]
# f1 = f1.v[0,idt]







##############
##############
# linew = 0.7
# # fig,(ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5, sharex=True, num=f'{shot} - Time trends')
# fig,(ax0,ax1) = plt.subplots(nrows=2, sharex=True, num=f'{shot} - Time trends')

# ax0.lw = 0.5
# ax0.plot(TimeTS,TempTS, lw = linew, label='HRTS')
# ax0.plot(TimeLID,TempLID, lw = 1, label='Lidar')
# ax0.plot(TimeEM,TempEM, lw = 1,label='ECE Michelson')
# ax0.axvline(x = tlim1,c='r',ls='--',lw=.5)
# ax0.axvline(x = tlim2,c='r',ls='--',lw=.5)
# ax0.legend(fontsize=8)
# ax0.set_title(f'Shot n.{shot} - Time trends')
# ax0.set_ylabel('Te (keV)')  


# ax2.plot(TimeICRH,PICRH, lw = linew, label='P_ICRH')
# ax2.legend(fontsize=8)
# ax2.set_ylabel('P ICRH (MW)')  

# ax1.plot(timeF1,f1, lw = linew, label='G1H-G101')
# ax1.set_ylabel('C1H-G101')

# ax4.plot(timeF2,f2, lw = linew, label='G1H-G102')
# ax4.set_ylabel('C1H-G102')
# ax4.set_xlabel('time (sec)')
##############
##############
