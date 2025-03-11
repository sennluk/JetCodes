"""
Created 11/10/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

Partendo da Pippo06 verso il confronto degli sapri per analisi correlazione fishbones
Provo inserendo le funzioni di egio, anzichè fare i passaggi a mano

OBJ: focalizziamo su DTE3
HRTS e ECM1 con errori (definito sper HRTS, per ECE stimiamo un 2%,
                        massimo errore riscontrato)
Coordinate di equilibrio - Psi
Media su un intervallo di Psi
V01: provo ciclando sulle colonne per il calcolo della media di psi nel valori desiderati
V02: uso 'slice' invece di 'get' rta le funzioni di Egio
V03: aggiungo le medie pesate con gli errori per i valori in psi dell'HRTS
     per l'ece, assegnado valori tutti con lo stesso errore, mantengo la media aritmetica
"""
import numpy as np
import ppfeg
import my_flush
from ppfeg import ppfs, jetdata, V2d
import matplotlib.pyplot as plt
from scipy.signal import resample

#plt.close('all')  
save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
# Faccio l'analisi sui tempi scalti, ma visualizzo + e - due secondio
shot = 96994 # 99950, 99971
tlim1 = 47.5    # limite inferiore selezione tempi - in secondi
tlim2 = 55 
delta = 1     # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.0     # raggio al quale viene fatta l'analisi (in metri)
psi1 = 0.004  # Intervallo in psi su cui mediare
psi2 = 0.084
eP = 0.02  # Relative error assigned to the Ece data 
 
w = ppfs(shot)    # Acquisizione del canale relativo allo shot selezionato
# Acquisisco anche i dati sul JPF dei modi MHD
# f1 = jetdata(shot, userid='JETJPF',jpf='C1H-G101')

tTs = w.hrts.te       # Chan Te HRTS  
errTs = w.hrts.dte    # Chann Errors HRTS  
psiTs = w.hrts.psi    # Channel psi hrts

tEce = w.ecm1.prfl    # Chan El Temp ECE Michelson 
###############################################################################
# Plot time trend at a specific position R --> To have a check

sTs = tTs.slice(r=rad)       # time slice for TS-HRTS --> the nearest R to rad
sPsiTs = tTs.slice(r=rad)
sErrTs = errTs.slice(r=rad)
idt = (sTs.t >= tlim1-delta) & (sTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = sTs.t[idt]
tempTs = sTs.v[0,idt]/1000   # in keV
errTs = sErrTs.v[0,idt]/1000  # in keV
errTs[errTs>100] = 1       # Controllo che l'errore non sia troppo elevato--> errato
posTs = sTs.r
print('TS position=',posTs)

sEce = tEce.slice(r=rad)    # time slice for TS-HRTS ò the nearest R to rad
idt = (sEce.t >= tlim1-delta) & (sEce.t <= tlim2+delta)
timeEce = sEce.t[idt]
tempEce = sEce.v[0,idt]/1000
errEce = tempEce*eP   # eP = relative errorassigned to the Ece Measurement
posEce = sEce.r
print('Ece position = ', posEce)

linew = 0.5
# fig,(ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5, sharex=True, num=f'{shot} - Time trends')
fig,(ax0) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -R- Time trend')

ax0.lw = 0.5
ax0.errorbar(timeTs,tempTs,lw = linew, yerr = errTs,ecolor='g', elinewidth=.2,label='Tts')
ax0.errorbar(timeEce,tempEce,lw = linew, yerr = errEce,ecolor='r', elinewidth=.2,label='Tece')
ax0.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax0.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax0.legend(fontsize=8)
ax0.set_title(f'{shot} - Te Profile at Rts={posTs} m and Rece={posEce}')
ax0.set_ylabel('Te (keV)')  
ax0.set_xlabel('Time (s)')

# Check figure for psi values
# plt.figure()
# plt.plot(psiTs.r,psiTs.v[:,500])

# plt.figure()
# plt.plot(timeTs,sPsiTs.v[0,idt])
plt.figure()
w.hrts.psi.plot(r=rad)
plt.xlabel('time(sec)')
plt.ylabel('Psi')
plt.title(f'Psi HRTS trend over time at R={rad}')
plt.xlim(tlim1,tlim2)
plt.ylim(-0.05,0.05)

tlim = (tlim2+tlim1)/2
g = w.hrts.psi.slice(t=tlim)
plt.figure('PSI HRTS profile')
plt.plot(g.r,g.v, color='g', linewidth=0.5)
idxPsi = ((g.v[:,0] >= psi1) & (g.v[:,0] <= psi2))
plt.plot(g.r[idxPsi],g.v[idxPsi,0],color='r', linewidth=1.5, label='average interval')
plt.xlabel('R(m)')
plt.ylabel('PSI HRTS')
plt.title(f'{shot} - PSI profile over the LoF of HRTS at t=t_mean and the PSI positions of the average')
plt.legend()

###############################################################################
# Plot in PSI
# Acq dati, selezione punto radiale e intervallo temporale

idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = tTs.t[idt]
tempTs_ = tTs.v[:,idt]/1000 # in keV
psiTs_ = psiTs.v[:,idt]
errTs_ = w.hrts.dte.v[:,idt]/1000  
# errTs = errTs[idt]

# Ciclo sulle Selezione intervallo psi e media rta psi1 e psi2: 
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
    
# Check plot of the PSI mean values considered for the evaluation of the temperature
plt.figure('Check plot of the averaged PSI values')
# plt.scatter(timeTs,psiTsM)
plt.plot(timeTs,psiTsM_)
plt.xlabel('Time (sec)')
plt.ylabel('PSI')
plt.title('Selected PSI mean values over time for HRTS')
plt.ylim(psi1,psi2)

############################################################# PSI ECE KK1 calc

time_ax = tEce.t
zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista, 
# la seconda è l'intersezione con l'asse y
psiKk1 = np.zeros(tEce.v.shape)

for i,time in enumerate(time_ax):
    r, te = tEce.get(t=time)  # Profilo a T=time
    z = np.full_like(r, zKk1) # Retta linea di vista parallela asse r
    ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
    # e Prendo il tempo più vicino a quello desiderato
    psi, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
    psiKk1[:,i] = psi

# Plot profilo di psi sulle LoSs del HRTS e dell'ECE
index = np.argmin(abs(time_ax-tlim))   # Indice dell'istante temporale tlim da selezionare
h = psiKk1[:,index]    # PSI profile at the instant tlim
idxPsiE = ((h >= psi1) & (h <= psi2))  # Indici dell'array di PsiEce

plt.figure('PSI Profiles')
plt.plot(g.r,g.v, color='g', linewidth=0.5, label='HRTS LoS')
plt.plot(g.r[idxPsi],g.v[idxPsi,0],color='r', linewidth=1.5, label='HRTS average interval')
plt.plot(r,psiKk1[:,index], linewidth=0.5, label='ECE LoS' )

plt.plot(r[idxPsiE],psiKk1[idxPsiE,index],linewidth=1.5, label='ECE average interval')
plt.xlabel('R(m)')
plt.ylabel('PSI')
plt.title(f'{shot} - PSI profile over the LoS of HRTS and ECE at t=t_mean and the PSI positions of the average')
plt.legend()

######  Average values calculation
idt = (time_ax >= tlim1-delta) & (time_ax <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeEce = time_ax[idt]
tempEce = w.ecm1.prfl.v[:,idt]/1000  # in keV
psiEce = psiKk1[:,idt]
errEce = tempEce*0.02
### Medie sull'intervallo di psi scelto
tempEceM = np.zeros(tempEce.shape[1])
psiEceM = np.zeros(psiEce.shape[1])
errEceM = np.zeros(tempEce.shape[1])

for i in range(0,(psiEce.shape[1])):
    mask1 = (psiEce[:,i] >= psi1) & (psiEce[:,i] <= psi2)    # Seleziono gli indici delle psi di interesse 
    tempEceM[i]  = np.mean(tempEce[mask1,i], axis=0)
    psiEceM[i] = np.mean(psiEce[mask1,i],axis=0)
    errEceM[i] = np.mean(errEce[mask1,i], axis=0)
    
####################### Plot in psi
fig,(ax0) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend')

ax0.lw = 0.5
ax0.errorbar(timeTs,tempTsM,  lw = linew, color = 'b', 
             yerr= errTsM, ecolor='g', elinewidth= 0.2, label='HRTS mean')    # arithmetic mean HRTS
ax0.errorbar(timeTs,xm,lw = linew, color='g',
              yerr = errXm, ecolor='g', elinewidth=.3, label='HRTS w-mean') # xm is the weighted mean for HRTS
ax0.errorbar(timeEce,tempEceM, lw = linew, color='orange', 
             yerr= errEceM, ecolor='r', elinewidth= 0.2, label='ECE mean')      # arithmetic mean ECE Michelson
ax0.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax0.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax0.legend(fontsize=8)
ax0.set_title(f'{shot} - Mean Profile for {psi1}<PSI<{psi2}')
ax0.set_ylabel('Te (keV)') 
ax0.set_xlabel('Time (s)')


################################### Resampling and direct comparison - HRTS vs ECE 


# dim = timeTs.size
# v,t = resample(tempEceM, dim, t=timeEce, axis=1)
# tempEceM = v
# timeEce = t

# v,t = resample(errEce, dim, axis=1)  #t=timeEce,
# errEce = v
# timeEce = t

# def retta(x):
#     return x
# x = np.linspace(-1,13,dim)   # Range in keV di dove tracciare la retta
# y = retta(x)  

# plt.figure(f'{shot}_Tts_vs_Tece: Mean values over {psi1}<PSI<{psi2}', clear=True)
# plt.errorbar(tempEceM,xm,lw = linew, xerr = errEceM, yerr = errXm, ecolor='g', elinewidth=.3, label='ECE vs TS')
# plt.plot(x,y,'g--', lw=.8)
# plt.xlim(left=4)
# plt.ylim(bottom=4)
# plt.legend()
# plt.title(f'Shot n.{shot} - Te HRTS vs Te Ece-Michelson (Mean values over psi interval)')
# plt.xlabel('Te Ece-Michelson (keV)')
# plt.ylabel('Te HRTS (keV)')


