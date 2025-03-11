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
V03: aggiungo le medie pesate con gli errori per i valori in psi dell'HRTS
     per l'ece, assegnado valori tutti con lo stesso errore, mantengo la media aritmetica
03Bis:evito di sovrascrivere le variabile epr riusarle alla fine
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
tlim1 = 47  # limite inferiore selezione tempi - in secondi
tlim2 = 54 
delta = 1   # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.1     # raggio al quale viene fatta l'analisi (in metri)
psi1 = 0.08   # Intervallo in psi su cui mediare
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
####################

sTs = tTs.slice(r=rad)    # time slice for TS-HRTS ò the nearest R to rad
sErrTs = errTs.slice(r=rad)
idxTs = (sTs.t >= tlim1-delta) & (sTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = sTs.t[idxTs]
tempTs = sTs.v[0,idxTs]/1000   # in keV
errTs = sErrTs.v[0,idxTs]
errTs = errTs/1000 # /1000 per keV
errTs[errTs>100] = 0       # Controllo che l'errore non sia troppo elevato--> errato
posTs = sTs.r
print('TS position=',posTs)

sEce = tEce.slice(r=rad)    # time slice for TS-HRTS ò the nearest R to rad
idxEce = (sEce.t >= tlim1-delta) & (sEce.t <= tlim2+delta)
timeEce = sEce.t[idxEce]
tempEce = sEce.v[0,idxEce]/1000
errEce = tempEce*0.02   # Assegno un errore generico del 2% sulla misura
posEce = sEce.r
print('Ece position = ', posEce)



##################### Plot andmento nel tempo ad R fissaato
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

# Ciclo sulle Selezione intervallo psi e media rta psi1 e psi2: 
# Per il HRTS
idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = tTs.t[idt]
tempTs = tTs.v[:,idt]/1000 # in keV
psiTs = psiTs.v[:,idt]
errTs = w.hrts.dte.v[:,idt]/1000 

tempTsM = np.zeros(tempTs.shape[1])  
psiTsM = np.zeros(psiTs.shape[1])
xm = []
errXm = []

for i in range(0,(psiTs.shape[1])):
    mask = (psiTs[:,i] >= psi1) & (psiTs[:,i] <= psi2)    # Seleziono gli indici dei tempi di interesse 
    tempTsM[i]  = np.mean(tempTs[mask,i], axis=0)   # Media aritmetica sui valori di psi selzionati
    psiTsM[i] = np.mean(psiTs[mask,i],axis=0)      # media aritmetica della posizione psi
    temp = tempTs[mask,i]                    # Valroi di temperaatura nelle psi selezionate
    sig = errTs[mask,i]       # Errore sui singoli valori di temperatura: sigma
    sig[sig>100] = 1  # Check valori errori sballati 
    ai = 1/sig**2            # Inverso del quadrato delle sigma
    xm_ = sum(ai*temp)/sum(ai)   # Media dei valori di temperatura pesata sui singoli errori
    xm.append(xm_)               # Valoro delle temperature 'attese'/poiù probabili, nel tempo
    errXm_ = np.sqrt(1/sum(ai)) 
    errXm.append(errXm_)
    
# Check plot of the PSI mean values considered for the evaluation of the temperature
plt.figure()
# plt.scatter(timeTs,psiTsM)
plt.plot(timeTs,psiTsM)
plt.xlabel('Time (sec)')
plt.ylabel('PSI')
plt.title('Selected PSI mean values over time for HRTS')
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
### Medie sull'intervallo di psi scelto
# tempEceM = np.zeros(tempEce.shape[1])
tempEceM = np.zeros_like(tempEce[0:1,:])
psiEceM = np.zeros(psiEce.shape[1])
errEceM = np.zeros(tempEce.shape[1])

for i in range(0,(psiEce.shape[1])):
    mask1 = (psiEce[:,i] >= psi1) & (psiEce[:,i] <= psi2)    # Seleziono gli indici delle psi di interesse 
    tempEceM[0,i]  = np.mean(tempEce[mask1,i], axis=0)
    psiEceM[i] = np.mean(psiEce[mask1,i],axis=0)
    errEceM[i] = np.mean(errEce[mask1,i], axis=0)
    
####################### Plot in psi
fig,(ax0) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend')

ax0.lw = 0.5
ax0.plot(timeTs,tempTsM, lw = linew, label='HRTS mean')
ax0.errorbar(timeTs,xm,lw = linew, yerr = errXm,ecolor='g', elinewidth=.3, label='Tts')
ax0.plot(timeEce,tempEceM, lw = linew, label='ECE mean')
plt.plot(timeTs,xm,lw = linew, label='HRTS w-mean')
ax0.legend(fontsize=8)
ax0.set_title(f'{shot} - Mean Profile for {psi1}<PSI<{psi2}')
ax0.set_ylabel('Te (keV)') 
ax0.set_xlabel('Time (s)')


################################### Resampling and direct comparison - HRTS vs ECE 

# idte = (tEce.t >= tlim1-delta) & (tEce.t <= tlim2+delta)
# TimeECE = tEce.t[idte]
# TempECE = tEce.v[:,idte]
# idre = np.argmin(abs(tEce.r - rad))
# RadECE = TempECE[idre,:][np.newaxis,:]

# # dim = TimeTS.size
# dim = timeTs.size
# v,t = resample(TempECE, dim, t=TimeECE, axis=1)
# TempECE = v
# TimeECE = t



##############

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


