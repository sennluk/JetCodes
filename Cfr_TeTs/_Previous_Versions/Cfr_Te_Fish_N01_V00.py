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
Media su un intervallo
NOTA: prendo due intervalli temporali diversi per i profili in R e in PSI

"""
import numpy as np
from ppfeg import ppfs, jetdata, V2d
import matplotlib.pyplot as plt
from scipy.signal import resample

plt.close('all')  
save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
# Faccio l'analisi sui tempi scalti, ma visualizzo + e - due secondio
shot = 99971 #99950
tlim1 = 48  # limite inferiore selezione tempi - in secondi
tlim2 = 53 
delta = 2   # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.3     # raggio al quale viene fatta l'analisi (in metri)
psi1 = 0.07   # Intervallo in psi su cui mediare
psi2 = 0.12

w = ppfs(shot)    # Acquisizione del canale relativo allo shot selezionato
# Acquisisco anche i dati sul JPF dei modi MHD
f1 = jetdata(104520, userid='JETJPF',jpf='C1H-G101')

DDAs = ['hrts','ecm1','C1H-G101','C1H-G102'] # Inutile
#########################
#########################
tTs = w.hrts.te       # Chan Te HRTS
errTs = w.hrts.dte    # Chann Errors HRTS
psiTs = w.hrts.psi    # Channel psi hrts

tEce = w.ecm1.prfl    # Chan El Temp ECE Michelson

#######################

# plt.figure(f'{shot} Te vs time') (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
ts,Ts = w.hrts.te.get(r=rad)    # time  and TS-HRTS electrons temperature
idt = (ts >= tlim1-delta) & (ts <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = ts[idt]
tempTs = Ts[idt]/1000   # in keV

te,Te = w.ecm1.prfl.get(r=rad)  # time  and ECE Mich(KK1) el. temp.
idt = (te >= tlim1-delta) & (te <= tlim2+delta)
timeEce = te[idt]
tempEce = Te[idt]/1000

linew = 0.7
# fig,(ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5, sharex=True, num=f'{shot} - Time trends')
fig,(ax0) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -R- Time trend')

ax0.lw = 0.5
ax0.plot(timeTs,tempTs, lw = linew, label='HRTS')
# ax0.plot(timeEce,tempEce, lw = 1,label='ECE Michelson')
ax0.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax0.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax0.legend(fontsize=8)
ax0.set_title(f'{shot} - Profile at R={rad} m')
ax0.set_ylabel('Te (keV)')  
ax0.set_xlabel('Time (s)')

# ax1.plot(timeF1,f1, lw = linew, label='G1H-G101')
# ax1.set_ylabel('C1H-G101')
#######################

# Acq dati, selzione punto radiale e intervallo temporale
idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
timeTs = tTs.t[idt]
tempTs = tTs.v[:,idt]
psiTs = psiTs.v[:,idt]
errTs = errTs.t[idt]

# Ciclo sulle Selezione intervallo psi e media rta psi1 e psi2:
mask = (psiTs[:,100] >= psi1) & (psiTs[:,100] <= psi2)    # Seleziono gli indici dei tempi di interesse 
pippo = tempTs[mask]
pluto = psiTs[mask]

tempTsM = np.mean(tempTs[mask], axis=0)
psiTsM = np.mean(psiTs[mask],axis=0)

plt.figure()
plt.scatter(timeTs,psiTsM)

####################### Plot in psi
fig,(ax0) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend')

ax0.lw = 0.5
ax0.plot(timeTs,tempTsM, lw = linew, label='HRTS')
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
