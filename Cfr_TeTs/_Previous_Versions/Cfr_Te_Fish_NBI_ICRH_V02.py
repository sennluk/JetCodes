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
tlim1 = 47  # limite inferiore selezione tempi - in secondi
tlim2 = 52 
delta = 2   # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3.1     # raggio al quale viene fatta l'analisi (in metri)

w = ppfs(shot)    # Acquisizione del canale relativo allo shot selezionato
# Acquisisco anche i dati sul JPF dei modi MHD
f1 = jetdata(99950, userid='JETJPF',jpf='C1H-G101')


DDAs = ['hrts','lidr','ecm1','NBI','ICRH','C1H-G101','C1H-G102'] # Inutile
#########################
#########################

# Acq dati, selzione punto radiale e intervallo temporale
ts,Ts = w.hrts.te.get(r=rad)    # time  and TS-HRTS electrons temperature
idt = (ts >= tlim1-delta) & (ts <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
TimeTS = ts[idt]
TempTS = Ts[idt]/1000   # in keV

tl,Tl = w.lidr.te.get(r=rad)    # time  and TS-Lidar el. temp.
idt = (tl >= tlim1-delta) & (tl <= tlim2+delta) 
TimeLID = tl[idt]
TempLID = Tl[idt]/1000

te,Te = w.ecm1.prfl.get(r=rad)  # time  and ECE Mich(KK1) el. temp.
idt = (te >= tlim1-delta) & (te <= tlim2+delta)
TimeEM = te[idt]
TempEM = Te[idt]/1000

idt = (f1.t >= tlim1-delta) & (f1.t <= tlim2+delta)
timeF1 = f1.t[idt]
f1 = f1.v[0,idt]

tn,Pn = w.nbi.ptot.get(r=rad)
idt = (tn >= tlim1-delta) & (tn <= tlim2+delta)
TimeNBI = tn[idt]
PNBI = Pn[idt]/10**6  # in MW

ti,Pi = w.icrh.ptot.get(r=rad)
idt = (ti >= tlim1-delta) & (ti <= tlim2+delta)
TimeICRH = ti[idt]
PICRH = Pi[idt]/10**6

# tn1,c1n1 = w.da.c1h-g101



##############
##############
linew = 0.7
# fig,(ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5, sharex=True, num=f'{shot} - Time trends')
fig,(ax0,ax1) = plt.subplots(nrows=2, sharex=True, num=f'{shot} - Time trends')

ax0.lw = 0.5
ax0.plot(TimeTS,TempTS, lw = linew, label='HRTS')
ax0.plot(TimeLID,TempLID, lw = 1, label='Lidar')
ax0.plot(TimeEM,TempEM, lw = 1,label='ECE Michelson')
ax0.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax0.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax0.legend(fontsize=8)
ax0.set_title(f'Shot n.{shot} - Time trends')
ax0.set_ylabel('Te (keV)')  


# ax2.plot(TimeICRH,PICRH, lw = linew, label='P_ICRH')
# ax2.legend(fontsize=8)
# ax2.set_ylabel('P ICRH (MW)')  

ax1.plot(timeF1,f1, lw = linew, label='G1H-G101')
ax1.set_ylabel('C1H-G101')

# ax4.plot(timeF2,f2, lw = linew, label='G1H-G102')
# ax4.set_ylabel('C1H-G102')
# ax4.set_xlabel('time (sec)')
##############
##############
