"""
Created 20/12/2023 - lsenni
Confronto tra i dati ti temperatura elettronica misurati 
con Thomson HRTS (canale hrts.te) e Ece Michelson (canale ecm1.pfrl) 

Figura di confronto riportando ogni diagnostica con il proprio rate di acquisizione
e 
Figura con il rapporto facendo il reshape dell'asse temporale dell'Ece, che e' piu veloce,
su quello dell'hrts che e' piu' lento'


V07 - due plot: 
    1 -Andamenti nel tempo. Ts con errore, Tece no. + Andamento rapporto nel tempo + Andamenti Riscaldamenti e Fishbones nel tempo
    2 - Andamento di Tts(con errore) VS Tece su coordinata radiale 
  --> next step: mettere tutto in coordinate di flusso  
    
OSS: Replicando tutta la prima sezione nella seconda ma levando il delta temporale, si ottiene l'andamento di TTs vs Tece
solamente all'interno dell'intervallo temporale desiderato: tra i due limiti fissati, 
mentre il primo plot riporta gli andamenti in un tempo maggiore, dato da tlim1-delta e tliom2 +delta
"""
import numpy as np
from ppfeg import ppfs 
import matplotlib.pyplot as plt
from scipy.signal import resample

plt.close('all')  
save = 0 # 1--> salva i grafici, 0 non li salva
# Scelgo lo shot, l'intervallo temporale, e il ragio di interesse
# Faccio l'analisi sui tempi scalti, ma visualizzo + e - due secondio
shot = 96994
tlim1 = 47  # limite inferiore selezione tempi - in secondi
tlim2 = 52 
delta = 2   # tempo in più di tmlim1 e in meno di tlim2 su cui viene fatta l'analisi (sec)
rad = 3     # raggio su ciu vien fatta l'analisi (in metri)

w = ppfs(shot)    # Acquisizione del canale relativom allo shot selezionato

Tts = w.hrts.te       # canale temp elettronica HRTS
Errts = w.hrts.dte    # Chann Errors HRTS
Tece = w.ecm1.prfl    # Can El Temp ECE Michelson

# Riporto i dati del Thomson all'interno dell'intervallo temporale
# e sul raggio scelto e gli errori --> commenti nelle versioni precedenti
idt = (Tts.t >= tlim1-delta) & (Tts.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
TimeTS = Tts.t[idt]
TempTS = Tts.v[:,idt]
Errts.t = Errts.t[idt]
Errts.v = Errts.v[:,idt]
idr = np.argmin(abs(Tts.r - rad))    # Seleziono gli indici dei raggi di interesse (il più vicino)
RadTS = Tts.r[idr][np.newaxis]
TempTS = TempTS[idr,:][np.newaxis,:]
Errts.r = Errts.r[idr][np.newaxis]
Errts.v = Errts.v[idr,:][np.newaxis,:]
NBI = w.nbi.ptot  # Total power of the NBI beam
ICRH = w.icrh.ptot # Total power of ICRH


# Riporto i dati del Ece Michelson all'interno dell'intervallo temporale
# e sul raggio scelto
idte = (Tece.t >= tlim1-delta) & (Tece.t <= tlim2+delta)
TimeECE = Tece.t[idte]
TempECE = Tece.v[:,idte]
idre = np.argmin(abs(Tece.r - rad))
RadECE = Tece.r[idre][np.newaxis]
TempECE = TempECE[idre,:][np.newaxis,:]

# Riporto i dati del Di Potenza NBI e ICRH all'interno dell'intervallo temporale (maggiore)
# e sul raggio scelto
idt = (NBI.t >= tlim1-delta) & (NBI.t <= tlim2+delta)
TimeNBI = NBI.t[idt]
Pnbi = NBI.v[:,idt]
# idr = np.argmin(abs(NBI.r - rad))
# RPnbi = NBI.r[idre][np.newaxis]
# Pnbi = NBI.v[idre,:][np.newaxis,:]
  
idt = (ICRH.t >= tlim1-delta) & (ICRH.t <= tlim2+delta)
TimeICRH = ICRH.t[idt]
Picrh = ICRH.v[:,idt]
# idr = np.argmin(abs(ICRH.r - rad))
# RPicrh = ICRH.r[idre][np.newaxis]
# Picrh = ICRH.v[idre,:][np.newaxis,:]

ErrY = Errts.v[0,:]/1000 # Divido per due perchè la barra è doppia, /1000 per keV
ErrY[ErrY>100] = 0       # Controllo che l'errore non sia troppo elevato--> errato

dim = TimeTS.size
v,t = resample(TempECE, dim, t=TimeECE, axis=1)
TempECE = v
TimeECE = t

Ratio = TempTS[0,:]/TempECE[0,:]

# plt.figure(f'{shot} Te vs time') (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
fig,(ax0,ax1,ax2) = plt.subplots(nrows=3, sharex=True, num=f'{shot} - Time trends')
ax0.errorbar(TimeTS,TempTS[0,:]/1000,yerr = ErrY,ecolor='g', elinewidth=.5,label='Tts')
ax0.plot(TimeECE,TempECE[0,:]/1000,label='Tece')
ax0.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax0.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax0.legend()
ax0.set_title(f'Shot n.{shot} - Time trends @ R=3m')
ax0.set_ylabel('Te (keV)')  

ax1.plot(TimeTS,Ratio,label='Ratio')
ax1.axhline(y = 1,ls='--',lw=.8)
ax1.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax1.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax1.legend()
ax1.set_title('Te-HRTS/ Te-ECE')
ax1.set_ylabel('Ratio Te-hrts / Te-ece')

ax2.plot(TimeNBI,Pnbi[0,:]/1e6,lw=1, label='Pnbi')
ax2.axvline(x = tlim1,c='r',ls='--',lw=.5)
ax2.axvline(x = tlim2,c='r',ls='--',lw=.5)
ax2.legend()
ax2.set_title('NBI Power (MW)')
ax2.set_xlabel('time(s)')
ax2.set_ylabel('P NBI (MW)')


if save == 1: 
   plt.savefig(f'{shot}_Te_vs_Time_Err',dpi=600)

def retta(x):
    return x
x = np.linspace(-1,13,dim)   # Range in keV di dove tracciare la retta
y = retta(x)  
#############################################################################
# Ri-seleziono i valori all'interno della zono (più stretta) di interesse
Errots = w.hrts.dte 

idt = (Tts.t >= tlim1) & (Tts.t <= tlim2) # Seleziono gli indici dei tempi di interesse
TimeTS = Tts.t[idt]
TempTS = Tts.v[:,idt]
Errots.t = Errots.t[idt]
Errots.v = Errots.v[:,idt]
idr = np.argmin(abs(Tts.r - rad))    # Seleziono gli indici dei raggi di interesse (il più vicino)
RadTS = Tts.r[idr][np.newaxis]
TempTS = TempTS[idr,:][np.newaxis,:]
Errots.r = Errots.r[idr][np.newaxis]
Errots.v = Errots.v[idr,:][np.newaxis,:]

# Riporto i dati del Ece Michelson all'interno dell'intervallo temporale
# e sul raggio scelto
idte = (Tece.t >= tlim1) & (Tece.t <= tlim2)
TimeECE = Tece.t[idte]
TempECE = Tece.v[:,idte]
idre = np.argmin(abs(Tece.r - rad))
RadECE = Tece.r[idre][np.newaxis]
TempECE = TempECE[idre,:][np.newaxis,:]

ErroY = Errots.v[0,:]/2000 # Divido per due perchè la barra è doppia, /1000 per keV
ErroY[ErroY>100] = 0       # Controllo che l'errore non sia troppo elevato--> errato

dim = TimeTS.size
v,t = resample(TempECE, dim, t=TimeECE, axis=1)
TempECE = v
TimeECE = t

q = TempTS[0,:]/1000
p = TempECE[0,:]/1000
#############################################################################
# Plot della Te TS vs Te ECE + retta x=y
plt.figure(f'{shot}_Tts_vs_Tece', clear=True)
# plt.scatter(Tece.v[0,:]/1000,Tts.v[0,:]/1000,s=5,marker='o',label='ECE vs TS')
plt.errorbar(p,q,yerr = ErroY,marker='o', markersize=3,ecolor='g',linestyle='none', 
             elinewidth=.5,label='ECE vs TS')
#plt.scatter(p,q,marker='o', markersize=3,ecolor='g',linestyle='none', label='ECE vs TS')      
plt.plot(x,y,'g--', lw=.8)
plt.xlim(left=4)
plt.ylim(bottom=4)
plt.legend()
plt.title(f'Shot n.{shot} - Te HRTS vs Te Ece-Michelson')
plt.xlabel('Te Ece-Michelson (keV)')
plt.ylabel('Te HRTS (keV)')
if save == 1: 
    plt.savefig(f'{shot}_Tts_vs_Tece_Err',dpi=600)



