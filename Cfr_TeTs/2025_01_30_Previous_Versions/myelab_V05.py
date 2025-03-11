"""
Created on Thu Mar 14 10:29:18 2024
@author: lsenni. Documentazione su file word: Note_Codice_Cfr_Te_ts

Introduco il calcolo di rho solamente nell'intervallo t1-Delta - t2 + Delta         
"""
import numpy as np
import my_flush
# from ppfeg import ppfs
import matplotlib.pyplot as plt
from scipy import signal,interpolate

# DTE3 shot list: 104990,991,994,995,999   
# 104520,521,522,523,524,526
  
##################################################
# Ricostriuscio la posizione dell'asse magnetico 
def magax(d,vars):
    shot = d['shot']
    tlim1 = d['tlim1'] 
    tlim2 = d['tlim2']
    w = vars['w']
    ####################
    plt.figure('Magnetic axis position over time')
    w.efit.rmag.plot()
    plt.axvline(x = tlim1,c='r',ls='--',lw=.5, label='t-lim1')
    plt.axvline(x = tlim2,c='r',ls='--',lw=.5, label = 't-lim2')
    plt.xlabel('time (sec)')
    plt.ylabel('R(m)')
    plt.title(f'JPN {shot} - Radial position of the magnetic axis from EFIT')
    plt.tight_layout()
    
    plt.figure()
    w.efit.zmag.plot()
    plt.axvline(x = tlim1,c='r',ls='--',lw=.5, label='t-lim1')
    plt.axvline(x = tlim2,c='r',ls='--',lw=.5, label = 't-lim2')
    plt.xlabel('time (sec)')
    plt.ylabel('z(m)')
    plt.title(f'JPN {shot} - Poloidal position of the magnetic axis from EFIT')
    plt.tight_layout()
    
    pippo=w.efit.rmag
    pluto = w.efit.zmag
    plt.figure()
    plt.plot(pippo.v[0,:],pluto.v[0,:],lw=0.7)
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Mag_Ax_Pos.pdf',dpi=300)
              
    return 1
##################################################
# Time trend at R=rad[']
def tprof(d,vars):
    shot = d['shot']
    tlim1 = d['tlim1'] 
    tlim2 = d['tlim2']
    rad = d['rad']
    delta = d['delta']
    eP = d['eP']
    tTs = vars['tTs']      # Chan Te HRTS  
    errTs = vars['errTs']    # Chann Errors HRTS
    tEce = vars['tEce']    # Chan El Temp ECE Michelson 
    ####################
    # Plot time trend at a specific position R --> To have a check
    
    sTs = tTs.slice(r=rad)       # time slice for TS-HRTS --> the nearest R to rad
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
    fig00,(ax00) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Time trend at R')
    ax00.lw = 0.5
    ax00.errorbar(timeTs,tempTs,lw = linew, yerr = errTs,ecolor='g', elinewidth=.2,label='Tts')
    ax00.errorbar(timeEce,tempEce,lw = linew, yerr = errEce,ecolor='r', elinewidth=.2,label='Tece')
    ax00.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax00.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax00.legend(fontsize=8)
    ax00.set_title(f'{shot} - Te time trend at Rts={posTs} m and Rece={posEce}')
    ax00.set_ylabel('Te (keV)')  
    ax00.set_xlabel('Time (s)')
    fig00.tight_layout()
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(t)_at_Rad.pdf',dpi=300)
            
    return 1
##################################################

def psicalc(d, vars):
    shot = d['shot']
    w = vars['w']
    ####################
    # Prove di calcolo e plots
    tEce = w.ecm1.prfl 
    zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista (in metri), 
    # la seconda è l'intersezione con l'asse y. 'Position of antenna aperture centre (z,majR) '
    rEce = tEce.r
    z = np.full_like(rEce, zKk1) # Retta linea di vista parallela asse r
    
    timeEce = tEce.t

    psiKk1 = np.zeros(tEce.v.shape)
    
    for i,time in enumerate(timeEce):
        ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
        # e Prendo il tempo più vicino a quello desiderato
        psi, _ = my_flush.Flush_getFlux(rEce*100, z*100) # Le coordinate vanne messe in cm!
        psiKk1[:,i] = psi
   
    vars['psiKk1'] = psiKk1      
    
    return 1
##################################################

def psifig(d,vars):
    shot = d['shot']
    psi1 = d['psi1']
    psi2 = d['psi2']
    tlim = vars['tlim']
    tTs = vars['tTs']       # Chan. Te HRTS  
    tEce = vars['tEce']
    psiTs = vars['psiTs']    # Channel psi hrts
    psiKk1 = vars['psiKk1']
   
    ####################
    timeEce = tEce.t
    timeTs = tTs.t
    # Calcolo del profilo a dato tempo 

    idts = np.argmin(abs(timeTs - tlim))
    psiTs_s = np.absolute(psiTs.v[:,idts]) # Abs del profilo psi hrts al tempo t=tlim   
    
    ide= np.argmin(abs(timeEce - tlim))
    psiEce_s = psiKk1[:,ide]     # profilo psi ece al tempo t=tlim          

    tempTs_s = tTs.slice(t=tlim)     # _s: slice at time..
    tempEce_s = tEce.slice(t=tlim)
    rTs = psiTs.r
    rEce = tEce.r
    # Seleziono gli indici corrispondenti all'intervallo in psi:psi1-psi2
    idxPsiTs = ((psiTs_s >= psi1) & (psiTs_s <= psi2))
    idxPsiE = ((psiEce_s >= psi1) & (psiEce_s <= psi2))  # Indici dell'array di PsiEce
    leftlim = min(min(rTs[idxPsiTs]), min(rEce[idxPsiE]) )
    rightlim = max(max(rTs[idxPsiTs]), max(rEce[idxPsiE]) )
    
    fig01,(ax001, ax01) = plt.subplots(nrows=2, sharex = False, num = f'JPN {shot} PSI profiles')
    ax001.plot(rTs, psiTs_s, linewidth=0.5, color='green', label= 'psi hrts')
    ax001.plot(rEce, psiEce_s, linewidth=0.5, color='blue', label='psi ece')
    ax001.set_ylabel('PSI')
    ax001.legend()
    ax001.axhline(y = psi1, c='r',ls='--',lw=.4)
    ax001.axhline(y = psi2,c='r',ls='--',lw=.4)
    ax01.plot(rTs, psiTs_s, linewidth=0.5, color='green', marker='o', ms=0.8, label= 'psi hrts')
    ax01.plot(rTs[idxPsiTs], psiTs_s[idxPsiTs], linewidth=1.5, color='red', marker='o', ms=3, label=f'{psi1}<psi hrts<{psi2}')
    ax01.plot(rEce, psiEce_s, linewidth=0.5, color='blue', marker='*', ms=0.8, label='psi ece')
    ax01.plot(rEce[idxPsiE], psiEce_s[idxPsiE], linewidth=1.5, color='orange', marker='*', ms=3,  label=f'{psi1}<psi ece<{psi2}')
    ax01.axhline(y = psi1, c='r',ls='--',lw=.4)
    ax01.axhline(y = psi2,c='r',ls='--',lw=.4)
    ax01.set_xlim(left = leftlim-0.1, right = rightlim+0.1)
    ax01.set_ylim(bottom = psi1-0.01, top = psi2 + 0.01)
    ax01.set_xlabel('R (m)')
    ax01.set_ylabel('PSI')
    ax01.legend()
    ax001.set_title(f'JPN {shot} PSI(R) for HRTS and ECE-KK1 at t={tlim} sec')
    fig01.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_PSI(R)_at_t.pdf',dpi=300)
            
        
    fig02, ax02 = plt.subplots(nrows=1, sharex=True, num = f'JPN {shot} Te vs PSI')
    ax02.scatter(psiTs_s, tempTs_s.v/1000, marker='o', lw = 0.7, facecolors='none', edgecolors='darkorange', label='Te hrts') #s= 5, 
    ax02.scatter(psiEce_s, tempEce_s.v/1000, marker='1', lw = 0.7, color = 'olive' ,label='Te ece - kk1') # s= 10,
    ax02.set_xlabel('PSI')
    ax02.set_ylabel('Electron Temperature (keV)')
    ax02.legend()
    ax02.set_title(f'JPN {shot} Te(PSI) profile at t={tlim} sec')
    fig02.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(psi)_at_t.pdf',dpi=300)
            
    print('PSI1 HRTS = ', round(psiTs_s[idxPsiTs][0],3), ' - R1 HRTS = ', round(rTs[idxPsiTs][0],3))
    print('PSI2 HRTS = ', round(psiTs_s[idxPsiTs][-1],3), ' - R2 HRTS = ', round(rTs[idxPsiTs][-1],3))
    print('Number of PSI-HRTS Values = ', psiTs_s[idxPsiTs].size)
    
    print('PSI1 ECE = ', round(psiEce_s[idxPsiE][0],3), ' - R1 ECE = ', round(rEce[idxPsiE][0],3))
    print('PSI2 ECE = ', round(psiEce_s[idxPsiE][-1],3), ' - R2 ECE = ', round(rEce[idxPsiE][-1],3))
    print('Number of PSI-ECE Values = ', psiEce_s[idxPsiE].size)

    return 1 

##################################################

def rhocalc(d, vars): 
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    delta = d['delta']
    w= vars['w']
    psiKk1 = vars['psiKk1']
    tTs = vars['tTs']       # Chan. Te HRTS  
    psiTs = vars['psiTs']    # Channel psi hrts
    tEce = vars['tEce']
    ####################
 
    idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ts = tTs.t[idt]
    temp_ts = tTs.v[:,idt]
    psi_ts = psiTs.v[:,idt]
    
    idt = (tEce.t >= tlim1-delta) & (tEce.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ece = tEce.t[idt] 
    temp_ece = tEce.v[:,idt]
    psi_ece = psiKk1[:,idt]
    
    rFlux = w.efit.ftor.r  # Posizioni dei valori di flusso
    vFlux = w.efit.ftor.v  # valori di flusso
    
    rho_ts = np.zeros(temp_ts.shape)
    rho_ece = np.zeros(temp_ece.shape)
     
    for i,time in enumerate(time_ts):
        psi_th = psi_ts[:,i]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux[:,i]/vFlux[-1,i]) # function to intertp the flux
        fl_int_hrts = fl_int(psi_th)
        rho_ts[:,i] = np.sqrt(fl_int_hrts)
        
    for i,time in enumerate(time_ece):
        psi = psi_ece[:,i]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux[:,i]/vFlux[-1,i])  # function to intertp the flux
        fl_int_ece = fl_int(psi)
        rho_ece[:,i] = np.sqrt(fl_int_ece)
    
    vars['time_ts'] = time_ts   # time interval t1-delta<t<t2+delta   
    vars['time_ece'] = time_ece 
    vars['rho_ts'] = rho_ts
    vars['rho_ece'] = rho_ece
        
    return  1
##################################################        

def rhofig(d, vars): 
    shot = d['shot']
    rho1 = d['rho1']
    rho2 = d['rho2']
    tlim = vars['tlim']
    tTs = vars['tTs']  
    tEce = vars['tEce'] 
    time_ts = vars['time_ts']
    time_ece = vars['time_ece']
    rho_ts = vars['rho_ts']
    rho_ece = vars['rho_ece']
    ####################
        
    # Calcolo del profilo a dato tempo - slices
    tempTs_s = tTs.slice(t=tlim)
    tempEce_s = tEce.slice(t=tlim)
        
    idts = np.argmin(abs(time_ts - tlim))
    ide = np.argmin(abs(time_ece - tlim))       
    
    fig08, ax08 = plt.subplots(nrows=1, sharex=True, num = f'JPN {shot} Te vs RHO')
    ax08.scatter(rho_ts[:,idts], tempTs_s.v/1000, marker='o', lw = 0.7, facecolors='none', edgecolors='darkorange', label='HRTS')
    ax08.scatter(rho_ece[:,ide],tempEce_s.v/1000, marker='1', lw = 0.7, color = 'olive', label='ECE-KK1')
    ax08.set_title(f'JPN. {shot} - HRTS and Ece Te vs rho - Profile at t={tlim}')
    ax08.set_xlabel('rho')    
    ax08.set_ylabel('Te (keV)')
    plt.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(rho).pdf',dpi=300)
            
    
    # Seleziono gli indici corrispondenti all'intervallo in rho: rho1-rho2
    idts = np.argmin(abs(time_ts - tlim))
    rho_ts_s = rho_ts[:,idts] # Profilo rho hrts al tempo t=tlim   
    
    ide= np.argmin(abs(time_ece - tlim))
    rho_ece_s = rho_ece[:,ide]     # profilo rho ece al tempo t=tlim          

    rTs = tTs.r
    rEce = tEce.r
    
    idxRhoTs = ((rho_ts_s >= rho1) & (rho_ts_s <= rho2))
    idxRhoE = ((rho_ece_s >= rho1) & (rho_ece_s <= rho2))  # Indici dell'array di PsiEce
    leftlim = min(min(rTs[idxRhoTs]), min(rEce[idxRhoE]))
    rightlim = max(max(rTs[idxRhoTs]), max(rEce[idxRhoE]))
    
    fig03, (ax003, ax03) = plt.subplots(nrows=2, sharex=False, num = f'JPN {shot} Normalized RHO profiles')
    ax003.plot(rTs, rho_ts_s, linewidth=0.5, color='green', label='Rho Hrts LoS')
    ax003.plot(rEce, rho_ece_s, linewidth=0.5, color='blue', label='Rho Ece LoS')
    ax003.axhline(y = rho1, c='r',ls='--',lw=.4)
    ax003.axhline(y = rho2,c='r',ls='--',lw=.4)
    ax003.set_ylabel('Normalized RHO')
    ax003.legend()
    ax03.plot(rTs, rho_ts_s, '-o', ms=1, linewidth=0.5, color='green', label='Rho Hrts LoS')
    ax03.plot(rTs[idxRhoTs], rho_ts_s[idxRhoTs], '-o', ms=3,  linewidth=1.5, color='red', label=f'{rho1}<rho hrts<{rho2}')
    ax03.plot(rEce, rho_ece_s, '-*', ms=1, linewidth=0.5, color='blue', label='Rho Ece LoS')
    ax03.plot(rEce[idxRhoE], rho_ece_s[idxRhoE], '-*', ms=3, linewidth=1.5, color='orange', label='{rho1}<Rho Ece<{rho2}')
    ax03.axhline(y = rho1, c='r',ls='--',lw=.4)
    ax03.axhline(y = rho2,c='r',ls='--',lw=.4)
    ax03.set_xlim(left = leftlim-0.1, right = rightlim+0.1)
    ax03.set_ylim(bottom = rho1-0.05, top = rho2 + 0.01)
    ax03.set_xlabel('R (m)')
    ax03.set_ylabel('Normalized RHO')
    ax03.legend()
    ax003.set_title(f'JPN {shot} RHO(R) for HRTS and ECE-KK1 at t={tlim} sec')
    fig03.tight_layout()
    
    if d['savefigs'] == 1: 
       plt.savefig(d['mypath']+f'{shot}_RHO(R)_at_t.pdf',dpi=300)
              
    print('RHO1 HRTS = ', round(rho_ts_s[idxRhoTs][0],3), ' - R1 HRTS = ', round(rTs[idxRhoTs][0],3))
    print('RHO2 HRTS = ', round(rho_ts_s[idxRhoTs][-1],3), ' - R2 HRTS = ', round(rTs[idxRhoTs][-1],3))
    print('Number of RHO-HRTS Values = ', rho_ts_s[idxRhoTs].size)
    
    print('RHO1 ECE = ', round(rho_ece_s[idxRhoE][0],3), ' - R1 ECE = ', round(rEce[idxRhoE][0],3))
    print('RHO2 ECE = ', round(rho_ece_s[idxRhoE][-1],3), ' - R2 ECE = ', round(rEce[idxRhoE][-1],3))
    print('Number of RHO-ECE Values = ', rho_ece_s[idxRhoE].size)
    
    rhoPPF = vars['w'].hrts.rho.slice(t=tlim)
    rhoPPF = np.abs(rhoPPF.v)
    plt.figure()
    plt.plot(rTs, rho_ts_s, label='RHO')
    plt.plot(rTs, rhoPPF, label='RHO FLUSH')
    plt.xlabel('R(m)')
    plt.ylabel('RHO')
    plt.legend()
    plt.title(f'JPN {shot} - Check plot: RHO vs RHO FLUSH')
  
    
    return 1
    
##################################################
# Calcolo delle medie nell'intervallo di psi scelto tra gli estremi psi1 e psi2

def meancalc(d, vars): 
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    delta = d['delta']
    psi1 = d['psi1']
    psi2 = d['psi2']
    eP = d['eP']

    tTs = vars['tTs']  
    errTs = vars['errTs']
    tEce = vars['tEce'] 
    psiTs = vars['psiTs']
    psiKk1 = vars['psiKk1']

    ####################
      
    timeEce = tEce.t
    
    idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ts = tTs.t[idt]
    temp_ts = tTs.v[:,idt]/1000 # in keV
    psi_ts = psiTs.v[:,idt]
    err_ts = errTs.v[:,idt]/1000  
    err_ts[err_ts>2] = 1 # metto a 1 keV l'errore sui punti dove diverge
                
    # Ciclo sulle Selezione intervallo psi e media tra psi1 e psi2: 
    # Per il HRTS
    
    temp_tsM = np.zeros(temp_ts.shape[1])  
    psi_tsM_ = np.zeros(psi_ts.shape[1])
    err_tsM = np.zeros(temp_ts.shape[1])
    xm = []
    err_xm = []
    
    for i in range(0,(psi_ts.shape[1])):
        mask = (psi_ts[:,i] >= psi1) & (psi_ts[:,i] <= psi2)    # Seleziono gli indici dei tempi di interesse 
        temp_tsM[i]  = np.mean(temp_ts[mask,i], axis=0)   # Media aritmetica sui valori di psi selzionati
        psi_tsM_[i] = np.mean(psi_ts[mask,i],axis=0)      # media aritmetica della posizione psi
        err_tsM[i] = np.mean(err_ts[mask,i], axis=0)
        temp = temp_ts[mask,i]                    # Valroi di temperaatura nelle psi selezionate
        sig = err_ts[mask,i]       # Errore sui singoli valori di temperatura: sigma
        ai = 1/sig**2            # Inverso del quadrato delle sigma
        somma = sum(ai)
        if somma==0:
            somma=1
        xm_ = sum(ai*temp)/somma   # Media dei valori di temperatura pesata sui singoli errori
        xm.append(xm_)               # Valore delle temperature 'attese'/più probabili, nel tempo
        err_xm_ = np.sqrt(1/somma) 
        err_xm.append(err_xm_)
            
    ### Medie sull'intervallo di psi scelto      
    
    idt = (timeEce >= tlim1-delta) & (timeEce <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ece = timeEce[idt]
    temp_ece = tEce.v[:,idt]/1000  # in keV
    psi_ece = psiKk1[:,idt]
    err_ece = temp_ece*eP
    
    temp_eceM = np.zeros(temp_ece.shape[1])
    psi_eceM = np.zeros(psi_ece.shape[1])
    err_eceM = np.zeros(temp_ece.shape[1])
    # Colcolo delle medie nell'intervallo di psi tra psi1 e psi2 nell'intervallo di tempo scelto
    for i in range(0,(psi_ece.shape[1])):
        mask1 = (psi_ece[:,i] >= psi1) & (psi_ece[:,i] <= psi2)    # Seleziono gli indici delle psi di interesse 
        temp_eceM[i]  = np.mean(temp_ece[mask1,i], axis=0)
        psi_eceM[i] = np.mean(psi_ece[mask1,i],axis=0)
        err_eceM[i] = np.mean(err_ece[mask1,i], axis=0)
    
    # Check plot of the PSI-HRTS and PSI-ECE mean values considered for the evaluation of the temperature
    plt.figure('Check plot of the averaged PSI values')
    # plt.scatter(timeTs,psiTsM)
    plt.plot(time_ts,psi_tsM_, label='Mean PSI - HRTS')
    plt.plot(time_ece,psi_eceM, label='Mean PSI - ECE KK1')
    plt.xlabel('Time (sec)')
    plt.ylabel('PSI')
    plt.title('Selected PSI mean values over time - {psi1}<PSI<{psi2}')
    # plt.ylim(psi1,psi2)    
    plt.legend()
    
    vars['temp_tsM'] = temp_tsM
    vars['err_tsM'] = err_tsM
    vars['xm'] = xm
    vars['err_xm'] = err_xm
    vars['temp_eceM'] = temp_eceM
    vars['err_eceM'] = err_eceM
    
    return 1 
 
###############################################################    
def fig_psi_mean(d, vars): 
    shot = d['shot']
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    delta = d['delta']
    psi1 = d['psi1']
    psi2 = d['psi2']
    tTs = vars['tTs']  
    tEce = vars['tEce'] 
    temp_tsM = vars['temp_tsM']
    err_tsM = vars['err_tsM']
    xm = vars['xm']
    err_xm = vars['err_xm']
    temp_eceM = vars['temp_eceM']
    err_eceM = vars['err_eceM']
    ####################
    
    idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ts = tTs.t[idt]
    
    idt = (tEce.t >= tlim1-delta) & (tEce.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ece = tEce.t[idt]
    # Plot in psi e rho
    
    # PLOT CON ERRORBARS e dato con media aritmetica e pesata
    fig04, ax04 = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend - Errorbars')
    linew = 0.7
    ax04.lw = 0.5
    ax04.errorbar(time_ts,temp_tsM,  lw = linew, color = 'b', 
                 yerr= err_tsM, ecolor='g', elinewidth= 0.2, label='HRTS mean')    # arithmetic mean HRTS
    ax04.errorbar(time_ts,xm,lw = linew, color='g',
                  yerr = err_xm, ecolor='g', elinewidth=.3, label='HRTS w-mean') # xm is the weighted mean for HRTS
    ax04.errorbar(time_ece,temp_eceM, lw = linew, color='orange', 
                  yerr= err_eceM, ecolor='r', elinewidth= 0.2, label='ECE mean')      # arithmetic mean ECE Michelson
    #ax03.plot(time_ece,temp_eceM,lw = linew, color='orange', label='ECE mean' )             
    ax04.fill_between(time_ece,temp_eceM+err_eceM/2,temp_eceM-err_eceM/2,color = 'darkorange',alpha=0.3)
    ax04.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax04.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax04.legend(fontsize=8)
    ax04.set_title(f'JPN {shot} - Mean (Ece and Hrts) and w-Mean(Hrts) for {psi1}<PSI<{psi2}')
    ax04.set_ylabel('Te (keV)') 
    ax04.set_xlabel('time (s)')
    fig04.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te_MEANS(psi).pdf',dpi=300)
               
    ########################## Plot in Psi seconda versione
    
    fig05, ax05 = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend - Bands')
    # err_xm va diviso per due o no?, cosa FA L'ERROR BAR DEL PLOT PRECEDENTE?
    err_xm = np.array(err_xm)
    linew = 0.7
    ax05.lw = 0.5
    ax05.plot(time_ts,xm,lw=linew,color='darkolivegreen',label='HRTS w-mean')
    ax05.fill_between(time_ts,xm+err_xm,xm-err_xm, color='g',alpha=0.3) # xm is the weighted mean for HRTS
    ax05.plot(time_ece,temp_eceM,lw = linew, color='orange', label='ECE mean' )             
    ax05.fill_between(time_ece,temp_eceM+err_eceM,temp_eceM-err_eceM/2,color = 'darkorange',alpha=0.3)
    ax05.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax05.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax05.legend(fontsize=8)
    ax05.set_title(f'JPN {shot} - Ece-Mean and HRTS w-Mean for {psi1}<PSI<{psi2} trends vs time ')
    ax05.set_ylabel('Te (keV)') 
    ax05.set_xlabel('time (s)')
    fig05.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te_MEANS(psi)_BANDS.pdf',dpi=300)
            
    return 1

##################################################
# Calcolo delle medie nell'intervallo di psi sceltoi tra gli estremi psi1 e psi2

def meancalc_12(d, vars): 
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    psi1 = d['psi1']
    psi2 = d['psi2']
    eP = d['eP']
    tTs = vars['tTs']  
    errTs = vars['errTs']
    tEce = vars['tEce'] 
    psiTs = vars['psiTs']
    psiKk1 = vars['psiKk1']
    ####################
      
    timeEce = tEce.t
    
    idt = (tTs.t >= tlim1) & (tTs.t <= tlim2) # Seleziono gli indici dei tempi di interesse
    temp_ts = tTs.v[:,idt]/1000 # in keV
    psi_ts = psiTs.v[:,idt]
    err_ts = errTs.v[:,idt]/1000  
    err_ts[err_ts>2] = 1 # metto a 1 keV l'errore sui punti dove diverge
                
    # Ciclo sulle Selezione intervallo psi e media tra psi1 e psi2: 
    # Per il HRTS
    
    xm12 = []
    err_xm12 = []
    
    for i in range(0,(psi_ts.shape[1])):
        mask = (psi_ts[:,i] >= psi1) & (psi_ts[:,i] <= psi2)    # Seleziono gli indici dei tempi di interesse 
        temp = temp_ts[mask,i]                    # Valori di temperaatura nelle psi selezionate
        sig = err_ts[mask,i]       # Errore sui singoli valori di temperatura: sigma
        ai = 1/sig**2            # Inverso del quadrato delle sigma
        somma = sum(ai)
        if somma==0:
            somma=1
        xm_ = sum(ai*temp)/somma   # Media dei valori di temperatura pesata sui singoli errori
        xm12.append(xm_)               # Valore delle temperature 'attese'/più probabili, nel tempo
        err_xm_ = np.sqrt(1/somma) 
        err_xm12.append(err_xm_)
            
    ### Medie sull'intervallo di psi scelto      
    
    idt = (timeEce >= tlim1) & (timeEce <= tlim2) # Seleziono gli indici dei tempi di interesse
    temp_ece = tEce.v[:,idt]/1000  # in keV
    psi_ece = psiKk1[:,idt]
    err_ece = temp_ece*eP
    
    temp_eceM12 = np.zeros(temp_ece.shape[1])
    psi_eceM12 = np.zeros(psi_ece.shape[1])
    err_eceM12 = np.zeros(temp_ece.shape[1])
    # Colcolo delle medie nell'intervallo di psi tra psi1 e psi2 nell'intervallo di tempo scelto
    for i in range(0,(psi_ece.shape[1])):
        mask1 = (psi_ece[:,i] >= psi1) & (psi_ece[:,i] <= psi2)    # Seleziono gli indici delle psi di interesse 
        temp_eceM12[i]  = np.mean(temp_ece[mask1,i], axis=0)
        psi_eceM12[i] = np.mean(psi_ece[mask1,i],axis=0)
        err_eceM12[i] = np.mean(err_ece[mask1,i], axis=0)

    vars['xm12'] = xm12
    vars['err_xm12'] = err_xm12
    vars['temp_eceM12'] = temp_eceM12
    vars['err_eceM12'] = err_eceM12
    
    return 1 
 
###############################################################   

def fig_cfr_psi(d, vars): #shot, w, tlim1, tlim2, delta, psi1, psi2, xm, err_xm, temp_eceM, err_eceM):  
    shot = d['shot']
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    psi1 = d['psi1']
    psi2 = d['psi2']
    xm12 = vars['xm12']
    err_xm12 = vars['err_xm12']
    temp_eceM12 = vars['temp_eceM12']
    err_eceM12 = vars['err_eceM12']
    ####################
    
    dimTs = len(xm12)
    dimEce = len(temp_eceM12)
    
    temp_eceM_R  = signal.resample_poly(temp_eceM12, up=dimTs, down=dimEce) # posizioni ricampionate dei valori TeEce
    err_eceM_R = signal.resample_poly(err_eceM12, up=dimTs, down=dimEce)    # posizioni ricampionate dei valori Err TeEce
 
    def retta(x):
        return x
    x = np.linspace(-1,13,dimTs)   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    q = xm12            # Valori TeHRTS media-pesata calcolata sull'intervallo di psi specificato
    p = temp_eceM_R   # Valori TeECE media calcolata sull'intervallo di psi specificato
    # Plot della Te TS vs Te ECE + retta x=y
    titlefig = f'{shot}_Tts_vs_Tece'
    
    fig06,ax06 = plt.subplots(nrows=1, sharex=True, num = titlefig)
    ax06.errorbar(p, q, xerr = err_eceM_R, yerr = err_xm12, 
                  marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5,label='ECE vs TS')   
    ax06.plot(x,y,'g--', lw=.8)
    ax06.set_xlim(left=4)
    ax06.set_ylim(bottom=4)
    ax06.legend()
    ax06.set_title(f'JPN {shot} - Te HRTS vs Te Ece-Michelson.{psi1}<PSI<{psi2} for {tlim1}<t<{tlim2} (s)')
    ax06.set_xlabel('Te Ece-Michelson (keV)')
    ax06.set_ylabel('Te HRTS (keV)')
    fig06.tight_layout()
    
    if d['savefigs'] == 1: 
       plt.savefig(d['mypath']+f'{shot}_Tts_vs_Tece_PSI.pdf',dpi=300)
            
    return 1

###############################################################################
# Calcolo le medie sull'intervallo di rho scelto, all'interno dell'intervallo temporale t1-t2
def mean_calc_rho(d, vars): 
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    rho1 = d['rho1']
    rho2 = d['rho2']
    eP = d['eP']
    tTs = vars['tTs']  
    errTs = vars['errTs']
    tEce = vars['tEce'] 
    time_ts = vars['time_ts']
    time_ece = vars['time_ece']
    rho_ts = vars['rho_ts']
    rho_ece = vars['rho_ece']
    ####################
    
    idt = (tTs.t >= tlim1) & (tTs.t <= tlim2) # Seleziono gli indici dei tempi di interesse
    temp_ts2 = tTs.v[:,idt]/1000 # in keV
    err_ts2 = errTs.v[:,idt]/1000  
    err_ts2[err_ts2>2] = 1 # metto a 1 keV l'errore sui punti dove diverge
    
    idt2 = (time_ts >= tlim1) & (time_ts <= tlim2)
    rho_ts2 = rho_ts[:,idt2] 
                 
    # Ciclo sulle Selezione intervallo rho e media tra rho1 e rho2: 
    # Per il HRTS
    
    temp_tsM_rho = []   # temp media in inttervallo rho
    err_tsM_rho = []   # errore medio in intervallo rho
    
    for i in range(0,(rho_ts2.shape[1])):
        mask = (rho_ts2[:,i] >= rho1) & (rho_ts2[:,i] <= rho2)    # Seleziono gli indici dei tempi di interesse 
        temp = temp_ts2[mask,i]                    # Valori di temperaatura nelle rho selezionate
        sig = err_ts2[mask,i]       # Errore sui singoli valori di temperatura: sigma
        ai = 1/sig**2            # Inverso del quadrato delle sigma
        somma = sum(ai)
        if somma==0:
            somma=1
        xm_ = sum(ai*temp)/somma   # Media dei valori di temperatura pesata sui singoli errori
        temp_tsM_rho.append(xm_)               # Valore delle temperature 'attese'/più probabili, nel tempo
        err_xm_ = np.sqrt(1/somma) 
        err_tsM_rho.append(err_xm_)
            
    ### Medie sull'intervallo di psi scelto      
    
    idt = (tEce.t >= tlim1) & (tEce.t <= tlim2) # Seleziono gli indici dei tempi di interesse
    temp_ece2 = tEce.v[:,idt]/1000  # in keV
    err_ece2 = temp_ece2*eP
    
    idt2 = (time_ece>=tlim1) & (time_ece <= tlim2)
    rho_ece2 = rho_ece[:,idt2]
    
    temp_eceM_rho = np.zeros(temp_ece2.shape[1])
    rho_eceM_rho = np.zeros(rho_ece2.shape[1])
    err_eceM_rho = np.zeros(temp_ece2.shape[1])
    
    # Colcolo delle medie nell'intervallo di rho tra rho1 e rho2 nell'intervallo di tempo scelto
    for i in range(0,(rho_ece2.shape[1])):
        mask1 = (rho_ece2[:,i] >= rho1) & (rho_ece2[:,i] <= rho2)    # Seleziono gli indici delle psi di interesse 
        temp_eceM_rho[i]  = np.mean(temp_ece2[mask1,i], axis=0)
        rho_eceM_rho[i] = np.mean(rho_ece2[mask1,i],axis=0)
        err_eceM_rho[i] = np.mean(err_ece2[mask1,i], axis=0)

    vars['temp_tsM_rho'] = temp_tsM_rho
    vars['err_tsM_rho'] = err_tsM_rho
    vars['temp_eceM_rho'] = temp_eceM_rho
    vars['err_eceM_rho'] = err_eceM_rho
    
    return 1 

###############################################################   

def fig_cfr_rho(d, vars): #shot, w, tlim1, tlim2, delta, psi1, psi2, xm, err_xm, temp_eceM, err_eceM):  
    shot = d['shot']
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    rho1 = d['rho1']
    rho2 = d['rho2']
    save = d['savefigs']
    temp_tsM_rho = vars['temp_tsM_rho']
    err_tsM_rho = vars['err_tsM_rho']
    temp_eceM_rho = vars['temp_eceM_rho']
    err_eceM_rho = vars['err_eceM_rho']
    ####################
    
    dimTs = len(temp_tsM_rho)
    dimEce = len(temp_eceM_rho)
    
    temp_eceM_rho_R  = signal.resample_poly(temp_eceM_rho, up=dimTs, down=dimEce) # posizioni ricampionate dei valori TeEce
    err_eceM_rho_R = signal.resample_poly(err_eceM_rho, up=dimTs, down=dimEce)    # posizioni ricampionate dei valori Err TeEce
 
    def retta(x):
        return x
    x = np.linspace(4,12,dimTs)   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    q = temp_tsM_rho             # Valori TeHRTS media-pesata calcolata sull'intervallo di psi specificato
    p = temp_eceM_rho_R   # Valori TeECE media calcolata sull'intervallo di psi specificato
    # Plot della Te TS vs Te ECE + retta x=y
    titlefig = f'{shot}_Tts_vs_Tece_rho'
    
    fig07,ax07 = plt.subplots(nrows=1, sharex=True, num = titlefig)
    ax07.errorbar(p, q, xerr = err_eceM_rho_R, yerr = err_tsM_rho, 
                  marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5,label='ECE vs TS')   
    ax07.plot(x,y,'g--', lw=.8)
    # ax07.set_xlim(left=4)
    # ax07.set_ylim(bottom=4)
    ax07.legend()
    ax07.set_title(f'JPN {shot} - Te HRTS vs Te Ece-Michelson. \n {rho1}<RHO<{rho2} for {tlim1}<t<{tlim2} (s)')
    ax07.set_xlabel('Te Ece-Michelson (keV)')
    ax07.set_ylabel('Te HRTS (keV)')
    fig07.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Tts_vs_Tece_RHO.pdf',dpi=300)
            
    return 1
    
 