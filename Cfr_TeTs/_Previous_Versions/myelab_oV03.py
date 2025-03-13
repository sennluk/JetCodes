#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:29:18 2024
@author: lsenni. Documentazione su file word: Note_Codice_Cfr_Te_ts
Ho cercato di uniformare i nomi delle variabili
I plot 'di servizio li faccio con plt.plot
I plot da salvare con fig e ax

Funzioni:
    magax: plot delle posizione dell'asse magnetico
    tprof(shot, rad, tlim1, tlim2, delta): restituisce l'andamenti temporale 
                       di Te per Ece e Ts per la coordinata radiale 'rad'
                       nell'intervallo temporale tra tlim1 e tlim2
    psicalc(shot,tim): restituisce andamento di PSI(R) e Te(PSI) 
                       per lo shoot 'shot'a dato tempo tim
    rhocalc                   
"""
import numpy as np
import my_flush
from ppfeg import ppfs
import matplotlib.pyplot as plt
from scipy import signal,interpolate

# Ricostriuscio la posizione dell'asse magnetico 
def magax(shot, w, tlim1, tlim2):
    plt.figure('Magnetic axis position over time')
    w.efit.rmag.plot()
    plt.axvline(x = tlim1,c='r',ls='--',lw=.5, label='t-lim1')
    plt.axvline(x = tlim2,c='r',ls='--',lw=.5, label = 't-lim2')
    plt.xlabel('time (sec)')
    plt.ylabel('R(m)')
    plt.title(f'JPN {shot} - Radial position of the magnetic axis from EFIT')
    
    return

##################################################

def tprof(shot, w, rad, tlim1, tlim2, delta, eP):
       
    tTs = w.hrts.te       # Chan Te HRTS  
    errTs = w.hrts.dte    # Chann Errors HRTS
    tEce = w.ecm1.prfl    # Chan El Temp ECE Michelson 
###############################################################################
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
    ax00.set_title(f'{shot} - Te trend over time at Rts={posTs} m and Rece={posEce}')
    ax00.set_ylabel('Te (keV)')  
    ax00.set_xlabel('Time (s)')
    
    return fig00, ax00
##################################################

def psicalc(shot, w):
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
          
    return psiKk1

##################################################
def psifig(shot, w, psiKk1, tlim, psi1, psi2):
    tTs = w.hrts.te       # Chan. Te HRTS  
    psiTs = w.hrts.psi    # Channel psi hrts
    tEce = w.ecm1.prfl 
    timeEce = tEce.t
    timeTs = tTs.t
    # Calcolo del profilo a dato tempo 

    idts = np.argmin(abs(timeTs - tlim))
    psiTs_s = np.absolute(psiTs.v[:,idts]) # Abs del profilo psi hrts al tempo t=tlim   
    
    ide= np.argmin(abs(timeEce - tlim))
    psiEce_s = psiKk1[:,ide]     # profilo psi ece al tempo t=tlim          

    tempTs_s = tTs.slice(t=tlim)
    tempEce_s = tEce.slice(t=tlim)
    rTs = psiTs.r
    rEce = tEce.r
    # Seleziono gli indici corrispondenti all'intervallo in psi:psi1-psi2
    idxPsiTs = ((psiTs_s >= psi1) & (psiTs_s <= psi2))
    idxPsiE = ((psiEce_s >= psi1) & (psiEce_s <= psi2))  # Indici dell'array di PsiEce
    
    fig01, ax01 = plt.subplots(nrows=1,  num = f'JPN {shot} PSI profiles')
    ax01.plot(rTs, psiTs_s, linewidth=0.5, color='green', label= f'{psi1}<psi hrts<{psi2}')
    ax01.plot(rTs[idxPsiTs], psiTs_s[idxPsiTs], linewidth=1.5, color='red', label='psi hrts')
    ax01.plot(rEce, psiEce_s, linewidth=0.5, color='blue', label='psi ece')
    ax01.plot(rEce[idxPsiE], psiEce_s[idxPsiE], linewidth=1.5, color='orange', label=f'{psi1}<psi ece<{psi2}')
    ax01.axhline(y = psi1, c='r',ls='--',lw=.4)
    ax01.axhline(y = psi2,c='r',ls='--',lw=.4)
    ax01.set_xlabel('R (m)')
    ax01.set_ylabel('PSI')
    ax01.legend()
    ax01.set_title(f'JPN {shot} PSI(R) for HRTS and ECE-KK1 at t={tlim} sec')

    fig02, ax02 = plt.subplots(nrows=1, sharex=True, num = f'JPN {shot} Te vs PSI')
    ax02.scatter(psiEce_s, tempEce_s.v/1000, label='Te ece - kk1')
    ax02.scatter(psiTs_s, tempTs_s.v/1000,label='Te hrts')
    ax02.set_xlabel('PSI')
    ax02.set_ylabel('Electron Temperature (keV)')
    ax02.legend()
    ax02.set_title(f'JPN {shot} Te(PSI) at t={tlim} sec')
    
    print('PSI1 ECE = ', psiEce_s[idxPsiE][0], ' - R1 ECE = ', rEce[idxPsiE][0])
    print('PSI2 ECE = ', psiEce_s[idxPsiE][-1], ' - R2 ECE = ', rEce[idxPsiE][-1])
    print('Number of PSI-ECE Values = ', psiEce_s[idxPsiE].size)
    
    print('PSI1 HRTS = ', psiTs_s[idxPsiTs][0], ' - R1 HRTS = ', rTs[idxPsiTs][0])
    print('PSI2 HRTS = ', psiTs_s[idxPsiTs][-1], ' - R2 HRTS = ', rTs[idxPsiTs][-1])
    print('Number of PSI-HRTS Values = ', psiTs_s[idxPsiTs].size)
    
    # return fig01, ax01, fig02, ax02 #, psi_ece, psi_ts, tempEce, tempTs 

##################################################

def rhocalc(shot, w, psiKk1, tlim):
    tTs = w.hrts.te       # Chan Te HRTS  
    psiTs = w.hrts.psi    # Channel psi hrts
    tEce = w.ecm1.prfl 
    
    timeTs = tTs.t
    timeEce = tEce.t
        
    rhoTs = np.zeros(tTs.v.shape)
    rhoEce = np.zeros(tEce.v.shape)
     
    for i,time in enumerate(timeTs):
        psi_th = psiTs.v[:,i]
        fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i]) # function to intertpolate the flux
        fl_int_hrts = fl_int(psi_th)
        rhoTs[:,i] = np.sqrt(fl_int_hrts)
        
    for i,time in enumerate(timeEce):
        psi = psiKk1[:,i]
        fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i]) # function to intertpolate the flux
        fl_int_ece = fl_int(psi)
        rhoEce[:,i] = np.sqrt(fl_int_ece)
        
        return rhoTs,rhoEce
##################################################        

def rhofig(shot, w, tlim, rho1, rho2, rhoTs, rhoEce): 
    tTs = w.hrts.te       # Chan Te HRTS  
    tEce = w.ecm1.prfl 
    
    timeTs = tTs.t
    timeEce = tEce.t
        
    # Calcolo del profilo a dato tempo - slices
    tempTs_s = tTs.slice(t=tlim)
    tempEce_s = tEce.slice(t=tlim)
        
    idts = np.argmin(abs(timeTs - tlim))
    ide = np.argmin(abs(timeEce - tlim))       
    
    plt.figure()
    plt.scatter(rhoTs[:,idts], tempTs_s.v/1000)
    plt.scatter(rhoEce[:,ide],tempEce_s.v/1000)
    plt.title(f'HRTS and Ece Te vs rho - Profile at t={tlim}')
    plt.xlabel('rho')    
    plt.ylabel('Te (keV)')
    
    # Seleziono gli indici corrispondenti all'intervallo in rho: rho1-rho2
    idts = np.argmin(abs(timeTs - tlim))
    rhoTs_s = rhoTs[:,idts] # Profilo rho hrts al tempo t=tlim   
    
    ide= np.argmin(abs(timeEce - tlim))
    rhoEce_s = rhoEce[:,ide]     # profilo rho ece al tempo t=tlim          

    tempTs_s = tTs.slice(t=tlim)
    tempEce_s = tEce.slice(t=tlim)
    rTs = tTs.r
    rEce = tEce.r
    
    idxRhoTs = ((rhoTs_s >= rho1) & (rhoTs_s <= rho2))
    idxRhoE = ((rhoEce_s >= rho1) & (rhoEce_s <= rho2))  # Indici dell'array di PsiEce
    
    fig01, ax01 = plt.subplots(nrows=1,  num = f'JPN {shot} Normalized RHO profiles')
    ax01.plot(rTs, rhoTs_s, linewidth=0.5, color='green', label='Rho Hrts LoS')
    ax01.plot(rTs[idxRhoTs], rhoTs_s[idxRhoTs], linewidth=1.5, color='red', label=f'{rho1}<rho hrts<{rho2}')
    ax01.plot(rEce, rhoEce_s, linewidth=0.5, color='blue', label='Rho Ece LoS')
    ax01.plot(rEce[idxRhoE], rhoEce_s[idxRhoE], linewidth=1.5, color='orange', label='{rho1}<Rho Ece<{rho2}')
    ax01.axhline(y = rho1, c='r',ls='--',lw=.4)
    ax01.axhline(y = rho2,c='r',ls='--',lw=.4)
    ax01.set_xlabel('R (m)')
    ax01.set_ylabel('Normalized RHO')
    ax01.legend()
    ax01.set_title(f'JPN {shot} RHO(R) for HRTS and ECE-KK1 at t={tlim} sec')

    # fig02, ax02 = plt.subplots(nrows=1, sharex=True, num = f'JPN {shot} Te vs PSI')
    # ax02.scatter(psiEce_s, tempEce_s.v/1000, label='Te ece - kk1')
    # ax02.scatter(psiTs_s, tempTs_s.v/1000,label='Te hrts')
    # ax02.set_xlabel('PSI')
    # ax02.set_ylabel('Electron Temperature (keV)')
    # ax02.legend()
    # ax02.set_title(f'JPN {shot} Te(PSI) at t={tlim} sec')
    
    # print('PSI1 ECE = ', psiEce_s[idxPsiE][0], ' - R1 ECE = ', rEce[idxPsiE][0])
    # print('PSI2 ECE = ', psiEce_s[idxPsiE][-1], ' - R2 ECE = ', rEce[idxPsiE][-1])
    # print('Number of PSI-ECE Values = ', psiEce_s[idxPsiE].size)
    
    # print('PSI1 HRTS = ', psiTs_s[idxPsiTs][0], ' - R1 HRTS = ', rTs[idxPsiTs][0])
    # print('PSI2 HRTS = ', psiTs_s[idxPsiTs][-1], ' - R2 HRTS = ', rTs[idxPsiTs][-1])
    # print('Number of PSI-HRTS Values = ', psiTs_s[idxPsiTs].size)
    
    
    
    
##################################################
# Calcolo delle medie nell'intervallo di psi sceltoi tra gli estremi psi1 e psi2

def meancalc(shot, w, tlim1, tlim2, delta, psi1, psi2, eP, psiKk1, rhoTs, rhoEce):
    tTs = w.hrts.te       # Chan Te HRTS
    psiTs = w.hrts.psi    # Channel psi hrts
    tEce = w.ecm1.prfl 
    
    timeEce = tEce.t
    
    idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ts = tTs.t[idt]
    temp_ts = tTs.v[:,idt]/1000 # in keV
    psi_ts = psiTs.v[:,idt]
    err_ts = w.hrts.dte.v[:,idt]/1000  
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
    temp_ece = w.ecm1.prfl.v[:,idt]/1000  # in keV
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
        
    return temp_tsM, err_tsM, xm, err_xm, temp_eceM, err_eceM
###############################################################    
def fig_psi_mean(shot, w, tlim1, tlim2, delta, psi1, psi2, eP, temp_tsM, err_tsM, xm, err_xm, temp_eceM, err_eceM):
    tTs = w.hrts.te       # Chan Te HRTS
    tEce = w.ecm1.prfl 
    
    idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ts = tTs.t[idt]
    
    idt = (tEce.t >= tlim1-delta) & (tEce.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ece = tEce.t[idt]
    # Plot in psi e rho
    
    # PLOT CON ERRORBARS e dato con media aritmetica e pesata
    fig03, ax03 = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Profile -PSI- Time trend - Errorbars')
    linew = 0.7
    ax03.lw = 0.5
    ax03.errorbar(time_ts,temp_tsM,  lw = linew, color = 'b', 
                 yerr= err_tsM, ecolor='g', elinewidth= 0.2, label='HRTS mean')    # arithmetic mean HRTS
    ax03.errorbar(time_ts,xm,lw = linew, color='g',
                  yerr = err_xm, ecolor='g', elinewidth=.3, label='HRTS w-mean') # xm is the weighted mean for HRTS
    ax03.errorbar(time_ece,temp_eceM, lw = linew, color='orange', 
                  yerr= err_eceM, ecolor='r', elinewidth= 0.2, label='ECE mean')      # arithmetic mean ECE Michelson
    #ax03.plot(time_ece,temp_eceM,lw = linew, color='orange', label='ECE mean' )             
    ax03.fill_between(time_ece,temp_eceM+err_eceM/2,temp_eceM-err_eceM/2,color = 'darkorange',alpha=0.3)
    ax03.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax03.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax03.legend(fontsize=8)
    ax03.set_title(f'JPN {shot} - Mean (Ece and Hrts) and w-Mean(Hrts)  for {psi1}<PSI<{psi2}')
    ax03.set_ylabel('Te (keV)') 
    ax03.set_xlabel('time (s)')
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
    
    return #fig03, ax03, fig05, ax05 # se metto gli outpu poi a volte me lo sovrascrive 
    
    #############################################
def fig_cfr_psi(shot, w, tlim1, tlim2, delta, psi1, psi2, xm, err_xm, temp_eceM, err_eceM):    
    tTs = w.hrts.te       # Chan Te HRTS
    tEce = w.ecm1.prfl 
    
    idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ts = tTs.t[idt]
    
    idt = (tEce.t >= tlim1-delta) & (tEce.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ece = tEce.t[idt]
    
    dim = time_ts.size
    temp_eceM_R  = signal.resample_poly(temp_eceM, up=150, down=402)
    err_eceM_R = signal.resample_poly(err_eceM, up=150, down=402)
 
    def retta(x):
        return x
    x = np.linspace(-1,13,dim)   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    q = xm
    p = temp_eceM_R
    # Plot della Te TS vs Te ECE + retta x=y
    
    fig06,ax06 = plt.subplots(nrows=1, sharex=True, num = f'{shot}_Tts_vs_Tece')
    ax06.errorbar(p,q,xerr = err_eceM_R, yerr = err_xm,marker='o', markersize=3,ecolor='g',linestyle='none', 
                 elinewidth=.5,label='ECE vs TS')   
    ax06.plot(x,y,'g--', lw=.8)
    ax06.set_xlim(left=4)
    ax06.set_ylim(bottom=4)
    ax06.legend()
    ax06.set_title(f'JPN {shot} - Te HRTS vs Te Ece-Michelson. {tlim1}<t<{tlim2} - {psi1}<PSI<{psi2}')
    ax06.set_xlabel('Te Ece-Michelson (keV)')
    ax06.set_ylabel('Te HRTS (keV)')
    # if save == 1: 
    #     plt.savefig(f'{shot}_Tts_vs_Tece_Err',dpi=600)
            
    return fig06, ax06 
