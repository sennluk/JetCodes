#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:29:18 2024
@author: lsenni
Funzioni:
    magax: plot delle posizione dell'asse magnetico
    tprof(shot, rad, tlim1, tlim2, delta): restituisce l'andamenti temporale 
                       di Te per Ece e Ts per la coordinata radiale 'rad'
                       nell'intervallo temporale tra tlim1 e tlim2
    psicalc(shot,tim): restituisce andamento di PSI(R) e Te(PSI) 
                       per lo shoot 'shot'a dato tempo tim
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

def psirhocalc(shot, w, tlim, psi1, psi2):
    # Prove di calcolo e plots

    tTs = w.hrts.te       # Chan Te HRTS  
    psiTs = w.hrts.psi    # Channel psi hrts

    tEce = w.ecm1.prfl 
    zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista (in metri), 
    # la seconda è l'intersezione con l'asse y. 'Position of antenna aperture centre (z,majR) '
    time_ece = tEce.t
    time_ts = tTs.t

    psiKk1 = np.zeros(tEce.v.shape)
    rho_ece_ = np.zeros(tEce.v.shape)
    rho_ts_ = np.zeros(tTs.v.shape)

    for i,time in enumerate(time_ece):
        r, te = tEce.get(t=time)  # Profilo a T=time
        z = np.full_like(r, zKk1) # Retta linea di vista parallela asse r
        ts, ier = my_flush.flushinit(15, shot, time) # Definisco il tipo di equilibrio da usare
        # e Prendo il tempo più vicino a quello desiderato
        psi, _ = my_flush.Flush_getFlux(r*100, z*100) # Le coordinate vanne messe in cm!
        psiKk1[:,i] = psi
            
        fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i]) # function to intertpolate the flux
        fl_int_ece = fl_int(psi)
        rho_ece_[:,i] = np.sqrt(fl_int_ece)

    for i,time in enumerate(time_ts):
        psi_th = psiTs.v[:,i]
        fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i]) # function to intertpolate the flux
        fl_int_hrts = fl_int(psi_th)
        rho_ts_[:,i] = np.sqrt(fl_int_hrts)

    # Profilo temperatura in rho al tempo di indice 350 - Check plot   
    plt.figure()
    plt.scatter(rho_ece_[:,350],tEce.v[:,350])    
        
    plt.figure()
    plt.scatter(rho_ts_[:,100],tTs.v[:,100])  

    # Calcolo del profilo a dato tempo 

    rEce = r
    idx= np.argmin(abs(time_ece - tlim))
    psi_ece = psiKk1[:,idx]     # profilo psi ece al tempo t=tlim
    rho_ece = rho_ece_[:,idx]
    
    idx = np.argmin(abs(time_ts-tlim))
    psi_ts = np.absolute(psiTs.v[:,idx]) # Abs del profilo psi hrts al tempo t=tlim 
    rho_ts = rho_ts_[:,idx]              

    temp_ece = tEce.slice(t=tlim)
    temp_ts = tTs.slice(t=tlim)
    # temp_ts = tempTs.v[idx]
    rad_ts = psiTs.r
    # Seleziono gli indicii corrispondenti all'intervallo in psi:psi1-psi2
    idxPsiE = ((psi_ece >= psi1) & (psi_ece <= psi2))  # Indici dell'array di PsiEce
    idxPsiTs = ((psi_ts >= psi1) & (psi_ts <= psi2))

    fig01, ax01 = plt.subplots(nrows=1,  num = f'JPN {shot} PSI profiles')
    ax01.plot(rEce,psi_ece, linewidth=0.5, color='blue', label='psi ece')
    ax01.plot(rEce[idxPsiE],psi_ece[idxPsiE], linewidth=1.5, color='orange', label='psi ece')
    ax01.plot(rad_ts,psi_ts, linewidth=0.5, color='green', label='psi hrts')
    ax01.plot(rad_ts[idxPsiTs],psi_ts[idxPsiTs], linewidth=1.5, color='red', label='psi hrts')
    ax01.set_xlabel('R (m)')
    ax01.set_ylabel('PSI')
    ax01.legend()
    ax01.set_title(f'JPN {shot} PSI(R) at t={tlim} sec')

    fig02, ax02 = plt.subplots(nrows=1, sharex=True, num = f'JPN {shot} Te vs PSI')
    ax02.scatter(psi_ece, temp_ece.v/1000, label='Te ece - kk1')
    ax02.scatter(psi_ts,temp_ts.v/1000,label='Te hrts')
    ax02.set_xlabel('PSI')
    ax02.set_ylabel('Electron Temperature (keV)')
    ax02.legend()
    ax02.set_title(f'JPN {shot} Te(PSI) at t={tlim} sec')
    
    print('PSI1 ECE = ', psi_ece[idxPsiE][0], ' - R1 ECE = ', rEce[idxPsiE][0])
    print('PSI2 ECE = ', psi_ece[idxPsiE][-1], ' - R2 ECE = ', rEce[idxPsiE][-1])
    print('Number of PSI-ECE Values = ', psi_ece[idxPsiE].size)
    
    print('PSI1 HRTS = ', psi_ts[idxPsiTs][0], ' - R1 HRTS = ', rad_ts[idxPsiTs][0])
    print('PSI2 HRTS = ', psi_ece[idxPsiE][-1], ' - R2 HRTS = ', rad_ts[idxPsiTs][-1])
    print('Number of PSI-HRTS Values = ', psi_ts[idxPsiTs].size)
    
    return fig01, ax01, fig02, ax02 #, psi_ece, psi_ts, tempEce, tempTs 

##################################################
# Calcolo delle medie nell'intervallo di psi sceltoi tra gli estremi psi1 e psi2

def meancalc(shot, w, tlim1, tlim2, delta, psi1, psi2, eP):
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
    errTs_[errTs_>2] = 1 # metto a 1 keV l'errore sui punti dove diverge
                
    # Ciclo sulle Selezione intervallo psi e media tra psi1 e psi2: 
    # Per il HRTS
    
    tempTsM = np.zeros(tempTs_.shape[1])  
    psiTsM_ = np.zeros(psiTs_.shape[1])
    errTsM = np.zeros(tempTs_.shape[1])
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
    errEce = tempEce*eP
    
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
    plt.title('Selected PSI mean values over time - {psi1}<PSI<{psi2}')
    # plt.ylim(psi1,psi2)    
    plt.legend()
    
    ####################### Plot in psi e rho
        
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
    ax03.set_title(f'JPN {shot} - Mean profile for {psi1}<PSI<{psi2}')
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
    ax05.set_title(f'JPN {shot} - Mean Profile for {psi1}<PSI<{psi2}')
    ax05.set_ylabel('Te (keV)') 
    ax05.set_xlabel('time (s)')
    
        
    #############################################
    
    dim = timeTs.size
    tempEceM_  = signal.resample_poly(tempEceM, up=150, down=402)
    errEceM_ = signal.resample_poly(errEceM, up=150, down=402)
 
    def retta(x):
        return x
    x = np.linspace(-1,13,dim)   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    q = xm
    p = tempEceM_
    # Plot della Te TS vs Te ECE + retta x=y
    
    fig06,ax06 = plt.subplots(nrows=1, sharex=True, num = f'{shot}_Tts_vs_Tece')
    ax06.errorbar(p,q,xerr = errEceM_, yerr = errXm,marker='o', markersize=3,ecolor='g',linestyle='none', 
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
            
    return fig03, ax03, fig05, ax05, fig06, ax06 
