"""
Created on 28/01/2025 
@author: lsenni. 
please refer to word file: Note_Codice_Cfr_Te_ts
Based on the myelab_V12.py library file
Contains the functions to do plots

Functions:
    multiplot: fornische una multiplot con le caratteristiche principali dello shot
    magax: fornisce l'andamento delle posizione dell'ase magnetico secondo EFIT sul piano poloidale
    tprof: mostra l'andamento del profilo di Te per le due ddiagnostiche
    rhofig: 

    '

"""
import numpy as np
# import my_flush
from ppfeg import ppfs
# from ppfeg import ppfs
import matplotlib.pyplot as plt
# from scipy import signal,interpolate
from scipy.signal import savgol_filter
# import my_manage_file as mym
# import time

######################## Multiplot
# Shot characterization 
def multiplot(d,vars): 
     
    shot = d['shot']      # 103117   # 99950, 99971
    # tlim1 = d['tlim1']    #43     # limite inferiore selezione tempi - in secondi
    # tlim2 = d['tlim2']    #51.4 
    tlim1 = vars['ti']
    tlim2 = vars['tf']
    delta = d['delta']
    rad = d['rad']        #3.00 
       
    w = ppfs(shot) 
    tlim1 = tlim1-delta
    tlim2 = tlim2  #+delta
    
    btor = w.magn.bvac   # toroidal magnetic field
    idt = (btor.t >= tlim1) & (btor.t <= tlim2) 
    tBtor = btor.t[idt]
    btor = abs(btor.v[0,idt])  # In MWatt  
    
    ipla = w.magn.ipla   # Total NBI Power
    idt = (ipla.t >= tlim1) & (ipla.t <= tlim2) 
    tIpla = ipla.t[idt]
    ipla = abs(ipla.v[0,idt]/10**6)  # In MAmps  
    
    pnbi = w.nbi.ptot   # Total NBI Power
    idt = (pnbi.t >= tlim1) & (pnbi.t <= tlim2) 
    tPnbi = pnbi.t[idt]
    pnbi = pnbi.v[0,idt]/10**6  # In MWatt   
    
    
    pbolo = w.bolo.topi     # Tot rad power (improved) - Bolometry
    idt = (pbolo.t >= tlim1) & (pbolo.t <= tlim2)
    tPrad = pbolo.t[idt]
    prad = pbolo.v[0,idt]/10**6
    prad[(prad>30) | (prad<0)] = 0
    
    
    picrh = w.icrh.ptot   # Total NBI Power
    idt = (picrh.t >= tlim1) & (picrh.t <= tlim2) 
    tPicrh = picrh.t[idt]
    picrh = picrh.v[0,idt]/10**6  # In MWatt  
    
    betan = w.efit.btnm      # Beta normalizaed - MHD
    idt = (betan.t >= tlim1) & (betan.t <= tlim2)
    tBetan = betan.t[idt]
    betan = betan.v[0,idt]
    
    tmax = w.hrtx.tmax      # max Te max hrts (hrtx channel)
    idt = (tmax.t >= tlim1) & (tmax.t <= tlim2)
    tTmax = tmax.t[idt]
    tmax = tmax.v[0,idt]/1000  # in keV
    
    ti = w.xcs.ti  # Ion temperature 
    idt = (ti.t >= tlim1) & (ti.t <= tlim2)
    tTi = ti.t[idt]
    ti = ti.v[0,idt]/1000  # in keV
    # ti = w.ks5
    
    dens = w.hrts.ne.slice(r=rad)  # density profile hrts at r=rad
    idt = (dens.t >= tlim1) & (dens.t <= tlim2)
    tNe = dens.t[idt]
    ne = dens.v[0,idt]/10**19
    
    qu = w.efit.qax       # Simulated q on axis
    idt = (qu.t >= tlim1) & (qu.t <= tlim2)
    tQu = qu.t[idt]
    qu = qu.v[0,idt]
    
    nr = w.tin.rnt # KN1 neutron rate
    idt = (nr.t >= tlim1) & (nr.t <= tlim2)
    tNr = nr.t[idt]
    nr = nr.v[0,idt]/10**16
    
    gas = w.gash.eler # KN1 neutron rate
    idt = (gas.t >= tlim1) & (gas.t <= tlim2)
    tGas = gas.t[idt]
    gas = gas.v[0,idt]/10**21
    
    ###############################
    linew = 0.5  # plot  lines dimension
    fonts = 8 # 5.5 for pdf   # size of the legend fonts
    fs = 8 # 6 for pdf       # size of the axes labels
    fst = 6     # size of the ticks
    
    fig, (ax1, ax2, ax3, ax4, ax5, ax6,  ax8) = plt.subplots(nrows=7, sharex=True, num=f'JPN {shot} Trends over time')
    plt.subplots_adjust(hspace=0)
    
    ax1.plot(tBtor,btor, lw = linew, label='$B_{tor}$')
    ax1.set_ylabel('T', fontsize= fs)
    ax1.yaxis.set_tick_params(labelsize=fst)
    ax1.legend(fontsize = fonts, loc="upper right")
    ax1.set_title(f'Traces JPN {shot}')
    # ax1.set_title(f'Traces over time')
    
    ax2.plot(tIpla, ipla, lw = linew, color='b', label= r'$I_{pla}$')
    ax2.set_ylabel('MA', fontsize= fs)
    ax2.yaxis.set_tick_params(labelsize=fst)
    ax2.legend(fontsize = fonts, loc="upper right")
    
    ax3.plot(tPnbi,pnbi, lw = linew, label='$P_{NBI}$')
    ax3.plot(tPrad, prad, lw = linew, label='$P_{rad}$')
    ax3.plot(tPicrh, picrh, lw = linew, label='$P_{icrh}$')
    ax3.set_ylabel('MW', fontsize= fs)
    ax3.yaxis.set_tick_params(labelsize=fst)
    ax3.legend(fontsize = fonts, loc="upper right")
    
    ax4.plot(tBetan, betan, lw = linew, color='b', label= r'$\beta_N$')
    ax4.yaxis.set_tick_params(labelsize=fst)
    ax4.legend(fontsize = fonts, loc="upper right")
    
    ax5.plot(tTmax, tmax, lw = linew, color='r', label='$T_{e}$ max')
    ax5.plot(tTi, ti, lw = linew, color='g', label='$T_{i}$ <$N_i$>26')
    ax5.set_ylabel('keV', fontsize= fs)
    ax5.yaxis.set_tick_params(labelsize=fst)
    ax5.legend(fontsize = fonts, loc="upper right")
    
    ax6.plot(tNe, ne, lw = linew, color='darkblue', label=f'Density at R={rad} m')
    ax6.set_ylabel('$10^{19} m^{-3}$', fontsize= fs)
    ax6.yaxis.set_tick_params(labelsize=fst)
    ax6.legend(fontsize = fonts, loc="upper right")
    
    ax8.plot(tNr, nr, lw = linew, color='darkblue', label='Neutron rate')
    ax8.set_ylabel('$10^{16} n/s$', fontsize = fs)
    ax8.yaxis.set_tick_params(labelsize=fst)
    ax8.set_xlabel('time(sec)', fontsize = fs)
    ax8.legend(fontsize = fonts, loc="upper right")
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Multiplot.pdf',dpi=300)
    ## Non usati:
  
##################################################    
    
# Ricostruisco la posizione dell'asse magnetico 
def magax(d,vars):
    shot = d['shot']
    # tlim1 = d['tlim1'] 
    # tlim2 = d['tlim2']
    tlim1 = vars['ti']
    tlim2 = vars['tf']
    w = vars['w']
    ####################

    pippo = w.efit.rmag
    pluto = w.efit.zmag
    zKk1 = w.ecm1.antp.v[1,0]
    zTs = w.hrts.z
    idt = (pippo.t >= tlim1) & (pippo.t <= tlim2)
    
    plt.figure('Mag ax pos vs time')
    plt.plot(pippo.v[0,:],pluto.v[0,:],lw=0.7)
    plt.plot(pippo.v[0,idt],pluto.v[0,idt],lw=1, label='t1<t<t2' )
    plt.axhline(y = zKk1, c='r',ls='--',lw=.4, label='KK1 LoS')
    plt.plot(zTs.r,zTs.v, c='g',ls='--',lw=.4, label = 'HRTS LoS')
    plt.xlim(right=max(pippo.v[0,:]+0.05))
    plt.xlabel('R(m)')
    plt.ylabel('z(m)')
    plt.legend()
    plt.title(f'JPN {shot} - Position of the magnetic axis on the poloidal plane - from EFIT')
    
    zTs = w.hrts.z.v[1,0]
    plt.figure('Vert Mag ax pos vs time')
    w.efit.zmag.plot()
    plt.axvline(x = tlim1,c='r',ls='--',lw=.5, label='ti')
    plt.axvline(x = tlim2,c='r',ls='--',lw=.5, label = 'tf')
    plt.axhline(y = zKk1, c='r',ls='--',lw=.4, label='KK1 LoS')
    plt.axhline(y = zTs, c='g',ls='--',lw=.4, label = 'HRTS LoS')
    plt.xlabel('time (sec)')
    plt.ylabel('z(m)')
    plt.title(f'JPN {shot} - Poloidal position of the magnetic axis from EFIT')
    plt.legend()
    plt.tight_layout()

    
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Mag_Ax_Pos.pdf',dpi=300)
              
    return 1
##################################################

##################################################
# Time trend at R=rad[']
def tprof(d,vars):
    shot = d['shot']
    # tlim1 = d['tlim1'] 
    # tlim2 = d['tlim2']
    tlim1 = vars['ti']
    tlim2 = vars['tf']
    # rad = d['rad']
    delta = d['delta']
    eP = d['eP']
    tTs = vars['tTs']      # Chan Te HRTS  
    errTs = vars['errTs']    # Chann Errors HRTS
    tEce = vars['tEce']    # Chan El Temp ECE Michelson 
    rhoTs = vars['rhoTs']
    ####################
    # Plot time trend at a specific position R --> To have a check
    
    idt = (tTs.t >= tlim1) & (tTs.t <= tlim2) 
    rho_ts = rhoTs[:,idt]
    minimo = np.min(np.nanmin(rho_ts))  # Minimo escudendo eventuali NaN
    indices = np.where(rho_ts == minimo)  # Come centro plasma prendo la posizione del minimo di Rho TS
    rad_ts = tTs.r[indices[0]]
    rad = rad_ts
    
    sTs = tTs.slice(r=rad)       # time slice for TS-HRTS --> the nearest R to rad
    sErrTs = errTs.slice(r=rad)
    
    idt = (sTs.t >= tlim1-delta) & (sTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    timeTs = sTs.t[idt]
    tempTs = sTs.v[0,idt]/1000   # in keV
    errTs = sErrTs.v[0,idt]/1000  # in keV
    errTs[errTs>5] = 1       # Controllo che l'errore non sia troppo elevato--> errato
    posTs = sTs.r
    print('TS position=',round(posTs[0],4))
    
    sEce = tEce.slice(r=rad)    # time slice for TS-HRTS ò the nearest R to rad
    idt = (sEce.t >= tlim1-delta) & (sEce.t <= tlim2+delta)
    timeEce = sEce.t[idt]
    tempEce = sEce.v[0,idt]/1000
    errEce = tempEce*eP   # eP = relative errorassigned to the Ece Measurement
    posEce = sEce.r
    print('Ece position = ', round(posEce[0],4))
    
    linew = 0.5
   
    fig00,(ax00) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Time trend at R')
    ax00.lw = 0.5
    ax00.errorbar(timeTs, tempTs, color = 'royalblue', lw = linew, yerr = errTs, ecolor='c', elinewidth=.2, label='Tts')
    ax00.errorbar(timeEce, tempEce, color = 'tomato', lw = linew, yerr = errEce, ecolor='r', elinewidth=.2, label='Tece')
    ax00.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax00.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax00.legend(fontsize=8)
    ax00.set_title(f'Te time trend (no averages) at R = {rad}')
    ax00.set_ylabel('Te (keV)')  
    ax00.set_xlabel('Time (s)')
    fig00.tight_layout()
  
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(t)_at_Rad.pdf',dpi=300)
        
            
    return 1
##################################################

def psifig(d,vars):
    shot = d['shot']
    # psi1 = d['psi1']
    # psi2 = d['psi2']
    tlim = vars['tlim']
    tTs = vars['tTs']       # Chan. Te HRTS  
    tEce = vars['tEce']
    psiTs = vars['psiTs']    # Channel psi hrts
    psiKk1 = vars['psiKk1']
    psi1 = vars['psi1']
    psi2 = vars['psi2']
   
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
    
    fig01,(ax001, ax01) = plt.subplots(nrows=2, sharex = False, num = 'PSI LoS profiles')
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
            
        
    fig02, ax02 = plt.subplots(nrows=1, sharex=True, num = 'Te vs PSI')
    ax02.scatter(psiTs_s, tempTs_s.v/1000, marker='o', lw = 0.7, facecolors='none', edgecolors='darkorange', label='Te hrts') #s= 5, 
    ax02.scatter(psiEce_s, tempEce_s.v/1000, marker='1', lw = 0.7, color = 'olive' ,label='Te ece - kk1') # s= 10,
    ax02.axvline(x = psi1,c='r',ls='--',lw=.5)
    ax02.axvline(x = psi2,c='r',ls='--',lw=.5)
    ax02.set_xlabel('PSI')
    ax02.set_ylabel('Electron Temperature (keV)')
    ax02.legend()
    ax02.set_title(f'JPN {shot} Te(PSI) profile at t={tlim} sec')
    fig02.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(psi)_at_t.pdf',dpi=300)
      
            # Trova l'ultimo indice di True

    # Trova l'ultimo indice di True nella prima sequenza di True
    ultimo_true_index = -1  # Valore di default se non trovati True
    for i, valore in enumerate(idxPsiTs):
        if valore:  # Se il valore è True
            ultimo_true_index = i  # Aggiorna l'indice
        elif ultimo_true_index != -1:  # Se trovi un False dopo aver trovato True
            break  # Esci dal ciclo se hai già trovato un True
    
    print("Ultimo indice di True nella prima sequenza:", ultimo_true_index)
      
    print('PSIi HRTS = ', round(psiTs_s[idxPsiTs][0],3), ' - Ri HRTS = ', round(rTs[idxPsiTs][0],3))
    print('PSIf HRTS = ', round(psiTs_s[idxPsiTs][-1],3), ' - Rf HRTS = ', round(rTs[idxPsiTs][-1],3))
    print('Number of PSI-HRTS Values = ', psiTs_s[idxPsiTs].size)
    print('Psi-hrts average over a length (cm): ', round((rTs[idxPsiTs][-1]-rTs[idxPsiTs][0])*100,3))
    print('PSIi ECE = ', round(psiEce_s[idxPsiE][0],3), ' - Ri ECE = ', round(rEce[idxPsiE][0],3))
    print('PSIf ECE = ', round(psiEce_s[idxPsiE][-1],3), ' - Rf ECE = ', round(rEce[idxPsiE][-1],3))
    print('Number of PSI-ECE Values = ', psiEce_s[idxPsiE].size)
    print('Psi-ece average over a length (cm): ', round((rEce[idxPsiE][-1]-rEce[idxPsiE][0])*100,3))

    return 1 


##################################################        

def rhofig(d, vars): 
    shot = d['shot']
    # rho1 = d['rho1']
    # rho2 = d['rho2']
    tlim = vars['tlim']
    tTs = vars['tTs']  
    tEce = vars['tEce'] 
    rhoTs = vars['rhoTs']
    rhoEce = vars['rhoEce']
    rho1 = vars['rho1']
    rho2 = vars['rho2']
    ####################
        
    # Calcolo del profilo a dato tempo - slices
    tempTs_s = tTs.slice(t=tlim)
    tempEce_s = tEce.slice(t=tlim)
        
    idts = np.argmin(abs(tTs.t - tlim))
    ide = np.argmin(abs(tEce.t - tlim))
           
    fig08, ax08 = plt.subplots(nrows=1, sharex=True, num = 'Te vs RHO')
    ax08.scatter(rhoTs[:,idts], tempTs_s.v/1000, marker='o', lw = 0.7, facecolors='none', edgecolors='darkorange', label='HRTS')
    ax08.scatter(rhoEce[:,ide], tempEce_s.v/1000, marker='1', lw = 0.7, color = 'olive', label='ECE-KK1')
    ax08.axvline(x = rho1, c='r',ls='--',lw=.4)
    ax08.axvline(x = rho2, c='r',ls='--',lw=.4)
    # ax08.set_title(f'JPN. {shot} - HRTS and Ece Te vs rho - Profile at t={tlim}')
    ax08.set_title(f'HRTS and Ece Te vs rho - Profile at t={tlim}')
    ax08.legend()
    ax08.set_xlabel('rho')    
    ax08.set_ylabel('Te (keV)')
    plt.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(rho).pdf',dpi=300)
            
    # Seleziono gli indici corrispondenti all'intervallo in rho: rho1-rho2
    idts = np.argmin(abs(tTs.t - tlim))
    rhoTs_s = rhoTs[:,idts] # Profilo rho hrts al tempo t=tlim   
    
    ide = np.argmin(abs(tEce.t - tlim))
    rhoEce_s = rhoEce[:,ide]     # profilo rho ece al tempo t=tlim          

    rTs = tTs.r
    rEce = tEce.r
    
    idxRhoTs = ((rhoTs_s >= rho1) & (rhoTs_s <= rho2))
    idxRhoE = ((rhoEce_s >= rho1) & (rhoEce_s <= rho2))  # Indici dell'array di PsiEce
    leftlim = min(min(rTs[idxRhoTs]), min(rEce[idxRhoE]))
    rightlim = max(max(rTs[idxRhoTs]), max(rEce[idxRhoE]))
    
    fig03, (ax003, ax03) = plt.subplots(nrows=2, sharex=False, num = 'RHO profiles')
    ax003.plot(rTs, rhoTs_s, linewidth=0.5, color='green', label='Rho Hrts LoS')
    ax003.plot(rEce, rhoEce_s, linewidth=0.5, color='blue', label='Rho Ece LoS')
    ax003.axhline(y = rho1, c='r',ls='--',lw=.4)
    ax003.axhline(y = rho2,c='r',ls='--',lw=.4)
    ax003.set_ylabel('Normalized RHO')
    ax003.legend()
    ax03.plot(rTs, rhoTs_s, '-o', ms=1, linewidth=0.5, color='green', label='Rho Hrts LoS')
    ax03.plot(rTs[idxRhoTs], rhoTs_s[idxRhoTs], '-o', ms=3,  linewidth=1.5, color='red', label=f'{rho1}<rho hrts<{rho2}')
    ax03.plot(rEce, rhoEce_s, '-*', ms=1, linewidth=0.5, color='blue', label='Rho Ece LoS')
    ax03.plot(rEce[idxRhoE], rhoEce_s[idxRhoE], '-*', ms=3, linewidth=1.5, color='orange', label=f'{rho1}<rho ece<{rho2}')
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
              
    print('RHOi HRTS = ', round(rhoTs_s[idxRhoTs][0],3), ' - Ri HRTS = ', round(rTs[idxRhoTs][0],3))
    print('RHOf HRTS = ', round(rhoTs_s[idxRhoTs][-1],3), ' - Rf HRTS = ', round(rTs[idxRhoTs][-1],3))
    print('Number of RHO-HRTS Values = ', rhoTs_s[idxRhoTs].size)
    print('Rho-hrts average over a lenght (cm): ', round((rTs[idxRhoTs][-1]-rTs[idxRhoTs][0])*100,3))
    print('RHOi ECE = ', round(rhoEce_s[idxRhoE][0],3), ' - Ri ECE = ', round(rEce[idxRhoE][0],3))
    print('RHOf ECE = ', round(rhoEce_s[idxRhoE][-1],3), ' - Rf ECE = ', round(rEce[idxRhoE][-1],3))
    print('Number of RHO-ECE Values = ', rhoEce_s[idxRhoE].size)
    print('Rho-ece average over a lenght (cm): ', round((rEce[idxRhoE][-1]-rEce[idxRhoE][0])*100,3))
        
    ################# Check plot for RHO computatios
    # rhoEce_s   # profilo rho ece al tempo t=tlim
    # rho mio: rhoMio
    rhoTs_s = rhoTs[:,idts]  # Profilo rho hrts al tempo t=tlim 
    selected_time_mio = tTs.t[idts]
    print('selected_time_mio:',selected_time_mio)
    
    # Rho Flush: rhoPPF
    # Rho flush idts: rhoFlush
    rhoPpf = vars['w'].hrts.rho
    idx = np.argmin(abs(rhoPpf.t - tlim))
    rhoPpf_s = rhoPpf.v[:,idx]
    rhoFlush = np.abs(rhoPpf_s)
    selected_time_slice = rhoPpf.t[idx]
    print('selected_time_slice:',selected_time_slice)
    
    # Plot di confronto tra i diversi modi di calcolare RHo: coord di flusso toroidale
    plt.figure()
    plt.plot(rTs, rhoTs_s, label='RHO-tor norm')
    plt.plot(rTs, rhoFlush, label='Norm Min Rad from Flush')
    plt.xlabel('R(m)')
    plt.ylabel('RHO')
    plt.legend()
    plt.title(f'JPN {shot} - Check plot: RHO-TORn HRTS vs Normalized Minor Radius from FLUSH(HRTS chann)')  
    if d['savefigs'] == 1: 
       plt.savefig(d['mypath']+f'{shot}_RHO-torN vs NomrMinRad.pdf',dpi=300)
       
    return fig03, ax03
##########################################
# plot for the direct comparison of Te ECE vs Te HRTS based oint he Ts time frame
def fig_cfr_rho(d, vars): #shot, w, tlim1, tlim2, delta, psi1, psi2, xm, err_xm, temp_eceM, err_eceM):  
    shot = d['shot']
    tlim1 = vars['ti']
    tlim2 = vars['tf']
    temp_tsM_rho = vars['temp_tsM_rho']
    err_tsM_rho = vars['err_tsM_rho']
    temp_eceM_rho = vars['temp_eceM_rho']
    err_eceM_rho = vars['err_eceM_rho']
    timeTs2 = vars['timeTs2']
    timeEce22 = vars['timeEce22']
    ranges = vars ['ranges']
    ####################
   
    fig, ax = plt.subplots(3, num = 'Check plots ECE vs HRTS')
    ax[0].plot(timeTs2,timeEce22, label='cfr tempi')
    ax[0].legend()
    ax[0].set_title(f'JPN {shot} - Check Plots')
    ax[1].plot(timeTs2,ranges[:,1], label='rho-up')
    ax[1].plot(timeTs2,ranges[:,0], label='rho-down')
    ax[1].legend()
    ax[2].plot(timeTs2, temp_tsM_rho, label='TS averaged')
    ax[2].plot(timeEce22,temp_eceM_rho, label='ECE averaged')
    ax[2].legend()
    
    left = min(np.nanmin(temp_eceM_rho), np.nanmin(temp_tsM_rho))
    right  = max(np.nanmax(temp_eceM_rho), np.nanmax(temp_tsM_rho)) 
    def retta(x):
        return x
    
    x = np.linspace(left-0.5,right+0.5,timeTs2.shape[0])   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    fig,ax=plt.subplots(1, num = 'ECE vs HRTS')
    # ax.scatter(temp_eceM_rho, temp_tsM_rho,label='Rho Averaged ECE vs TS ') 
    ax.plot(x,y,'g--', lw=.8)
    ax.errorbar(temp_tsM_rho, temp_eceM_rho, xerr = err_tsM_rho, yerr = err_eceM_rho, 
                  marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5, label='ECE vs TS')
    ax.set_title(f'JPN {shot} - Te Ece-Michelson vs Te HRTS  for {tlim1:.2f}<t<{tlim2:.2f} (s)') #:.2f per avere 2 cifre decimali
    ax.set_xlabel('Te HRTS (keV)')
    ax.set_ylabel('Te Ece-Michelson (keV)')
    ax.legend()
 
 
 
 
 