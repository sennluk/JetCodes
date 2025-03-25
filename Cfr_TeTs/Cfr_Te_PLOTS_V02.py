"""
Created on 28/01/2025 
@author: lsenni. 
please refer to word file: Note_Codice_Cfr_Te_ts
Based on the myelab_V12.py library file
Contains the functions to do plots

Functions:
    multiplot: fornische una multiplot con le caratteristiche principali dello shot
    magax: fornisce l'andamento delle posizione dell'ase magnetico secondo EFIT sul piano poloidale
    tprof: mostra l'andamento del profilo di Te per le due diagnostiche
    rhofig: 
Tolgo la delta dagli intervalli temporali nel multiplot- inserisco degli 'antispike' laddove servano
V02: elimino i PPFS dalle vars, e richiamo quindi direttamente i dati contenuti del dizionario
    

OSS: molte soluzioni nelle figure per l'articolo 'standardazing
Aggiungo Ratio e difference e relativi plot''

"""
# Branch new_plots - 24/03/2025
import numpy as np
from ppfeg import ppfs
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.signal import savgol_filter
######################## Multiplot
# Shot characterization 
def multiplot(d,vars): 
    shot = d['shot']      
    tlim1 = vars['ti'] # Tolgo i 40 secondi di norma del JET 
    tlim2 = vars['tf']
    rad = vars['rad']
    delta = d['delta']
    ##################
       
    w = ppfs(shot) 
    tlim1 = tlim1-40-delta # Tolgo i 40 secondi di norma del JET #-delta
    tlim2 = tlim2-40+delta # +1 #  Tolgo i 40 secondi di norma del JET  #+delta
    
    btor = w.magn.bvac   # toroidal magnetic field
    idt = (btor.t-40 >= tlim1) & (btor.t-40 <= tlim2) 
    tBtor = btor.t[idt]-40
    btor = abs(btor.v[0,idt])  # In MWatt  
    
    ipla = w.magn.ipla   # Total NBI Power
    idt = (ipla.t-40 >= tlim1) & (ipla.t-40 <= tlim2) 
    tIpla = ipla.t[idt]-40
    ipla = abs(ipla.v[0,idt]/10**6)  # In MAmps  
    
    pnbi = w.nbi.ptot   # Total NBI Power
    idt = (pnbi.t-40 >= tlim1) & (pnbi.t-40 <= tlim2) 
    tPnbi = pnbi.t[idt]-40
    pnbi = pnbi.v[0,idt]/10**6  # In MWatt   
    
    pbolo = w.bolo.topi     # Tot rad power (improved) - Bolometry
    idt = (pbolo.t-40 >= tlim1) & (pbolo.t-40 <= tlim2)
    tPrad = pbolo.t[idt]-40
    prad = pbolo.v[0,idt]/10**6
    prad[(prad>30) | (prad<0)] = 0
    
    picrh = w.icrh.ptot   # Total NBI Power
    idt = (picrh.t-40 >= tlim1) & (picrh.t-40 <= tlim2) 
    tPicrh = picrh.t[idt]-40
    picrh = picrh.v[0,idt]/10**6  # In MWatt  
    
    betan = w.efit.btnm      # Beta normalizaed - MHD
    idt = (betan.t-40 >= tlim1) & (betan.t-40 <= tlim2)
    tBetan = betan.t[idt]-40
    betan = betan.v[0,idt]
    
    tmax = w.hrtx.tmax      # max Te max hrts (hrtx channel)
    idt = (tmax.t-40 >= tlim1) & (tmax.t-40 <= tlim2)
    tTmax = tmax.t[idt]-40
    tmax = tmax.v[0,idt]/1000  # in keV
    
    ti = w.xcs.ti  # Ion temperature 
    idt = (ti.t-40 >= tlim1) & (ti.t-40 <= tlim2)
    tTi = ti.t[idt]-40
    ti = ti.v[0,idt]/1000  # in keV
    # ti = w.ks5
    
    dens = w.hrts.ne.slice(r=rad)  # density profile hrts at r=rad
    idt = (dens.t-40 >= tlim1) & (dens.t-40 <= tlim2)
    tNe = dens.t[idt]-40
    ne = dens.v[0,idt]/10**19
    
    qu = w.efit.qax       # Simulated q on axis
    idt = (qu.t-40 >= tlim1) & (qu.t-40 <= tlim2)
    tQu = qu.t[idt]-40
    qu = qu.v[0,idt]
    
    nr = w.tin.rnt # KN1 neutron rate
    idt = (nr.t-40 >= tlim1) & (nr.t-40 <= tlim2)
    tNr = nr.t[idt]-40
    nr = nr.v[0,idt]/10**16
    
    gas = w.gash.eler # KN1 neutron rate
    idt = (gas.t-40 >= tlim1) & (gas.t-40 <= tlim2)
    tGas = gas.t[idt]-40
    gas = gas.v[0,idt]/10**21
    
    ###############################
    linew = 0.5  # plot  lines dimension
    fonts = 6 # 5.5 for pdf   # size of the legend fonts
    fs = 8 # 6 for pdf       # size of the axes labels
    fst = 6     # size of the ticks
    
    fig, (ax1, ax2, ax3, ax4, ax5, ax6,  ax8) = plt.subplots(nrows=7, sharex=True, num=f'JPN {shot} Trends over time')
    plt.subplots_adjust(hspace=0)
    
    ax1.plot(tBtor,btor, lw = linew, label='$B_{tor}$')
    ax1.set_ylabel('T', fontsize= fs)
    ax1.yaxis.set_tick_params(labelsize=fst)
    ax1.legend(fontsize = fonts, loc="upper left").set_draggable(True)
    # ax1.set_title(f'Traces JPN {shot}') 
    # ax1.set_title(f'Traces over time')
    
    ax2.plot(tIpla, ipla, lw = linew, color='b', label= r'$I_{pla}$')
    ax2.set_ylabel('MA', fontsize= fs)
    ax2.yaxis.set_tick_params(labelsize=fst)
    ax2.legend(fontsize = fonts, loc="upper left").set_draggable(True) 
    
    ax3.plot(tPnbi,pnbi, lw = linew, label='$P_{NBI}$')
    ax3.plot(tPrad, prad, lw = linew, label='$P_{rad}$')
    ax3.plot(tPicrh, picrh, lw = linew, label='$P_{icrh}$')
    ax3.set_ylabel('MW', fontsize= fs)
    ax3.yaxis.set_tick_params(labelsize=fst)
    ax3.legend(fontsize = fonts, loc="upper left").set_draggable(True)
    
    ax4.plot(tBetan, betan, lw = linew, color='b', label= r'$\beta_N$')
    ax4.yaxis.set_tick_params(labelsize=fst)
    ax4.legend(fontsize = fonts, loc="upper left").set_draggable(True) 
    
    ax5.plot(tTmax, tmax, lw = linew, color='r', label='$T_{e}$ max')
    ax5.plot(tTi, ti, lw = linew, color='g', label='$T_{i}$ <$N_i$>26')
    ax5.set_ylabel('keV', fontsize= fs)
    ax5.yaxis.set_tick_params(labelsize=fst)
    ax5.legend(fontsize = fonts, loc="upper left").set_draggable(True)
    
    # posne = round(rad[0],2)
    ax6.plot(tNe, ne, lw = linew, color='darkblue', label=f'Density at R={rad[0]:.2f} m')
    ax6.set_ylabel('$10^{19} m^{-3}$', fontsize= fs)
    ax6.yaxis.set_tick_params(labelsize=fst)
    ax6.legend(fontsize = fonts, loc="upper left").set_draggable(True)
    
    ax8.plot(tNr, nr, lw = linew, color='darkblue', label='Neutron rate')
    ax8.set_ylabel('$10^{16} n/s$', fontsize = fs)
    ax8.yaxis.set_tick_params(labelsize=fst)
    ax8.set_xlabel('time(sec)', fontsize = fs)
    ax8.legend(fontsize = fonts, loc="upper left").set_draggable(True) 
    # plt.savefig('Fig_0_Traces.pdf', bbox_inches="tight", dpi=100)
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Multiplot.pdf',bbox_inches="tight", dpi=300)
    ## Non usati:

##################################################
# Time trend at R=rad: cenrtro plasma =  minimo di RHO
def tprof(d,vars):
    shot = d['shot']
    # delta = d['delta']
    tlim1 = vars['ti']-40 # -delta
    tlim2 = vars['tf']-40 # +delta
    tTs_v = vars['tTs_v']      # Chan Te HRTS  
    tTs_t = vars['tTs_t']-40 
    tTs_r = vars['tTs_r'] 
    errTs = vars['errTs']    # Chann Errors HRTS
    tEce_v = vars['tEce_v']    # Chan El Temp ECE Michelson 
    tEce_t = vars['tEce_t']-40 
    tEce_r = vars['tEce_r'] 
    errEce = vars['errEce']
    rad = vars['rad']
    ####################
    # Plot time trend at a specific position R --> To have a check
    
    indice = np.argmin(abs(tTs_r-rad)) # indice della posizione del max
    sTs = tTs_v[indice,:]           # Slice TS data
    sErrTs = errTs[indice,:]

    idt = (tTs_t >= tlim1) & (tTs_t <= tlim2) # Seleziono gli indici dei tempi di interesse
    timeTs = tTs_t[idt]
    tempTs = sTs[idt]/1000   # in keV
    errTs = sErrTs[idt]/1000  # in keV
    errTs[errTs>5] = 1       # Controllo che l'errore non sia troppo elevato--> errato
    print('TS position=',round(rad[0],4))
    
    indice = np.argmin(abs(tEce_r-rad)) # indice della posizione del max
    sEce = tEce_v[indice,:]
    sErrEce = errEce[indice,:]

    idt = (tEce_t >= tlim1) & (tEce_t <= tlim2) # Seleziono gli indici dei tempi di interesse
    timeEce = tEce_t[idt]
    tempEce = sEce[idt]/1000   # in keV
    errEce = sErrEce[idt]/1000  # in keV
    errEce[errEce>5] = 1       # Controllo che l'errore non sia troppo elevato--> errato
    posEce = tEce_r[indice]
    print('ECE position=',round(posEce,4))
    
    linew = 0.5
    posizione = round(rad[0],2) 
    fig00,(ax00) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Time trend at R')
    ax00.lw = 0.5
    # ax00.errorbar(timeTs, tempTs, color = 'royalblue', lw = linew, yerr = errTs, ecolor='g', elinewidth=.2, label='Tts')
    # ax00.errorbar(timeEce, tempEce, color = 'tomato', lw = linew, yerr = errEce, ecolor='r', elinewidth=.2, label='Tece')
    ax00.plot(timeTs, tempTs, color = 'royalblue', lw = linew, label='Te - HRTS')
    ax00.plot(timeEce, tempEce, color = 'tomato', lw = linew, label='Te - ECE')
    ax00.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax00.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax00.xaxis.set_major_locator(MultipleLocator(1))
    ax00.legend(fontsize=10, loc = "upper left").set_draggable(True)
    ax00.set_title(f'Te time trend (no averages) at R = {posizione:.2f}')
    ax00.set_ylabel('Te (keV)')  
    ax00.set_xlabel('Time (s)')
    fig00.tight_layout()
    # plt.savefig('Fig_0_Te_vs_Time.pdf', bbox_inches="tight", dpi=100)
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(t)_at_Rad.pdf',dpi=300)
                   
    return 1
##################################################        

##################################################    
# Ricostruisco la posizione dell'asse magnetico 
def magax(d,vars):
    shot = d['shot']
    tlim1 = vars['ti']
    tlim2 = vars['tf']
    pippo = vars['rmag']
    pluto = vars['zmag'] 
    paperino = vars['tmag']
    zKk1 = vars['zKk1'] 
    zTs_v = vars['zTs_v']
    zTs_r = vars['zTs_r']
    ####################
    idt = (paperino >= tlim1) & (paperino <= tlim2)
    
    # cmap='YlGnBu', zored: regola la sovrapposizione
    plt.figure(f'JPN {shot} Mag ax pos vs time')
    plt.plot(pippo[0,:],pluto[0,:],lw=0.7, color = 'blue', zorder=1)
    sc = plt.scatter(pippo[0,idt], pluto[0,idt], c=paperino[idt]-40, cmap='viridis', s=4, zorder=2) 
    plt.colorbar(sc, label="Time (s)")
    plt.axhline(y = zKk1, c='r',ls='--',lw=.4, label='ECE LoS')
    plt.plot(zTs_r, zTs_v, c='g',ls='--',lw=.4, label = 'HRTS LoS')
    plt.xlim(right=max(pippo[0,:]+0.05))
    plt.xlabel('R(m)')
    plt.ylabel('z(m)')
    plt.legend()
    plt.title(f'JPN {shot} - Position of the magnetic axis on the poloidal plane - from EFIT')
    # plt.savefig('Fig_04_Mag_Ax_Pos.pdf',bbox_inches="tight", dpi=300) # Per art
    
    zTs_t = np.linspace(paperino[0], paperino[-1], zTs_v.shape[0])
    w = ppfs(shot)
    plt.figure(f'JPN {shot} Vert Mag ax pos vs time')
    w.efit.zmag.plot()
    plt.axvline(x = tlim1,c='r',ls='--',lw=.5, label='ti')
    plt.axvline(x = tlim2,c='r',ls='--',lw=.5, label = 'tf')
    plt.axhline(y = zKk1, c='r',ls='--',lw=.4, label='KK1 LoS')
    plt.plot(zTs_t, zTs_v, c='g',ls='--',lw=.4, label = 'HRTS LoS')
    plt.xlabel('time (sec)')
    plt.ylabel('z(m)')
    plt.title(f'JPN {shot} - Poloidal position of the magnetic axis from EFIT')
    plt.legend()
    plt.tight_layout()
 
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Mag_Ax_Pos.pdf',format = "pdf", dpi=300, bbox_inches="tight")
              
    return 1
 
##################################################        
# Plot profili di Te in RHO ad un dato tempo: di dafault t=tlim : tempo un cui è max la Tmax
def rho_fig(d, vars): 
    shot = d['shot']
    # rho1 = d['rho1']
    # rho2 = d['rho2']
    tlim = vars['tlim']
    tTs_v = vars['tTs_v']       # Chan. Te HRTS 
    tTs_t = vars['tTs_t']       # Chan. Te HRTS 
    tTs_r = vars['tTs_r']
    tEce_v = vars['tEce_v']
    tEce_t = vars['tEce_t']
    tEce_r = vars['tEce_r']
    # psiTs_r = vars['psiTs_r'] 
    rhoTs = vars['rhoTs']
    rhoEce = vars['rhoEce']
    timeTs2 = vars['timeTs2']
    ranges = vars['ranges']
    ####################
    timeEce = tEce_t
    timeTs = tTs_t
    # Calcolo del profilo a dato tempo: seleziono l'oindice temporale
    idts = np.argmin(abs(timeTs - tlim))    
    ide= np.argmin(abs(timeEce - tlim))

    tempTs_s = tTs_v[:,idts]     # _s: slice at time..
    tempEce_s = tEce_v[:,ide]
    # rTs = psiTs_r
    rTs = tTs_r
    rEce = tEce_r
    ##################    
    
    # Devo selezionare l'indice corrispondente a tlim in time Ts2
    # per poter determinare il range di RHO a tale tempo
    indice = np.argmin(abs(timeTs2-tlim))
    rho1 = ranges[indice,0] # che sarebbe rhodown
    rho2 = ranges[indice,1] # che sarebbe rhoup
    
    tempo = round(tlim,3)        
    fig08, ax08 = plt.subplots(nrows=1, sharex=True, num = f'JPN {shot} Te vs RHO')
    ax08.scatter(rhoTs[:,idts], tempTs_s/1000, marker='o', lw = 0.7, facecolors='none', edgecolors='darkorange', label='HRTS')
    ax08.scatter(rhoEce[:,ide], tempEce_s/1000, marker='1', lw = 0.7, color = 'olive', label='ECE-KK1')
    ax08.axvline(x = rho1, c='r',ls='--',lw=.4)
    ax08.axvline(x = rho2, c='r',ls='--',lw=.4)
    # ax08.set_title(f'JPN. {shot} - HRTS and Ece Te vs rho - Profile at t={tlim}')
    ax08.set_title(f'HRTS and Ece Te profiles vs RHO - Profile at t={tempo:.2f} - Max Te')
    ax08.legend()
    ax08.set_xlabel('rho')    
    ax08.set_ylabel('Te (keV)')
    plt.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(rho).pdf',dpi=300)
            
    # Seleziono gli indici corrispondenti all'intervallo in rho: rho1-rho2
    rhoTs_s = rhoTs[:,idts] # Profilo rho hrts al tempo t=tlim   
    rhoEce_s = rhoEce[:,ide]     # profilo rho ece al tempo t=tlim          
    
    idxRhoTs = ((rhoTs_s >= rho1) & (rhoTs_s <= rho2))
    idxRhoE = ((rhoEce_s >= rho1) & (rhoEce_s <= rho2))  # Indici dell'array di PsiEce
    leftlim = min(min(rTs[idxRhoTs]), min(rEce[idxRhoE]))
    rightlim = max(max(rTs[idxRhoTs]), max(rEce[idxRhoE]))
    
    rhoi = round(rho1,2)
    rhof = round(rho2,2)

    fig03, (ax003, ax03) = plt.subplots(nrows=2, sharex=False, num = f'JPN {shot} RHO profiles')
    ax003.plot(rTs, rhoTs_s, linewidth=0.5, color='green', label=r'HRTS LoS($\rho$)')
    ax003.plot(rEce, rhoEce_s, linewidth=0.5, color='blue', label=r'ECE LoS($\rho$)')
    ax003.axhline(y = rho1, c='r',ls='--',lw=.4)
    ax003.axhline(y = rho2,c='r',ls='--',lw=.4)
    ax003.set_ylabel(r'$\rho_N$')
    ax003.legend()
    ax03.plot(rTs, rhoTs_s, '-o', ms=1, linewidth=0.5, color='green', label=r'HRTS LoS($\rho$)')
    ax03.plot(rTs[idxRhoTs], rhoTs_s[idxRhoTs], 'o', ms=3,  linewidth=1.5, color='red', label=fr'{rhoi}<$\rho_N$-hrts<{rhof}')
    ax03.plot(rEce, rhoEce_s, '-*', ms=1, linewidth=0.5, color='blue', label=r'ECE LoS($\rho$)')
    ax03.plot(rEce[idxRhoE], rhoEce_s[idxRhoE], 'o', ms=3, linewidth=1.5, color='black', label=fr'{rhoi}<$\rho_N$-ece<{rhof}')
    ax03.axhline(y = rho1, c='r',ls='--',lw=.4)
    ax03.axhline(y = rho2,c='r',ls='--',lw=.4)
    ax03.set_xlim(left = leftlim-0.1, right = rightlim+0.1)
    # ax03.set_xlim(2.95, right = 3.20)
    ax03.set_ylim(bottom = rho1-0.05, top = rho2 + 0.01)
    ax03.set_xlabel('R (m)')
    ax03.set_ylabel(r'$\rho_N$')
    ax03.legend() #fontsize=6, loc = "lower left")
    ax003.set_title(f'JPN {shot} RHO(R) for HRTS and ECE-KK1 at t={tlim} sec')
    fig03.tight_layout()
    # plt.savefig('Fig_07_2_LoSs_vs_Rho.pdf',bbox_inches="tight", dpi=100) # Per art
    
    if d['savefigs'] == 1: 
       plt.savefig(d['mypath']+f'{shot}_RHO(R)_at_t.pdf',dpi=300)
       
    # print('RHOi HRTS = ', round(rhoTs_s[idxRhoTs][0],3), ' - Ri HRTS = ', round(rTs[idxRhoTs][0],3))
    # print('RHOf HRTS = ', round(rhoTs_s[idxRhoTs][-1],3), ' - Rf HRTS = ', round(rTs[idxRhoTs][-1],3))
    # print('Number of RHO-HRTS Values = ', rhoTs_s[idxRhoTs].size)
    # print('Rho-hrts average over a lenght (cm): ', round((rTs[idxRhoTs][-1]-rTs[idxRhoTs][0])*100,3))
    # print('RHOi ECE = ', round(rhoEce_s[idxRhoE][0],3), ' - Ri ECE = ', round(rEce[idxRhoE][0],3))
    # print('RHOf ECE = ', round(rhoEce_s[idxRhoE][-1],3), ' - Rf ECE = ', round(rEce[idxRhoE][-1],3))
    # print('Number of RHO-ECE Values = ', rhoEce_s[idxRhoE].size)
    # print('Rho-ece average over a lenght (cm): ', round((rEce[idxRhoE][-1]-rEce[idxRhoE][0])*100,3))
       
    return fig03, ax03
##########################################
# Plot of the Temeprature time trends (ECE and HRTS) - smoothed -
# of the averages over the selected RHO ranges - With errorbars

def te_trends(d,vars):
    shot = d['shot']
    win_len = d['win_len']
    deg_pol = d['deg_pol']
    tlim1 = vars['ti']
    tlim2 = vars['tf']
    temp_tsM_rho = vars['temp_tsM_rho']
    err_tsM_rho = vars['err_tsM_rho']
    temp_eceM_rho = vars['temp_eceM_rho']
    err_eceM_rho = vars['err_eceM_rho']
    timeTs2 = vars['timeTs2']

    # PLOT cfr dati smootati CON ERRORBARS ece aritmetica vs ts pesata
    linew = 0.5
   
    fig00,(ax00) = plt.subplots(nrows=1, sharex=True, num=f'JPN {shot} - Time trends averaged')
    ax00.lw = linew
    ax00.errorbar(timeTs2, temp_tsM_rho, color = 'royalblue', lw = linew, yerr = err_tsM_rho, ecolor='c', elinewidth=.2, label='Tts')
    ax00.errorbar(timeTs2, temp_eceM_rho, color = 'tomato', lw = linew, yerr = err_eceM_rho, ecolor='r', elinewidth=.2, label='Tece')
    ax00.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax00.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax00.legend(fontsize=8)
    ax00.set_title('Te time trend around plasma center - average on RHO Ranges')
    ax00.set_ylabel('Te (keV)')  
    ax00.set_xlabel('Time (s)')
    fig00.tight_layout()
    
    # Same plot with a smooothing algorithm
    fig005, ax005 = plt.subplots(nrows=1, sharex=True, num= f'JPN {shot} Smooted - PSI- Time trend - Bands')
    linew = 0.7
    ax005.lw = 0.5
    
    time = timeTs2
    te_hrts = savgol_filter(temp_tsM_rho, win_len, deg_pol) # smooth media aritmetica
    err_hrts = savgol_filter(err_tsM_rho, win_len, deg_pol)   
    te_ece = savgol_filter(temp_eceM_rho, win_len, deg_pol) 
    err_ece = savgol_filter(err_eceM_rho, win_len, deg_pol) 
    
    ax005.plot(time, te_hrts, lw=linew, color='darkolivegreen',label='HRTS w-mean')
    ax005.fill_between(time, te_hrts+err_hrts, te_hrts-err_hrts, color='g',alpha=0.3) # xm is the weighted mean for HRTS
    ax005.plot(time,te_ece,lw = linew, color='orange', label='ECE mean' )             
    ax005.fill_between(time,te_ece+err_ece,te_ece-err_ece,color = 'darkorange',alpha=0.3)
    ax005.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax005.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax005.legend(fontsize=8)
    ax005.set_title(f'JPN {shot} -Smooted Ece-Mean and Smooted HRTS w-Mean trends vs time ')
    # ax005.set_title(f'Ece and HRTS Te for {psi1}<PSI<{psi2} vs time ')
    ax005.set_ylabel('Te (keV)') 
    ax005.set_xlabel('time (s)')
    fig005.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te_MEANS(psi)_BANDS_Smooted.pdf',dpi=300)

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
   
    fig, ax = plt.subplots(3, num = f'JPN {shot} Check plots ECE vs HRTS')
    ax[0].plot(timeTs2,timeEce22, label='cfr tempi')
    ax[0].legend()
    ax[0].set_title(f'JPN {shot} - Check Plots')
    ax[1].plot(timeTs2,ranges[:,1], label='rho-up')
    ax[1].plot(timeTs2,ranges[:,0], label='rho-down')
    ax[1].legend()
    ax[2].plot(timeTs2, temp_tsM_rho, label='TS averaged')
    ax[2].plot(timeEce22,temp_eceM_rho, label='ECE averaged')
    ax[2].legend()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Check plots ECE vs HRTS.pdf',dpi=300)
              
    left = min(np.nanmin(temp_eceM_rho), np.nanmin(temp_tsM_rho))
    right  = max(np.nanmax(temp_eceM_rho), np.nanmax(temp_tsM_rho)) 
    def retta(x):
        return x
    
    x = np.linspace(left-0.5,right+0.5,timeTs2.shape[0])   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    fig,ax=plt.subplots(1, num = f' JPN {shot} ECE vs HRTS') # 
    # ax.scatter(temp_eceM_rho, temp_tsM_rho,label='Rho Averaged ECE vs TS ') 
    ax.plot(x,y,'g--', lw=.8, label=r'$T_e$-ECE = $T_e$-HRTS')
    ax.errorbar(temp_tsM_rho, temp_eceM_rho, xerr = err_tsM_rho, yerr = err_eceM_rho, 
                  marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5) #, label='ECE vs TS')
    ax.set_title(f'JPN {shot} - Te Ece-Michelson vs Te HRTS  for {tlim1:.2f}<t<{tlim2:.2f} (s)') #:.2f per avere 2 cifre decimali
    ax.set_xlabel(r'$T_e$-HRTS (keV)')
    ax.set_ylabel(r'$T_e$-ECE (keV)')
    ax.legend()
    fig.tight_layout()
    # plt.savefig('Fig_08_Tece_vs_Tts.pdf', bbox_inches="tight", dpi=100) # Per art
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_ECE vs HRTS.pdf',dpi=300)
 
 #### Figura per articolo, cercando di rispettare le proporzioni
    
    import tkinter as tk
    root = tk.Tk()
    screen_width = root.winfo_screenwidth() / root.winfo_fpixels('1i')  # Larghezza in pollici
    screen_height = root.winfo_screenheight() / root.winfo_fpixels('1i')  # Altezza in pollici
    root.destroy()  # Chiudi la finestra Tkinter
    
    # Crea la figura con la stessa dimensione dello schermo
    fig, ax = plt.subplots(figsize=(screen_width, screen_height), dpi=100)  
    
    # Plotta i dati
    ax.plot(x, y, 'g--', lw=.8, label=r'$T_e$-ECE = $T_e$-HRTS')
    ax.errorbar(temp_tsM_rho, temp_eceM_rho, xerr=err_tsM_rho, yerr=err_eceM_rho, 
                marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5)
    fsize = 22
    # Label e legenda
    ax.tick_params(axis='both',labelsize=20)
    ax.set_xlabel(r'$T_e$-HRTS (keV)', fontsize = fsize)
    ax.set_ylabel(r'$T_e$-ECE (keV)', fontsize = fsize)
    ax.legend(fontsize = fsize, loc="upper left")
    fig.tight_layout()
    fig.canvas.draw()
    plt.savefig('Fig_08_2_Tece_vs_Tts.pdf', bbox_inches="tight", dpi=100)  # Usa lo stesso dpi per mantenere le proporzioni
    plt.show()
    
##########################################
# plot ratio vs time and difference vs time
def fig_rat_dist(d, vars): #shot, w, tlim1, tlim2, delta, psi1, psi2, xm, err_xm, temp_eceM, err_eceM):  
    shot = d['shot']
    ratio = vars['ratio']
    distance = vars['distance']
    err_dist = vars['err_dist']
    err_dist_perc = vars['err_dist_perc']
    time = vars['timeTs2']
    te_ece = vars['temp_eceM_rho']
    ####################         
    # fig,ax=plt.subplots(1, num = 'Ratio vs time') # 
    # ax.plot(time,ratio, lw=.8, label=r'$T_e$-ECE / $T_e$-HRTS')
    # ax.axhline(y=1, c='r', ls='--', lw = 0.6)
    # ax.set_title(f'JPN {shot} - Ratio vs time') 
    # ax.set_xlabel('time(s)')
    # ax.set_ylabel(r'$T_e$-ECE / $T_e$-TS')
    # ax.legend()
    # fig.tight_layout()
    
    # if d['savefigs'] == 1: 
    #     plt.savefig(d['mypath']+f'{shot}_Ratio_vs_time.pdf',dpi=300)
 
    # fig,ax=plt.subplots(1, num = r'Ratio vs T_e ECE') # 
    # ax.scatter(te_ece, ratio, marker='o', s=1, label=r'$T_e$-ECE / $T_e$-HRTS vs $T_e$-ECE')
    # ax.axhline(y=1, c='r', ls='--', lw = 0.6)
    # ax.set_title(f'JPN {shot} - Ratio vs time') 
    # ax.set_xlabel(r'$T_e$-ECE (keV)')
    # ax.set_ylabel(r'$T_e$-ECE / $T_e$-TS')
    # ax.legend()
    # fig.tight_layout()
    
    # if d['savefigs'] == 1: 
    #     plt.savefig(d['mypath']+f'{shot}_Ratio_vs_Te.pdf',dpi=300)
    ############################## 
    fig,ax=plt.subplots(1, num = f' JPN {shot} Difference vs time') # 
    ax.plot(time,distance, lw=.8, label=r'$T_e$-ECE - $T_e$-HRTS')
    ax.axhline(y=0, c='r', ls='--', lw = 0.6)
    ax.set_title(f'JPN {shot} - Difference vs time') 
    ax.set_xlabel('time(s)')
    ax.set_ylabel(r'$T_e$-ECE - $T_e$-TS')
    ax.legend()
    fig.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Diff_vs_time.pdf',dpi=300)
 
    fig,ax=plt.subplots(1, num = f' JPN {shot} Difference vs T_e ECE') # 
    ax.scatter(te_ece, distance, marker='o', s=1, label=r'$T_e$-ECE - $T_e$-HRTS vs $T_e$-ECE')
    ax.axhline(y=0, c='r', ls='--', lw = 0.6)
    ax.set_title(f'JPN {shot} - Difference vs time') 
    ax.set_xlabel(r'$T_e$-ECE (keV)')
    ax.set_ylabel(r'$T_e$-ECE - $T_e$-TS')
    ax.legend()
    fig.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Diff_vs_Te.pdf',dpi=300)
        
    fig,ax=plt.subplots(1, num = f' JPN {shot} Diff and err vs T_e ECE') # 
    ax.errorbar(te_ece, distance, yerr=err_dist, 
                marker='o', markersize=1, ecolor='g', linestyle='none', elinewidth=.5, label=r'Difference vs $T_e$-ECE') 
    #marker='o', s=1, )
    ax.axhline(y=0, c='r', ls='--', lw = 0.6)
    ax.set_title(f'JPN {shot} - Difference (av. values) with errorbars vs Te-Ece') 
    ax.set_xlabel(r'$T_e$-ECE (keV)')
    ax.set_ylabel(r'$T_e$-ECE - $T_e$-TS')
    ax.legend()
    fig.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Diff_err_vs_Te.pdf',dpi=300)    
        
        
    fig,ax=plt.subplots(1, num = f' JPN {shot} Diff/T_ece and err vs T_e ECE') # 
    ax.errorbar(te_ece, distance/te_ece, yerr=err_dist_perc, 
                marker='o', markersize=1, ecolor='g', linestyle='none', elinewidth=.5, label=r'Difference vs $T_e$-ECE') 
    #marker='o', s=1, )
    ax.axhline(y=0, c='r', ls='--', lw = 0.6)
    ax.set_title(f'JPN {shot} - Difference (av. values)/Te-ece with errorbars vs Te-Ece') 
    ax.set_xlabel(r'$T_e$-ECE (keV)')
    ax.set_ylabel(r'$T_e$-ECE - $T_e$-TS')
    ax.legend()
    fig.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Diff_Perc_err_vs_Te.pdf',dpi=300)        