"""
Created on Thu Mar 14 10:29:18 2024
@author: lsenni. Documentazione su file word: Note_Codice_Cfr_Te_ts

Introduco il calcolo di rho solamente nell'intervallo t1-Delta - t2 + Delta 
V07: cambio il calcolo del rho considerando la selezione rpecisa dei temnpi ai quali viene calcolato.   
V08: introduco la stima automatica dei range di psi e rho 
che vengono poin richiamati nelle medie: obiettivo è quello di relizzare le medie sugli stessi intervalli 
in psi e in rho    
I calcolati sono in vars, i preimpostati in d
"""
import numpy as np
import my_flush
from ppfeg import ppfs
# from ppfeg import ppfs
import matplotlib.pyplot as plt
from scipy import signal,interpolate
from scipy.signal import savgol_filter
import my_manage_file as mym
import time

# DTE3 shot list: 104990,991,994,995,999   
# 104520,521,522,523,524,526
##################################################
# Deinisco il dizionario delle variabili

def dictvar(d):
    shot = d['shot']
    timestr = time.strftime("%Y%m%d-%H%M%S")
    
    w = ppfs(shot)   
    
    if d['savefigs'] == 1: 
        mym.create_folder(d)
        mym.save_param(d)
        # mym.save_print_output_to_file(d,vars)
    else:
        print('You are not saving the plots')
     
    vars = {
        'tlim' : (d['tlim1']+d['tlim2'])/2,          # tempo medio al quale vengono calcolati i profili per plot di controllo
        'w' : w,                           # canale dati del ppfeg con tutti i dati relativi allo sparo in oggetto
        'psi1' : d['psi1'],    # inizializzo i valori di psi1,2 e rho1,2
        'psi2' : d['psi2'],
        'rho1' : d['rho1'],
        'rho2' : d['rho2'],
        'tTs' : w.hrts.te,                 # Chann Te HRTS  
        'errTs' : w.hrts.dte,              # Chann Errors HRTS.dte
        'tEce' : w.ecm1.prfl,              # Chann Te Ece-Kk1
        'errEce' : (w.ecm1.prfl)*d['eP'],  # Error associated to the Kk1 meas: ab 2%
        'psiTs' : w.hrts.psi,              # Chann PSI HRTS 
        'psiKk1' : None,
        'rhoKk1' : None,
        'time_ts' : None,
        'time_ece' : None,
        'rhoTs' : None,
        'rhoEce' : None, 
        'temp_tsM' : None, 
        'err_tsM': None, 
        'xm': None, 
        'err_xm' : None, 
        'temp_eceM' : None, 
        'err_eceM' : None, 
        'xm12': None, 
        'err_xm12' : None, 
        'temp_eceM12' : None, 
        'err_eceM12' : None, 
        'temp_tsM_rho' : None,
        'err_tsM_rho' : None,
        'temp_eceM_rho' : None,
        'err_eceM_rho' : None} 
    return vars

##################################################
# Define the time interval of the analysis for the single shot
# based on the temperature values higher than Tref

def tdef(d,vars):
    shot = d['shot']
    Tref = d['Tref']
    rad = d['rad']
    w = vars['w']
  
    ##############
    data = w.hrts.te.get(r=rad)
    timelim1 = data[0][np.where(data[1]/1000>Tref)][0]
    timelim2 = data[0][np.where(data[1]/1000>Tref)][-1]
    
    # temp = data[1]/1000
    # range = temp[np.where(temp>Tref)]
    # ind1 = range[0]
    # ind2 = range[-1]
    # timelim1 = temp[0,ind1]
    # timelim2 = temp[0,ind2]
    
    d['tlim1'] = round(timelim1,2) 
    d['tlim2'] = round(timelim2,2)
    
    
##################################################    
    
# Ricostruisco la posizione dell'asse magnetico 
def magax(d,vars):
    shot = d['shot']
    tlim1 = d['tlim1'] 
    tlim2 = d['tlim2']
    w = vars['w']
    ####################
    plt.figure('Rad Mag ax pos vs time')
    w.efit.rmag.plot()
    plt.axvline(x = tlim1,c='r',ls='--',lw=.5, label='t-lim1')
    plt.axvline(x = tlim2,c='r',ls='--',lw=.5, label = 't-lim2')
    plt.xlabel('time (sec)')
    plt.ylabel('R(m)')
    plt.title(f'JPN {shot} - Radial position of the magnetic axis from EFIT')
    plt.tight_layout()
    
    plt.figure('Vert Mag ax pos vs time')
    w.efit.zmag.plot()
    plt.axvline(x = tlim1,c='r',ls='--',lw=.5, label='t-lim1')
    plt.axvline(x = tlim2,c='r',ls='--',lw=.5, label = 't-lim2')
    plt.xlabel('time (sec)')
    plt.ylabel('z(m)')
    plt.title(f'JPN {shot} - Poloidal position of the magnetic axis from EFIT')
    plt.tight_layout()

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
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Mag_Ax_Pos.pdf',dpi=300)
              
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

def rhocalc(d, vars): 
    # tlim1 = d['tlim1']
    # tlim2 = d['tlim2']
    # delta = d['delta']
    w= vars['w']
    psiKk1 = vars['psiKk1']
    tTs = vars['tTs']       # Chan. Te HRTS  
    psiTs = vars['psiTs']    # Channel psi hrts
    tEce = vars['tEce']
    ####################
    
    torFlux = w.efit.ftor
    rFlux = w.efit.ftor.r  # Posizioni dei valori di flusso
    vFlux = w.efit.ftor.v  # valori di flusso
    
    rhoTs = np.zeros(tTs.v.shape)
    rhoEce = np.zeros(psiKk1.shape)
         
    for i,time in enumerate(tTs.t):
        psiTs_ = np.abs(psiTs.v[:,i])
        index = np.argmin(abs(torFlux.t - time))
        vFlux_ = vFlux[:,index]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1]) # function to interp the flux
        fl_int_hrts = fl_int(psiTs_, extrapolate=False) # Interpolation without extrapolating
        rhoTs[:,i] = np.sqrt(fl_int_hrts)
        
    for i,time in enumerate(tEce.t):
        psi = psiKk1[:,i]
        index = np.argmin(abs(torFlux.t - time))
        vFlux_ = vFlux[:,index]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1])  # function to intertp the flux
        fl_int_ece = fl_int(psi, extrapolate=False)
        rhoEce[:,i] = np.sqrt(fl_int_ece)
    
    vars['rhoTs'] = rhoTs
    vars['rhoEce'] = rhoEce
        
    return  1
##################################################
# Definisco i range di valori di PSI e RHO su cui fare le analisi:
# Prendo gli andamenti di psi per HRTS e ECE.
# Vedo qual'è piu basso e lo prendo come psi1 (tolgo qualcosa in più)
# prendo la posizione del minimo, aggiungo un intervallo(interv = 5 cm),
# e prendo il valore della psi della diagnostiche che non era con il minimo
# in corrispondenza della posizione del minimo +-interv   
def def_range(d,vars):

    shot = d['shot']
    interv = d['interv']
    tlim1 = d['tlim1'] 
    tlim2 = d['tlim2']
    tlim = vars['tlim']
    tTs = vars['tTs']
    tEce = vars['tEce']    # Chan Electron Temp ECE Michelson 
    psiTs = vars['psiTs']    # Channel psi hrts
    psiKk1 = vars['psiKk1']
    rhoTs = vars['rhoTs']
    rhoEce = vars['rhoEce']
    #########
    # interv = 0.1 # intervallo sul quale fare la media/2 (m)
    fix = 0.0001 # termine aggiuntivo per sicurezza valori psi e rho
    time = tlim
    psi_ts_ = psiTs.v[:,np.argmin(abs(tTs.t-time))] # valori di psi hrts all'istante scelto(_ finale per 'provvisorio')
    mp_ts = np.nanmin(abs(psi_ts_))                 # min coord psi-hrts
    pos_mp = psiTs.r[np.where(psi_ts_ == mp_ts)] # posizione del minimo di psi hrts in R
    psi_ece_ = psiKk1[:,np.argmin(abs(tEce.t-time))] # valori di psi ece all'istante scelto(_ finale per 'provvisorio')
    mp_ece = np.nanmin(abs(psi_ece_))                 # min coord psi-ece
    pos_mpe = tEce.r[np.where(psi_ece_ == mp_ece)] # posizione del minimo di psi ece in R
     
    if mp_ts<mp_ece:           # se il minimo della rho-hrts è più piccolo del minimo della rho-ece
        psi1_ = mp_ts - fix   # prendo come psi1 il min di psi-hrts a cui sottraggo un piccolo valore
        pos = pos_mp+interv     # posizione del minimo più interv (10 centimetri?)
        idx = np.argmin(abs(tEce.r-pos))
        psi2_ = psi_ece_[idx]+fix
        
        psi1 = round(max(psi1_, 0),4)
        psi2 = round(psi2_,4)
       
    else:
        psi1_ = mp_ece - fix   # prendo come psi1 le psi-ece a cui sottraggo un piccolo valore
        pos = pos_mpe+interv     # posizione del minimo più interv (10 centimetri?)
        idx = np.argmin(abs(tTs.r-pos))
        psi2_ = psi_ts_[idx]+fix
        
        psi1 = round(max(psi1_, 0),4)
        psi2 = round(psi2_,4)
            
    rho_ts_ = rhoTs[:,np.argmin(abs(tTs.t-time))]  # valore di rho hrts all'istante scelto
    mr_ts = np.nanmin(rho_ts_)               # min coord rho hrts
    pos_mr = tTs.r[np.where(rho_ts_ == mr_ts)]
    rho_ece_ = rhoEce[:,np.argmin(abs(tEce.t-time))]  # valore di rho ece all'istante scelto
    mr_ece = np.nanmin(rho_ece_)# min coord rho ece
    pos_mre = tEce.r[np.where(rho_ece_ == mr_ece)]
    
    if mr_ts<mr_ece:           # se il minimo della rho-hrts è più piccolo del minimo della rho-ece
        rho1_ = mr_ts - fix   # prendo come psi1 il min di psi-hrts a cui sottraggo un piccolo valore
        pos = pos_mre+interv     # posizione del minimo più interv (5 centimetri)
        idx = np.argmin(abs(tEce.r-pos))
        rho2_ = rho_ece_[idx]+fix
        
        rho1 = round(max(rho1_, 0),4)
        rho2 = round(rho2_,4)
       
    else:
        rho1_ = mr_ece - fix   # prendo come psi1 le psi-ece a cui sottraggo un piccolo valore
        pos = pos_mr+interv     # posizione del minimo più interv (5 centimetri)
        idx = np.argmin(abs(tTs.r-pos))
        rho2_ = rho_ts_[idx]+fix
        
        rho1 = round(max(rho1_, 0),4)
        rho2 = round(rho2_,4)
       
    print('Computed psi1 = ', psi1, 'in pos =', pos_mp)
    print('Computed psi2 = ', psi2, 'Approx range = ', 2*interv*100, 'cm')
    print('Computed rho1 = ', rho1, 'in pos =', pos_mr)
    print('Computed rho2 = ', rho2, 'Approx range = ', 2*interv*100, 'cm')
    
    
    vars['psi1'] = psi1
    vars['psi2'] = psi2
    vars['rho1'] = rho1
    vars['rho2'] = rho2
    
    return  1

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
    rhoTs = vars['rhoTs']
    ####################
    # Plot time trend at a specific position R --> To have a check
    
    idt = (tTs.t >= tlim1) & (tTs.t <= tlim2) 
    rho_ts = rhoTs[:,idt]
    minimo = np.min(np.nanmin(rho_ts))  # Minimo escudendo eventuali NaN
    indices = np.where(rho_ts == minimo)
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
   
    # fig00,(ax00) = plt.subplots(nrows=1, sharex=True, num=f'{shot} - Time trend at R')
    # ax00.lw = 0.5
    # ax00.errorbar(timeTs, tempTs, color = 'royalblue', lw = linew, yerr = errTs, ecolor='c', elinewidth=.2, label='Tts')
    # ax00.errorbar(timeEce, tempEce, color = 'tomato', lw = linew, yerr = errEce, ecolor='r', elinewidth=.2, label='Tece')
    # ax00.axvline(x = tlim1,c='r',ls='--',lw=.5)
    # ax00.axvline(x = tlim2,c='r',ls='--',lw=.5)
    # ax00.legend(fontsize=8)
    # ax00.set_title('Te time trend around plasma center - Raw data')
    # ax00.set_ylabel('Te (keV)')  
    # ax00.set_xlabel('Time (s)')
    # fig00.tight_layout()
    
# Figura per seminario:    
    fig00,(ax00) = plt.subplots(nrows=1, sharex=True, num='Time trend at R')
    ax00.lw = 0.5
    ax00.plot(timeTs, tempTs, color = 'forestgreen', lw = linew, label='Tts')
    ax00.errorbar(timeEce, tempEce, color = 'darkorange', lw = linew, label='Tece')
    ax00.legend(fontsize=10)
    # ax00.set_title('T$_e$ time trend around plasma center - Raw data - No errorbars', fontsize = 10)
    ax00.set_title(f'JPN{shot} T$_e$ time trend at R = {rad}- Raw data - No errorbars', fontsize = 10)
    ax00.set_ylabel(r'Te (keV)')  
    ax00.set_xlabel('Time (s)')
    # plt.savefig(d['mypath']+f'{shot}_Te(t)_at_Rad_X_SEMINARIO.pdf',dpi=300)
 
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(t)_at_Rad_X_SEMINARIO.pdf',dpi=300)
        
            
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
    
    fig01,(ax001, ax01) = plt.subplots(nrows=2, sharex = False, num = 'PSI profiles')
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
    ax02.set_xlabel('PSI')
    ax02.set_ylabel('Electron Temperature (keV)')
    ax02.legend()
    ax02.set_title(f'JPN {shot} Te(PSI) profile at t={tlim} sec')
    fig02.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te(psi)_at_t.pdf',dpi=300)
            
    print('PSIi HRTS = ', round(psiTs_s[idxPsiTs][0],3), ' - Ri HRTS = ', round(rTs[idxPsiTs][0],3))
    print('PSIf HRTS = ', round(psiTs_s[idxPsiTs][-1],3), ' - Rf HRTS = ', round(rTs[idxPsiTs][-1],3))
    print('Number of PSI-HRTS Values = ', psiTs_s[idxPsiTs].size)
    print('Psi-hrts average over a lenght (cm): ', round((rTs[idxPsiTs][-1]-rTs[idxPsiTs][0])*100,3))
    print('PSIi ECE = ', round(psiEce_s[idxPsiE][0],3), ' - Ri ECE = ', round(rEce[idxPsiE][0],3))
    print('PSIf ECE = ', round(psiEce_s[idxPsiE][-1],3), ' - Rf ECE = ', round(rEce[idxPsiE][-1],3))
    print('Number of PSI-ECE Values = ', psiEce_s[idxPsiE].size)
    print('Psi-ece average over a lenght (cm): ', round((rEce[idxPsiE][-1]-rEce[idxPsiE][0])*100,3))

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
    
    ide= np.argmin(abs(tEce.t - tlim))
    rhoEce_s = rhoEce[:,ide]     # profilo rho ece al tempo t=tlim          

    rTs = tTs.r
    rEce = tEce.r
    
    idxRhoTs = ((rhoTs_s >= rho1) & (rhoTs_s <= rho2))
    idxRhoE = ((rhoEce_s >= rho1) & (rhoEce_s <= rho2))  # Indici dell'array di PsiEce
    leftlim = min(min(rTs[idxRhoTs]), min(rEce[idxRhoE]))
    rightlim = max(max(rTs[idxRhoTs]), max(rEce[idxRhoE]))
    
    fig03, (ax003, ax03) = plt.subplots(nrows=2, sharex=False, num = 'Norm RHO profiles')
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
       
    return 1

##################################################        
# Utilizzo questa funzione per fare i  plot di controllo con il main: Eq_Coords_V00.py
def rhofig2(d, vars): 
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
                 
    # Seleziono gli indici corrispondenti all'intervallo in rho: rho1-rho2
    rhoTs_s = rhoTs[:,idts] # Profilo rho hrts al tempo t=tlim   
    rhoEce_s = rhoEce[:,ide]     # profilo rho ece al tempo t=tlim          

    rTs = tTs.r
    rEce = tEce.r
    
    idxRhoTs = ((rhoTs_s >= rho1) & (rhoTs_s <= rho2))
    idxRhoE = ((rhoEce_s >= rho1) & (rhoEce_s <= rho2))  # Indici dell'array di PsiEce
    leftlim = min(min(rTs[idxRhoTs]), min(rEce[idxRhoE]))
    rightlim = max(max(rTs[idxRhoTs]), max(rEce[idxRhoE]))
    
    fig03, (ax003, ax03) = plt.subplots(nrows=2, sharex=False, num = 'Norm RHO profiles')
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
    plt.plot(rTs, rhoTs_s, linewidth=0.5, color='green', label='Rho Hrts LoS')
    plt.plot(rEce, rhoEce_s, linewidth=0.5, color='blue', label='Rho Ece LoS')
    plt.plot(rTs, rhoFlush, label='Norm Min Rad from Flush')
    plt.xlabel('R(m)')
    plt.ylabel('RHO')
    plt.legend()
    plt.title(f'JPN {shot} - Check plot: RHO-TORn HRTS vs Normalized Minor Radius from FLUSH(HRTS chann)')  
    if d['savefigs'] == 1: 
       plt.savefig(d['mypath']+f'{shot}_RHO-torN vs NomrMinRad.pdf',dpi=300)
       
    return 1
    
    
##################################################
# Calcolo delle medie nell'intervallo di psi scelto tra gli estremi psi1 e psi2

def meancalc(d, vars): 
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    delta = d['delta']
    # psi1 = d['psi1']
    # psi2 = d['psi2']
    eP = d['eP']

    tTs = vars['tTs']  
    errTs = vars['errTs']
    tEce = vars['tEce'] 
    psiTs = vars['psiTs']
    psiKk1 = vars['psiKk1']
    psi1 = vars['psi1']
    psi2 = vars['psi2']

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
    # psi1 = d['psi1']
    # psi2 = d['psi2']
    win_len = d['win_len']
    deg_pol = d['deg_pol']
    tTs = vars['tTs']  
    tEce = vars['tEce'] 
    temp_tsM = vars['temp_tsM']
    err_tsM = vars['err_tsM']
    xm = vars['xm']
    err_xm = vars['err_xm']
    temp_eceM = vars['temp_eceM']
    err_eceM = vars['err_eceM']
    psi1 = vars['psi1']
    psi2 = vars['psi2']
    ####################
    
    idt = (tTs.t >= tlim1-delta) & (tTs.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ts = tTs.t[idt]
    
    idt = (tEce.t >= tlim1-delta) & (tEce.t <= tlim2+delta) # Seleziono gli indici dei tempi di interesse
    time_ece = tEce.t[idt]
    
    ##  Plot in psi e rho
    
    temp_tsM_sm = savgol_filter(temp_tsM, win_len, deg_pol) # smooth media aritmetica
    xm_sm = savgol_filter(xm, win_len, deg_pol)             # smooth media pesata
    err_xm_sm = savgol_filter(err_xm, win_len, deg_pol)
    temp_eceM_sm = savgol_filter(temp_eceM, win_len,deg_pol)
    err_eceM_sm = savgol_filter(err_eceM, win_len, deg_pol)
    
    # PLOT CON ERRORBARS e dato con media aritmetica e pesata
    fig04, ax04 = plt.subplots(nrows=1, sharex=True, num='PSI- Time trend - Errorbars')
    linew = 0.7
    ax04.lw = 0.5
    ax04.errorbar(time_ts,temp_tsM,  lw = linew, color = 'b', 
                 yerr= err_tsM, ecolor='g', elinewidth= 0.2, label='HRTS mean')    # arithmetic mean HRTS
    ax04.plot(time_ts,temp_tsM_sm, color='r')
    ax04.plot(time_ece, temp_eceM_sm)
    ax04.errorbar(time_ts,xm,lw = linew, color='g',
                  yerr = err_xm, ecolor='g', elinewidth=.3, label='HRTS w-mean') # xm is the weighted mean for HRTS
    ax04.errorbar(time_ece,temp_eceM, lw = linew, color='orange', 
                  yerr= err_eceM, ecolor='r', elinewidth= 0.2, label='ECE mean')      # arithmetic mean ECE Michelson            
    ax04.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax04.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax04.legend(fontsize=8)
    ax04.set_title(f'JPN {shot} - Mean (Ece and Hrts) and w-Mean(Hrts) for {psi1}<PSI<{psi2}')
    ax04.set_ylabel('Te (keV)') 
    ax04.set_xlabel('time (s)')
    fig04.tight_layout()
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te_MEANS(psi)_Comparison.pdf',dpi=300)
               
    ## Plot in Psi seconda versione
    
    fig05, ax05 = plt.subplots(nrows=1, sharex=True, num='PSI- Time trend - Bands')
    # err_xm va diviso per due o no?, cosa FA L'ERROR BAR DEL PLOT PRECEDENTE?
    err_xm = np.array(err_xm)
    linew = 0.7
    ax05.lw = 0.5
    ax05.plot(time_ts,xm,lw=linew,color='darkolivegreen',label='HRTS w-mean')
    ax05.fill_between(time_ts,xm+err_xm,xm-err_xm, color='g',alpha=0.3) # xm is the weighted mean for HRTS
    ax05.plot(time_ece,temp_eceM,lw = linew, color='orange', label='ECE mean' )             
    ax05.fill_between(time_ece,temp_eceM+err_eceM,temp_eceM-err_eceM,color = 'darkorange',alpha=0.3)
    ax05.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax05.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax05.legend(fontsize=8)
    ax05.set_title(f'JPN {shot} - Ece-Mean and HRTS w-Mean for {psi1}<PSI<{psi2} trends vs time ')
    ax05.set_ylabel('Te (keV)') 
    ax05.set_xlabel('time (s)')
    fig05.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te_MEANS(psi)_BANDS.pdf',dpi=300)
    
    # PLOT cfr dati smootati CON ERRORBARS ece aritmetica vs ts pesata
    fig005, ax005 = plt.subplots(nrows=1, sharex=True, num=' Smooted - PSI- Time trend - Bands')
    # err_xm va diviso per due o no?, cosa FA L'ERROR BAR DEL PLOT PRECEDENTE?
    err_xm = np.array(err_xm)
    linew = 0.7
    ax005.lw = 0.5
    ax005.plot(time_ts,xm_sm,lw=linew,color='darkolivegreen',label='HRTS w-mean')
    ax005.fill_between(time_ts,xm_sm+err_xm_sm,xm_sm-err_xm_sm, color='g',alpha=0.3) # xm is the weighted mean for HRTS
    ax005.plot(time_ece,temp_eceM_sm,lw = linew, color='orange', label='ECE mean' )             
    ax005.fill_between(time_ece,temp_eceM_sm+err_eceM_sm,temp_eceM_sm-err_eceM_sm,color = 'darkorange',alpha=0.3)
    ax005.axvline(x = tlim1,c='r',ls='--',lw=.5)
    ax005.axvline(x = tlim2,c='r',ls='--',lw=.5)
    ax005.legend(fontsize=8)
    ax005.set_title(f'JPN {shot} -Smooted Ece-Mean and Smooted HRTS w-Mean for {psi1}<PSI<{psi2} trends vs time ')
    # ax005.set_title(f'Ece and HRTS Te for {psi1}<PSI<{psi2} vs time ')
    ax005.set_ylabel('Te (keV)') 
    ax005.set_xlabel('time (s)')
    fig005.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Te_MEANS(psi)_BANDS_Smooted.pdf',dpi=300)
         
    return 1

##################################################
# Calcolo delle medie nell'intervallo di psi scelto tra gli estremi psi1 e psi2

def meancalc_12(d, vars): 
    tlim1 = d['tlim1']
    tlim2 = d['tlim2']
    # psi1 = d['psi1']
    # psi2 = d['psi2']
    eP = d['eP']
    tTs = vars['tTs']  
    errTs = vars['errTs']
    tEce = vars['tEce'] 
    psiTs = vars['psiTs']
    psiKk1 = vars['psiKk1']
    psi1 = vars['psi1']
    psi2 = vars['psi2']
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
    # Calcolo delle medie nell'intervallo di psi tra psi1 e psi2 nell'intervallo di tempo scelto
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
    # psi1 = d['psi1']
    # psi2 = d['psi2']
    xm12 = vars['xm12']
    err_xm12 = vars['err_xm12']
    temp_eceM12 = vars['temp_eceM12']
    err_eceM12 = vars['err_eceM12']
    psi1 = vars['psi1']
    psi2 = vars['psi2']
    ####################
    
    dimTs = len(xm12)
    dimEce = len(temp_eceM12)
    
    temp_eceM_R  = signal.resample_poly(temp_eceM12, up=dimTs, down=dimEce) # posizioni ricampionate dei valori TeEce
    err_eceM_R = signal.resample_poly(err_eceM12, up=dimTs, down=dimEce)    # posizioni ricampionate dei valori Err TeEce
 
    left = min(np.nanmin(temp_eceM_R),np.nanmin(xm12))
    right  = max(np.nanmax(temp_eceM_R),np.nanmax(xm12))
    def retta(x):
        return x
    x = np.linspace(left-0.5,right+0.5,dimTs)   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    q = xm12            # Valori TeHRTS media-pesata calcolata sull'intervallo di psi specificato
    p = temp_eceM_R   # Valori TeECE media calcolata sull'intervallo di psi specificato
    # Plot della Te TS vs Te ECE + retta x=y
    titlefig = 'Tts_vs_Tece - PSI'
    
    fig06,ax06 = plt.subplots(nrows=1, sharex=True, num = titlefig)
    ax06.errorbar(p, q, xerr = err_eceM_R, yerr = err_xm12, 
                  marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5,label='ECE vs TS')   
    ax06.plot(x,y,'g--', lw=.8)
    # ax06.set_xlim(left=4)
    # ax06.set_ylim(bottom=4)
    ax06.legend()
    ax06.set_title(f'JPN {shot} - Te HRTS vs Te Ece-Michelson.{psi1}<PSI<{psi2} for {tlim1:.2f}<t<{tlim2:.2f} (s)')
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
    # rho1 = d['rho1']
    # rho2 = d['rho2']
    eP = d['eP']
    tTs = vars['tTs']  
    errTs = vars['errTs']
    tEce = vars['tEce'] 
    rhoTs = vars['rhoTs']
    rhoEce = vars['rhoEce']
    rho1 = vars['rho1']
    rho2 = vars['rho2']
    ####################
    
    idt2 = (tTs.t >= tlim1) & (tTs.t <= tlim2) # Seleziono gli indici dei tempi di interesse
    temp_ts2 = tTs.v[:,idt2]/1000 # in keV
    err_ts2 = errTs.v[:,idt2]/1000  
    err_ts2[err_ts2>2] = 1 # metto a 1 keV l'errore sui punti dove diverge
    rho_ts2 = rhoTs[:,idt2] 
                 
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
    
    idtE2 = (tEce.t >= tlim1) & (tEce.t <= tlim2) # Seleziono gli indici dei tempi di interesse
    temp_ece2 = tEce.v[:,idtE2]/1000  # in keV
    err_ece2 = temp_ece2*eP
    rho_ece2 = rhoEce[:,idtE2]
    
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
    # rho1 = d['rho1']
    # rho2 = d['rho2']
    temp_tsM_rho = vars['temp_tsM_rho']
    err_tsM_rho = vars['err_tsM_rho']
    temp_eceM_rho = vars['temp_eceM_rho']
    err_eceM_rho = vars['err_eceM_rho']
    rho1 = vars['rho1']
    rho2 = vars['rho2']
    ####################
    
    dimTs = len(temp_tsM_rho)
    dimEce = len(temp_eceM_rho)
    
    temp_eceM_rho_R  = signal.resample_poly(temp_eceM_rho, up=dimTs, down=dimEce) # posizioni ricampionate dei valori TeEce
    err_eceM_rho_R = signal.resample_poly(err_eceM_rho, up=dimTs, down=dimEce)    # posizioni ricampionate dei valori Err TeEce
    
    left = min(np.nanmin(temp_eceM_rho), np.nanmin(temp_tsM_rho))
    right  = max(np.nanmax(temp_eceM_rho), np.nanmax(temp_tsM_rho))
    def retta(x):
        return x
    x = np.linspace(left-0.5,right+0.5,dimTs)   # Range in keV di dove tracciare la retta
    y = retta(x)  
    
    q = temp_tsM_rho             # Valori TeHRTS media-pesata calcolata sull'intervallo di psi specificato
    p = temp_eceM_rho_R   # Valori TeECE media calcolata sull'intervallo di psi specificato
    # Plot della Te TS vs Te ECE + retta x=y
    titlefig = 'Tts_vs_Tece - RHO'
    
    fig07,ax07 = plt.subplots(nrows=1, sharex=True, num = titlefig)
    ax07.errorbar(p, q, xerr = err_eceM_rho_R, yerr = err_tsM_rho, 
                  marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5,label='ECE vs TS')   
    ax07.plot(x,y,'g--', lw=.8)
    # ax07.set_xlim(left=4)
    # ax07.set_ylim(bottom=4)
    ax07.legend()
    ax07.set_title(f'JPN {shot} - Te HRTS vs Te Ece-Michelson. \n {rho1}<RHO<{rho2} for {tlim1:.2f}<t<{tlim2:.2f} (s)') #:.2f per avere 2 cifre decimali
    # ax07.set_title(f'Te HRTS vs Te Ece-Michelson. \n {rho1}<RHO<{rho2} for {tlim1}<t<{tlim2} (s)')
    ax07.set_xlabel('Te Ece-Michelson (keV)')
    ax07.set_ylabel('Te HRTS (keV)')
    # fig07.tight_layout()
    
    if d['savefigs'] == 1: 
        plt.savefig(d['mypath']+f'{shot}_Tts_vs_Tece_RHO.pdf',dpi=300)
            
    return 1
    
 ##############################
 ######################## Multiplot
 
def multiplot(d,vars): 
     
    shot = d['shot']      # 103117   # 99950, 99971
    tlim1 = d['tlim1']    #43     # limite inferiore selezione tempi - in secondi
    tlim2 = d['tlim2']    #51.4 
    rad = d['rad']        #3.00 
       
    w = ppfs(shot) 
    
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
        
