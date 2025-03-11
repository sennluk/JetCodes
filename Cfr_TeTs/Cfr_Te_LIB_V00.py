"""
Created on 28/01/2025 
@author: lsenni. 
please refer to word file: Note_Codice_Cfr_Te_ts
Based on the myelab_V12.py library file
no more PSI averages
Functions:
    dictvar : crea il dizionario della variabili
    tdef: calcola gli estremi ti e tf per la finestra temporale
    multiplot: fornische una multiplot con le caratteristiche principali dello shot
    magax: fornisce l'andamento delle posizione dell'ase magnetico secondo EFIT sul piano poloidale
    psicalc: calcola le coordinate PSI per le linee di vista delle due diagnostiche
    rhocalc: calcola le coordinate RHO per le linee di vista delle due diagnostiche
    tprof: mostra l'andamento del profilo di Te per le due ddiagnostiche
    rhofig
N.B.: Lascio i valori di OPSI come definiti in d--> per controlli e calcoli possibili    
    '

"""
import numpy as np
import my_flush
from ppfeg import ppfs
# from ppfeg import ppfs
# import matplotlib.pyplot as plt
from scipy import interpolate
# from scipy.signal import savgol_filter
import my_manage_file as mym
import matplotlib.pyplot as plt
# import time

# DTE3 shot list: 104990,991,994,995,999   
# 104520,521,522,523,524,526
##################################################
# Deinisco il dizionario delle variabili

def dictvar(d):
    shot = d['shot']
    # timestr = time.strftime("%Y%m%d-%H%M%S")
    
    w = ppfs(shot)   
    
    if d['savefigs'] == 1: 
        mym.create_folder(d)
        mym.save_param(d)
        # mym.save_print_output_to_file(d,vars)
    else:
        print('You are not saving the plots')
     
    vars = {
        'w' : w,                           # canale dati del ppfeg con tutti i dati relativi allo sparo in oggetto
        'ti' : None, # Tempo iniziale finestra analisi
        'tf' : None, # Tempo finale foinestra analisi
        'tlim' : None,          # tempo medio al quale vengono calcolati i profili per plot di controllo
        'Ti' : None,  # Temp HRTS/max dell'istante iniziale
        'Tf' : None,  # Temp HRTS/max dell'istante finale
        'psi1' : d['psi1'],    # inizializzo i valori di psi1,2 e rho1,2
        'psi2' : d['psi2'],    # questi possono essere usati se non si vuole usare
        'rho1' : None,    # quelli calcolati con porocedura automatica
        'rho2' : None,
        'rhodown' : None,     # rho
        'rhoup' : None,
        'tTs' : w.hrts.te,                 # Chann Te HRTS  
        'errTs' : w.hrts.dte,              # Chann Errors HRTS.dte
        'tEce' : w.ecm1.prfl,              # Chann Te Ece-Kk1
        'errEce' : (w.ecm1.prfl)*d['eP'],  # Error associated to the Kk1 meas: ab 2%
        'psiTs' : w.hrts.psi,              # Chann PSI HRTS 
        'psiKk1' : None,
        'psiTscalc' : None,
        'rhoKk1' : None,
        'time_ts' : None,
        'time_ece' : None,
        'rhoTs' : None,
        'rhoEce' : None, 
        'temp_tsM' : None, 
        'err_tsM': None,
        'timeTs2' : None,  # TS time instants in the ti-tf window
        'timeEce22' : None, # ECE time instants closest to the TS ones in the ti-tf window     '
        'ranges' : None,  # RHO ranges as calculated by the automatic procedure - TS interpolated in t1<t<t2
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
# based on the temperature raising

def tdef(d,vars):
    # shot = d['shot']
    w = vars['w']
    window_size = d['window_size']
    
    #####################################
    tmax = w.hrtx.tmax/1000     # max Te max hrts (hrtx channel) in keV
    data = tmax.v[0,:]

    def moving_average_np(data, window_size):
        """Calcola la media mobile semplice usando NumPy."""
        return np.convolve(data, np.ones(window_size) / window_size, mode='valid')

    def find_rising_point_np(data, window_size=3):
        """Trova il primo punto in cui i dati iniziano a salire."""
        # Calcola la media mobile
        smoothed_data = moving_average_np(data, window_size)
        
        # Trova il primo punto in cui il valore aumenta
        for i in range(1, len(smoothed_data)):
            if smoothed_data[i] > smoothed_data[i - 1]:
                return i + window_size - 1  # Indice relativo all'array originale
        
        return None  # Nessun punto di crescita trovato
        
    def find_return_point_np(data, start_index, target_value):
        """Trova il punto in cui i dati tornano al valore target dopo l'indice dato."""
        for i in range(start_index + 1, len(data)):
            if data[i] <= target_value:  # Modifica qui se vuoi confrontare con <= o ==
                return i
        return None  # Nessun punto di ritorno trovato

    rising_index = find_rising_point_np(data, window_size)

    # if rising_index is not None:
    #     print(f"I dati iniziano a salire all'indice {rising_index}, valore: {data[rising_index]}")
    # else:
    #     print("Non ci sono punti in cui i dati iniziano a salire.")
    if rising_index is not None:
        rising_value = data[rising_index]
        print(f"I dati iniziano a salire all'indice {rising_index}, valore: {rising_value}")
        
        # Trova il punto di ritorno
        return_index = find_return_point_np(data, rising_index, rising_value-0.2)
        if return_index is not None:
            print(f"I dati tornano al valore {rising_value} all'indice {return_index}")
        else:
            print(f"I dati non tornano più al valore {rising_value}.")
    else:
        print("Non ci sono punti in cui i dati iniziano a salire.")

    Ti = rising_value          # Temp istante inizale
    Tf = data[return_index]    # Temp istante finale
       
    ti = tmax.t[rising_index]  # Istante iniziale
    tf= tmax.t[return_index]   # Istante finale

    tlim1 = vars['ti']
    tlim2 = vars['tf']
    # prendo come tlim l'istante in cui la Tmax è massima
    tlim = tmax.t[np.nanargmax(tmax.v)] 
    
    vars['ti'] = ti
    vars['tf'] = tf
    vars['Ti'] = Ti
    vars['Tf'] = Tf
    vars['tlim'] = tlim
    
    print('JPN = ', d['shot'])
    print("Automated tlim1 =", ti)
    print("Automated tlim2 =", tf)
    print('t lim1 = ', tlim1)
    print('t lim2 = ', tlim2)
    print('Instant of the max Tmax=', tlim)
    
    return 1
##################################################
# Compute PSI on ECE and HRTS LoSs
def psicalc(d, vars):
    shot = d['shot']
    w = vars['w']
  
    tEce = w.ecm1.prfl 
    zKk1 = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista (in metri), 
    # la seconda è l'intersezione con l'asse y. 'Position of antenna aperture centre (z,majR) '
    rEce = tEce.r
    z = np.full_like(rEce, zKk1) # Retta linea di vista parallela asse r
    
    timeEce = tEce.t

    psiKk1 = np.zeros(tEce.v.shape)
    
    for i,tempo in enumerate(timeEce):
        ts, ier = my_flush.flushinit(15, shot, tempo) # Definisco il tipo di equilibrio da usare
        # e Prendo il tempo più vicino a quello desiderato
        psi, _ = my_flush.Flush_getFlux(rEce*100, z*100) # Le coordinate vanne messe in cm!
        psiKk1[:,i] = psi
   
    vars['psiKk1'] = psiKk1      
    
    return 1
##################################################
# Compute RHO on ECE and HRTS LoSs
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
         
    for i,tempo in enumerate(tTs.t):
        psiTs_ = np.abs(psiTs.v[:,i])
        index = np.argmin(abs(torFlux.t - tempo))
        vFlux_ = vFlux[:,index]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1]) # function to interp the flux
        fl_int_hrts = fl_int(psiTs_, extrapolate=False) # Interpolation without extrapolating
        rhoTs[:,i] = np.sqrt(fl_int_hrts)
        
    for i,tempo in enumerate(tEce.t):
        psi = psiKk1[:,i]
        index = np.argmin(abs(torFlux.t - tempo))
        vFlux_ = vFlux[:,index]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1])  # function to intertp the flux
        fl_int_ece = fl_int(psi, extrapolate=False)
        rhoEce[:,i] = np.sqrt(fl_int_ece)
    
    vars['rhoTs'] = rhoTs
    vars['rhoEce'] = rhoEce
        
    return  1
    
################################################## Ranges and averages
# Computes the RHO ranges for each of the slowest Diagnostics (HRTS)
# compute the values of ECE nearest to the same instants
# Perform the averages on the selctes RHO ranges, instant by instnat 
# of the Te values and errors for both ECE and HRTS
# il 2 finale indica la grandezza nell'intervallo ti-tf
# il 22 indica che la medesima grandezza è stata calcolate nei punti piuù vicini
# a quella più lenta, sempre nellintervallo ti-tf
 
def def_range_av(d,vars):
   #######################################################
   # Function to compute the RHO limits
   #def rhorange(d, vars): 
   shot = d['shot']
   eP = d['eP']
   ti = vars['ti']
   tf = vars['tf'] 
   w = vars['w']
   tTs = vars['tTs']  
   errTs = vars['errTs']
   tEce = vars['tEce'] 
   # psiTs = vars['psiTs']
   # psiEce = vars['psiKk1']
   rhoTs = vars['rhoTs']
   rhoEce = vars['rhoEce']

   tlim1 = ti
   tlim2 = tf

   idt1 = (tTs.t >= tlim1) & (tTs.t <= tlim2) # Seleziono gli indici dei tempi di interesse
   rhoTs2 = rhoTs[:,idt1] 
   timeTs2 = tTs.t[idt1]
   temp_ts2 = tTs.v[:,idt1]/1000 # in keV
   err_ts2 = errTs.v[:,idt1]/1000  
   err_ts2[err_ts2>2] = 1 # metto a 1 keV l'errore sui punti dove diverge

   idt2 = (tEce.t>= tlim1) & (tEce.t <= tlim2)
   rhoEce2 = rhoEce[:,idt2]
   timeEce2 = tEce.t[idt2]
   temp_ece2 = tEce.v[:,idt2]/1000  # in keV
   err_ece2 = temp_ece2*eP

   print('Dimensioni nell intervallo temporale di rhoTs = ',rhoTs2.shape)
   print('Dimensioni nell intervallo temporale di rhoEce = ',rhoEce2.shape)

   dimTs = rhoTs2.shape[1]    # number of Time points TS
   dimEce = rhoEce2.shape[1]  # number of Time points TS
   rTs = tTs.r
   rEce = tEce.r

   #######################################
   # Inserisco ciclo per il calcolo dell'intervallo di RHo in tutti gli istanti di TS tra t1-t2
   # Rho range estimation - Consider only the ECE time scale

   d2 = 0.08 # mezzo intervallo intorno al centro- in centimetri
   ranges = np.zeros([timeTs2.shape[0],2])  # inizializzo la matrice dei range in rho, uno per istante

   temp_tsM_rho = np.zeros_like(timeTs2)      #inizializzo matrice Tts valori temp in t1-t2
   err_tsM_rho = np.zeros_like(timeTs2)   #inizializzo matrice errTs valori Err-temp in t1-t2 

   timeEce22 = np.zeros_like(timeTs2)

   # valori tra t1 e t2 più vicini a quelli di TS
   rhoEce22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])
   temp_ece22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])
   err_ece22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])

   # inizializzo le mettrici delle medie, una per ogni tempo tra t1 e t2 basandosi sul campionamenti di time-TS
   temp_eceM_rho = np.zeros_like(timeTs2) 
   rho_eceM_rho = np.zeros_like(timeTs2)
   err_eceM_rho = np.zeros_like(timeTs2) 


   for i,inst  in enumerate(timeTs2):
       timeTs = inst                         # valore tempo per TS
       timeEce = np.min(abs(timeEce2-inst))  # valore tempo per Ece pià vicino
       iEce = np.argmin(abs(timeEce2-inst))  # valore indice del tempo Ece
       timeEce22[i] = timeEce2[iEce] 
       rhoEce22[:,i] = rhoEce2[:,iEce] 
       temp_ece22[:,i] = temp_ece2[:,iEce] 
       err_ece22[:,i] = err_ece2[:,iEce] 
       
       mTs= np.nanmin(rhoTs2[:,i])     # min rho Ts at the ith indice    
       indT = np.nanargmin(rhoTs2[:,i])   
       posT = rTs[indT]  
       
       mEce = np.nanmin(rhoEce2[:,iEce])   # min rho Ece at the i-th indice
       indE = np.nanargmin(rhoEce2[:,iEce]) # indice del minimo sopra calcolato
       posE = rEce[indE]                  # Position of the rho ECE minimum val
       
       M = np.max([mEce,mTs])        # Trovo il massimo tra i due minimi
       ind = np.argmax(([mEce,mTs])) # per definire la diagnostica sulla Los della quale si trova il max:
                                     # ind = 0 il masssimo dei due minimi appartiene all'ECE...
       diag = [rhoEce2, rhoTs2]      #vettore con le due posizioni dei minimi 
       rDiag = [rEce,rTs]            # Costriusco un vettore con i due vettori delle posizioni
       pos = [posE,posT]             
       # RHo - UP : calcolo d2 cm in meno e in più
       # rispetto alla poszione del massimo dei due minimi
       # e vedo il punto più vicino della medesima curva
       # Come RHO-UP prendo poi il massimo dei due valori trovati + 0.001
       
       diagM = diag[ind]  # La diagnostica che ha il massimo tra i due minimi
       radii = rDiag[ind]     # Le posiszioni della diagnostica di cui sopra
       
       pos1 = pos[ind]- d2  # Posizione a distanza d/2 dal centro plasma
       pos2 = pos[ind] + d2
       
       idt1 = np.nanargmin(abs(radii-pos1)) # indice1 più vicino sulla curva più alta
       idt2 = np.nanargmin(abs(radii-pos2))
       # vettore dei valori più vicino alla posizione estrema dell'intervallo 
       #(per la diagn che ha il max tra i due min)
       temp = [iEce,i]
       indice1 = temp[ind]
       pippo = [v1,v2] = [diagM[idt1,indice1],diagM[idt2,indice1]] 
        
       siaM = np.max(pippo) # Seleziono il pun to più alto tra i due, da scegliere come rho-up
       rhoup = siaM + 0.0001 # aggiungo un piccolo delta di sicurezza
       
       #############################
       # RHO down
       # Trovo il punto più vicino della diag non considerata prima
       # seleziono l'altra diagnostica: minimo rta i due massimi
       indx = np.argmin(([mEce,mTs]))   # Trovo il minimo tra i due minimi
       rm = rDiag[indx]     # Le posiszioni della diagnostica di cui sopra
       diagm = diag[indx]      # seleziono la diagnostica 'piu bassa'
       # trovo l'indice del punto con rho piu vicina al minimo dell'altra:
       # avendo assi temporali diversi devo secegliere indici diversi a secondo del caso    
       temp = [iEce,i]
       indice2 = temp[indx]
       id_pluto = np.nanargmin(abs(diagm[:,indice2]-M)) 
       pluto = diagm[id_pluto,i]
       
       rm = rDiag[indx]     # posizioni radiali della diagostica 'minore'
       pos_pluto = rm[id_pluto]
       
       # se il punto è prima del minimo aggiungo due punto per determinare la rho down
       # se invece è dopo, ne tolgo due
       # in ambedue i casi tolgo poi un piccolo margine
       if id_pluto < indT:
           rhodown = diagm[id_pluto+2,indice2] - 0.001
       else:
           rhodown = diagm[id_pluto-2,indice2] - 0.001
       ranges[i,:] = [rhodown,rhoup]       
       
       mask = (rhoTs2[:,i] >= rhodown) & (rhoTs2[:,i] <= rhoup)    # Seleziono gli indici dei tempi di interesse 
       temp = temp_ts2[mask,i]                    # Valori di temperaatura nelle rho selezionate
       sig = err_ts2[mask,i]       # Errore sui singoli valori di temperatura: sigma
       ai = 1/sig**2            # Inverso del quadrato delle sigma
       somma = sum(ai)
       if somma==0:
           somma=1
       xm_ = sum(ai*temp)/somma   # Media dei valori di temperatura pesata sui singoli errori
       temp_tsM_rho[i] = (xm_)               # Valore delle temperature 'attese'/più probabili, nel tempo
       err_xm_ = np.sqrt(1/somma) 
       err_tsM_rho[i] = (err_xm_)

       mask1 = (rhoEce22[:,i] >= rhodown) & (rhoEce22[:,i] <= rhoup)    # Seleziono gli indici delle psi di interesse 
       temp_eceM_rho[i]  = np.mean(temp_ece22[mask1,i], axis=0)
       rho_eceM_rho[i] = np.mean(rhoEce22[mask1,i],axis=0)
       err_eceM_rho[i] = np.mean(err_ece22[mask1,i], axis=0)

      
   print('Dimensioni temp_eceM_rho = ',temp_eceM_rho.shape)
   print('Dimensioni temp_TsM_rho = ',temp_tsM_rho.shape)    
   
   vars['temp_tsM_rho'] = temp_tsM_rho
   vars['err_tsM_rho'] = err_tsM_rho
   vars['temp_eceM_rho'] = temp_eceM_rho
   vars['err_eceM_rho'] = err_eceM_rho
   vars['timeTs2'] = timeTs2
   vars['timeEce22'] = timeEce22
   vars['ranges'] = ranges
   # Si trova la matrice 'ranges' con tutti i valoori di rhoup/down a tutti gli istanti del TS
   # Si trovano qiundi le due matrici: temp_eceM_rho e temp_tsM_rho 
   # che contengono i valori di Te a tutte le rho e a tutte gli istanti tra t1 e t2
   # assi temporali: timeTs2 e timeEce2
   # le posizioni: rTs e rEce
    
   # fig, ax = plt.subplots(3)
   # ax[0].plot(timeTs2,timeEce22, label='cfr tempi')
   # ax[0].legend()
   # ax[0].set_title('Check Plot')
   # ax[1].plot(timeTs2,ranges[:,1], label='rho-up')
   # ax[1].plot(timeTs2,ranges[:,0], label='rho-down')
   # ax[1].legend()
   # ax[2].plot(timeTs2, temp_tsM_rho, label='TS averaged')
   # ax[2].plot(timeEce22,temp_eceM_rho, label='ECE averaged')
   # ax[2].legend()

   # left = min(np.nanmin(temp_eceM_rho), np.nanmin(temp_tsM_rho))
   # right  = max(np.nanmax(temp_eceM_rho), np.nanmax(temp_tsM_rho)) 
   # def retta(x):
   #     return x

   # x = np.linspace(left-0.5,right+0.5,timeTs2.shape[0])   # Range in keV di dove tracciare la retta
   # y = retta(x)  

   # fig,ax=plt.subplots(1)
   # # ax.scatter(temp_eceM_rho, temp_tsM_rho,label='Rho Averaged ECE vs TS ') 
   # ax.plot(x,y,'g--', lw=.8)
   # ax.errorbar(temp_tsM_rho, temp_eceM_rho, xerr = err_tsM_rho, yerr = err_eceM_rho, 
   #               marker='o', markersize=3, ecolor='g', linestyle='none', elinewidth=.5, label='ECE vs TS')
   # ax.set_title(f'JPN {shot} - Te Ece-Michelson vs Te HRTS  for {tlim1:.2f}<t<{tlim2:.2f} (s)') #:.2f per avere 2 cifre decimali
   # ax.set_xlabel('Te HRTS (keV)')
   # ax.set_ylabel('Te Ece-Michelson (keV)')
   # ax.legend()
    
    
    