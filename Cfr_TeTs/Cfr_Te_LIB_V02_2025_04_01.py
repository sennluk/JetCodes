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
V02: cambio l'importazione dei dati per non avere nessun ppfs nel dizionario
e poterlo salvare in json'
"""
import numpy as np
import my_flush
from ppfeg import ppfs
from scipy import interpolate
import my_manage_file as mym
##################################################
# Definisco il dizionario delle variabili

def dictvar(d):
    print('------------ Start new shot!!! ------------')
 
    vars = {
        'tTs_v' : None,     # HRTS Te for all the times and positions 
        'tTs_t' : None,     # HRTS time instant values
        'tTs_r' : None,     # HRTS positions of the measurements in Major Radius coord 
        'errTs' : None,     # HRTS errors for each time and pos
        'psiTs_v' : None,   # HRTS PSI values for the LoS at each time
        'psiTs_t' : None,   # time for HRTS PSI values for the LoS at each time
        'psiTs_r' : None,   # position for HRTS PSI values for the LoS at each time
        'zTs_v' : None,     # HRTS z Coords of the LoS
        'zTs_r' : None,     # HRTS R Coords of the LoS
        'tmax_v' : None,      # HRTS max Te for each time (Tmax channel)
        'tmax_t' : None,      # HRTS max Te for each time (Tmax channel) - time
        'tEce_v' : None,    # ECE-KK1 Te for all the times and positions 
        'tEce_t' : None,    # ECE-KK1time instant values
        'tEce_r' : None,    # ECE-KK1 positions of the measurements in Major Radius coord 
        'errEce' : None,    # ECE-KK1 estimated errors for each time and pos
        'flux_v': None,     # Toroidal Flux coordinates from EFIT 
        'flux_t' : None,    # Instants of the  Flux coordinates
        'flux_r' : None,    # Major radius Positions of the  Flux coordinate
        'rmag' : None,      # R-Posisitions of the magnetic axis ovet time (EFIT)
        'zmag' : None,      # z-Posisitions of the magnetic axis ovet time (EFIT)
        'tmag' : None,      # Time of the magnetic axixs ovet time (EFIT) computing
        'ti' : None,        # Initial time instant considered for the analisys
        'tf' : None,        # Final time instant considered for the analisys
        'tlim' : None,      # Time at which Tmax HRTS is max
        'maxT_hrts' : None, #  Max Te of the Tmax (HRTS) channel over the whole shot
        'Ti' : None,        # Temp HRTS/max dell'istante iniziale
        'Tf' : None,        # Temp HRTS/max dell'istante finale
        'rad' : None,       # Position of the plasma center a
        'psi1' : d['psi1'], # PSI values manually chosen 
        'psi2' : d['psi2'], 
        'rhodown' : None,   # RHO lower value considered for the RHO range - Average
        'rhoup' : None,     # RHO upper value considered for the RHO range - Average
        'zKk1' : None,      # z Coord of the ECE-KK1 LoS 
        'psiKk1' : None,    # PSI coordinates on the ECE-KK1 data points
        'psiTscalc' : None, # HRTS PSI Computed to check the correctnes of the data Channel iused
        'rhoKk1' : None,    # ECE-KK1 LoS PSI coordinates 
        'rhoTs' : None,     #RHO coordinates on the HRTS data points
        'rhoEce' : None,    # RHO coordinates on the ECE-KK1 data points
        'timeTs2' : None,     # TS time instants in the ti-tf window
        'timeEce22' : None,   # ECE time instants closest to the TS ones in the ti-tf window     '
        'ranges' : None,      # RHO ranges as calculated by the automatic procedure at each HRTS time- TS interpolated in t1<t<t2
        'temp_tsM_rho' : None,   # Te values averaged over the rho1-rho2 interval for HRTS
        'err_tsM_rho' : None,    # Error on HRTS Te values averaged over the rho1-rho2 interval
        'temp_eceM_rho' : None,  # Te values averaged over the rho1-rho2 interval for ECE
        'err_eceM_rho' : None,   # Error on ECE Te values averaged over the rho1-rho2 interval
        'ratio' : None,          # Ratio between temperaure values (averaged) - TO BE double-checked!
        'distance' : None,       # Difference between temperaure values (averaged)
        'err_dist' : None,       # Computed error to be associated to the difference
        'err_dist_perc' : None   # Computed error to be associated to the difference/Tece
        } 
    
    shot = d['shot']
    # timestr = time.strftime("%Y%m%d-%H%M%S")
    # Import data of HRTS and ECE KK1 and put it in the Dictionary
    w = ppfs(shot)        # PPFS channel containing all data for the shot
    
    vars['tTs_v'] = w.hrts.te.v     # HRTS Te for all the times and positions 
    vars['tTs_t'] = w.hrts.te.t     # HRTS time instant values
    vars['tTs_r'] = w.hrts.te.r     # HRTS positions of the measurements in Major Radius coord 
    vars['errTs'] = w.hrts.dte.v    # HRTS errors for each time and pos
    vars['psiTs_v'] = w.hrts.psi.v      # HRTS PSI values for the LoS at each time
    vars['psiTs_t'] = w.hrts.psi.t      # HRTS PSI values for the LoS at each time
    vars['psiTs_r'] = w.hrts.psi.r      # HRTS PSI values for the LoS at each time
    vars['zTs_v'] = w.hrts.z.v
    vars['zTs_r'] = w.hrts.z.r
    vars['tmax_v'] = w.hrtx.tmax.v      # HRTS max temperature for each time - values
    vars['tmax_t'] = w.hrtx.tmax.t      # HRTS max temperature for each time - time instants
      
    vars['tEce_v'] = w.ecm1.prfl.v  # ECE-KK1 Te for all the times and positions 
    vars['tEce_t'] = w.ecm1.prfl.t  # ECE-KK1time instant values
    vars['tEce_r'] = w.ecm1.prfl.r  # ECE-KK1 positions of the measurements in Major Radius coord 
    vars['errEce'] = vars['tEce_v']*d['eP'] # ECE-KK1 estimated errors for each time and pos
    vars['zKk1'] = w.ecm1.antp.v[1,0]   # antp è il canale con le coordinate della linea di vista (in metri), 
    # la seconda è l'intersezione con l'asse y. 'Position of antenna aperture centre (z,majR) '
    
    vars['flux_v'] = w.efit.ftor.v  # Toroidal Flux values (EFIT) for each time and position
    vars['flux_t'] = w.efit.ftor.t  # Tor Flux times
    vars['flux_r'] = w.efit.ftor.r  # Tor Flux positions
    vars['rmag'] = w.efit.rmag.v    # r Magnetic Axe position
    vars['tmag'] = w.efit.rmag.t
    vars['zmag'] = w.efit.zmag.v    #  z Magn Ax pos
    
    if d['savefigs'] == 1: 
        mym.create_folder(d)
        mym.save_param(d)
        # mym.save_print_output_to_file(d,vars)
    else:
        print('You are not saving the plots')
        
    print(f'---- Dict vars JPN {shot} created! ----')
    return vars

##################################################
# Define the time interval of the analysis for the single shot
# based on the temperature raising

def tdef(d,vars):
    Tref = d['Tref']
    window_size = d['window_size']
    min_increase = d['min_increase']
    tmax_v = vars['tmax_v']
    tmax_t = vars ['tmax_t']
    #####################################
    # tmax = vars['tmax_v']/1000     # max Te max hrts (hrtx channel) in keV
    data = tmax_v[0,:]/1000
    # prendo come tlim l'istante in cui la Tmax è massima
    # valid_indices = np.where((tmax_t > 42) & (tmax_t <= 60))[0]
    valid_indices = np.where(tmax_t <= 60)[0]
    tlim = tmax_t[np.nanargmax(tmax_v[0,valid_indices])] 
    ind_tlim = np.nanargmax(tmax_v[0,valid_indices])
    TMAX = data[ind_tlim]
    
    def moving_average_np(data, window_size):
        """Calcola la media mobile semplice usando NumPy."""
        return np.convolve(data, np.ones(window_size) / window_size, mode='valid')

    def find_rising_point_np(data, window_size, threshold, min_increase):
        """Trova il primo punto in cui i dati iniziano a salire."""
        # Calcola la media mobile
        smoothed_data = moving_average_np(data, window_size)
        
        # Trova il primo punto in cui il valore aumenta
        for i in range(1, len(smoothed_data)):
            if smoothed_data[i] > smoothed_data[i - 1]+ min_increase and smoothed_data[i] > threshold:
                return i + window_size - 1  # Indice relativo all'array originale
        
        return None  # Nessun punto di crescita trovato
        
    def find_return_point_np(data, start_index, target_value):
        """Trova il punto in cui i dati tornano al valore target dopo l'indice dato."""
        for i in range(start_index + 1, len(data)):
            if data[i] <= target_value and i > i: #+20:  # Modifica qui se vuoi confrontare con <= o ==
                return i
            data_subset = data[ind_tlim + 1:]  # Prendi solo la parte dopo start_index
            closest_index = np.argmin(np.abs(data_subset - target_value))  # Trova il valore più vicino
            closest_index += start_index + 1  # Riporta l'indice nel riferimento originale
       
        return closest_index  # Restituisce l'indice del valore più vicino
    #####################
    rising_index = find_rising_point_np(data, window_size, Tref, min_increase)
    
    if rising_index is not None:
        rising_value = data[rising_index]
        print(f"I dati iniziano a salire all'indice {rising_index}, valore: {rising_value}")
        
        # Trova il punto di ritorno
        return_index = find_return_point_np(data, ind_tlim, rising_value-0.2)
        Tfin = data[return_index]
        if return_index is not None:
            print(f"L'indice finale è: {return_index}, corrispondente ad una valore: {Tfin}")
        else:
            print(f"I dati non tornano più al valore {rising_value}.")

    #####################
    Ti = rising_value          # Temp istante inizale
    Tf = data[return_index]    # Temp istante finale
       
    ti = tmax_t[rising_index]  # Istante iniziale
    tf= tmax_t[return_index]   # Istante finale

    tlim1 = vars['ti']
    tlim2 = vars['tf']
    
    vars['ti'] = ti
    vars['tf'] = tf
    vars['Ti'] = Ti
    vars['Tf'] = Tf
    vars['tlim'] = tlim
    vars['maxT_hrts'] = TMAX
    
    print('JPN = ', d['shot'])
    print("Automated tlim1 =", ti)
    print("Automated tlim2 =", tf)
    print('t lim1 = ', tlim1, ' None = Automatically computed')
    print('t lim2 = ', tlim2, ' None = Automatically computed')
    print('Instant of the max Tmax=', tlim)
    
    return 1
##################################################
# Compute PSI on ECE LoS coords (HRTS has his channel available) 
def psicalc(d, vars):
    shot = d['shot']
    
    timeEce = vars['tEce_t']
    rEce = vars['tEce_r']
    zKk1 = vars['zKk1']  # antp è il canale con le coordinate della linea di vista (in metri), 
    # la seconda è l'intersezione con l'asse y. 'Position of antenna aperture centre (z,majR) '

    z = np.full_like(rEce, zKk1) # Retta linea di vista parallela asse r
    psiKk1 = np.zeros(vars['tEce_v'].shape)
    
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
    tTs_v = vars['tTs_v']    # Chan. Te HRTS  
    tTs_t = vars['tTs_t']
    tTs_r = vars['tTs_r']
    psiTs_v = vars['psiTs_v']    # Channel psi hrts
    # psiTs_t = vars['psiTs_t']    # Channel psi hrts
    # psiTs_r = vars['psiTs_r']    # Channel psi hrts
    tEce_t = vars['tEce_t']
    psiKk1 = vars['psiKk1']
    ####################
    
    vFlux = vars['flux_v']     # valori di flusso
    rFlux = vars['flux_r']     # Posizioni dei valori di flusso
    tFlux = vars['flux_t']
    
    rhoTs = np.zeros(tTs_v.shape)
    rhoEce = np.zeros(psiKk1.shape)
         
    for i,tempo in enumerate(tTs_t):
        psiTs_ = np.abs(psiTs_v[:,i])
        index = np.argmin(abs(tFlux - tempo))
        vFlux_ = vFlux[:,index]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1]) # function to interp the flux
        fl_int_hrts = fl_int(psiTs_, extrapolate=False) # Interpolation without extrapolating
        rhoTs[:,i] = np.sqrt(fl_int_hrts)
        
    for i,tempo in enumerate(tEce_t):
        psi = psiKk1[:,i]
        index = np.argmin(abs(tFlux - tempo))
        vFlux_ = vFlux[:,index]
        fl_int = interpolate.make_interp_spline(rFlux, vFlux_/vFlux_[-1])  # function to intertp the flux
        fl_int_ece = fl_int(psi, extrapolate=False)
        rhoEce[:,i] = np.sqrt(fl_int_ece)
    
    tlim = vars['tlim']
    idt = np.argmin(abs(tTs_t[:]-tlim))  # index of the maximum Tmax
    rho_ts = rhoTs[:,idt]               # RHO coords at the time of Max Tmax
    minimo = np.min(np.nanmin(rho_ts))  # Minimo escudendo eventuali NaN
    indices = np.where(rho_ts == minimo)  # Come centro plasma prendo la posizione del minimo di Rho TS
    rad_ts = tTs_r[indices[0]]
    rad = rad_ts
    
    vars['rhoTs'] = rhoTs
    vars['rhoEce'] = rhoEce
    vars['rad'] = rad
    
    return  1
    
################################################## Ranges and averages
# Computes the RHO ranges for each of the slowest Diagnostics (HRTS)
# compute the values of ECE nearest to the same instants
# Perform the averages on the seleted RHO ranges, instant by instnat 
# of the Te values and errors for both ECE and HRTS
# il 2 finale indica la grandezza nell'intervallo ti-tf
# il 22 indica che la medesima grandezza è stata calcolate nei punti più vicini
# a quella più lenta, sempre nell'intervallo ti-tf
 
def def_range_av(d,vars):
   eP = d['eP']
   d2 = d['d2']
   numero_punti = d['np']
   ti = vars['ti']
   tf = vars['tf'] 
   tTs_v = vars['tTs_v']
   tTs_t = vars['tTs_t']
   tTs_r = vars['tTs_r']
   errTs = vars['errTs']
   tEce_v = vars['tEce_v']
   tEce_t = vars['tEce_t']
   tEce_r = vars['tEce_r']
   rhoTs = vars['rhoTs']
   rhoEce = vars['rhoEce']
   ###############
   tlim1 = ti
   tlim2 = tf

   idt1 = (tTs_t >= tlim1) & (tTs_t <= tlim2) # Seleziono gli indici dei tempi di interesse
   rhoTs2 = rhoTs[:,idt1]    # 2 stands for : 'in the time range ti-tf
   timeTs2 = tTs_t[idt1]      
   temp_ts2 = tTs_v[:,idt1]/1000 # in keV
   err_ts2 = errTs[:,idt1]/1000  
   err_ts2[err_ts2>2] = 1 # metto a 1 keV l'errore sui punti dove diverge

   idt2 = (tEce_t>= tlim1) & (tEce_t <= tlim2)
   rhoEce2 = rhoEce[:,idt2]
   timeEce2 = tEce_t[idt2]
   temp_ece2 = tEce_v[:,idt2]/1000  # in keV
   err_ece2 = temp_ece2*eP

   print('Dimensioni nell intervallo temporale di rhoTs = ',rhoTs2.shape)
   print('Dimensioni nell intervallo temporale di rhoEce = ',rhoEce2.shape)

   dimTs = rhoTs2.shape[1]    # number of Time points TS
   dimEce = rhoEce2.shape[1]  # number of Time points TS
   rTs = tTs_r
   rEce = tEce_r

   #######################################
   # Inserisco ciclo per il calcolo dell'intervallo di RHo in tutti gli istanti di TS tra t1-t2
   # Rho range estimation - Consider only the ECE time scale

   ranges = np.zeros([timeTs2.shape[0],2]) # inizializzo la matrice dei range in rho, uno per istante

   temp_tsM_rho = np.zeros_like(timeTs2)   #inizializzo matrice Tts valori temp in t1-t2
   err_tsM_rho = np.zeros_like(timeTs2)    #inizializzo matrice errTs valori Err-temp in t1-t2 

   timeEce22 = np.zeros_like(timeTs2)     # 22 stands for:' in the ti-tf range based on the time of HRTS 
                                          # same number of pints- the closest one will be chosen
   # valori tra t1 e t2 più vicini a quelli di TS
   rhoEce22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])
   temp_ece22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])
   err_ece22 = np.zeros([rEce.shape[0],timeTs2.shape[0]])

   # inizializzo le mettrici delle medie, una per ogni tempo tra t1 e t2 basandosi sul campionamenti di time-TS
   temp_eceM_rho = np.zeros_like(timeTs2) 
   rho_eceM_rho = np.zeros_like(timeTs2)
   err_eceM_rho = np.zeros_like(timeTs2) 


   for i,inst  in enumerate(timeTs2):                      # valore tempo per TS
       timeEce = np.min(abs(timeEce2-inst))  # valore tempo per Ece più vicino
       iEce = np.argmin(abs(timeEce2-inst))  # valore indice del tempo Ece
       timeEce22[i] = timeEce2[iEce] 
       rhoEce22[:,i] = rhoEce2[:,iEce] 
       temp_ece22[:,i] = temp_ece2[:,iEce] 
       err_ece22[:,i] = err_ece2[:,iEce] 
       
       mTs = np.nanmin(rhoTs2[:,i])         # min rho Ts at the ith indice    
       indT = np.nanargmin(rhoTs2[:,i])   
       posT = rTs[indT]  
       
       mEce = np.nanmin(rhoEce2[:,iEce])     # min rho Ece at the i-th indice
       indE = np.nanargmin(rhoEce2[:,iEce])  # indice del minimo sopra calcolato
       posE = rEce[indE]                     # Position of the rho ECE minimum val
       
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
       
       diagM = diag[ind]      #Le coord RHo della diagnostica che ha il massimo tra i due minimi
       radii = rDiag[ind]     # Le posiszioni della diagnostica di cui sopra
       
       pos1 = pos[ind] - d2    # Posizione a distanza d/2 dal centro plasma (minimo rho)
       pos2 = pos[ind] + d2
       
       idt1 = np.nanargmin(abs(radii-pos1))    # indice 1 più vicino alla pos 1 sulla curva più alta
       idt2 = np.nanargmin(abs(radii-pos2))    # indice 2 più vicino alla pos 2 sulla curva più alta
       
       # Creo 'temp': il vettore con gli indici dei tempi per Ts e ECE 
       # NO: un vettore dei valori più vicino alla posizione estrema dell'intervallo 
       # NO:  (per la diagn che ha il max tra i due min)
       temp = [iEce,i]
       indice1 = temp[ind]  # Prendo come indice quello per la diagnostica che ho selezionato come curva 'più alta'
       # Creo pippo: contiene i due valori agli estremi della diagnostica 'superiore', 
       # per vedere quale dei due è più alto all'istante di tempo in analisi
       pippo = [v1,v2] = [diagM[idt1,indice1],diagM[idt2,indice1]] # array con 
       #  
       siaM = np.max(pippo) # Seleziono il punto più alto tra i due, da scegliere come rho-up
       rhoup = siaM + 0.0001 # aggiungo un piccolo delta di sicurezza
       
       #############################
       # RHO down
       # Trovo il punto più vicino della diag non considerata prima
       # seleziono l'altra diagnostica: minimo tra i due massimi
       indx = np.argmin(([mEce,mTs]))   # Trovo il minimo tra i due minimi
       rm = rDiag[indx]                 # Le posiszioni della diagnostica di cui sopra
       diagm = diag[indx]               # seleziono la diagnostica 'piu bassa'
       # Trovo l'indice del punto con rho piu vicina al minimo dell'altra:
       # avendo assi temporali diversi devo secegliere indici diversi a secondo del caso    
       # ripentedo il criterio già usato prima,
       temp = [iEce,i]
       indice2 = temp[indx]
       id_pluto = np.nanargmin(abs(diagm[:,indice2]-M)) 
       # Creo 'pluto' che contiene il valore della diagnostica 'piu bassa' all'istante selezionato, più
       # vicino al valore del massimo tra i due minimi:M
       pluto = diagm[id_pluto,indice2]
       
       rm = rDiag[indx]             # selezionop la diagnostica 'più bassa'
       pos_pluto = rm[id_pluto]     # Posizioni radiali della diagostica 'più bassa'
       
       # Se la posizione del punto è prima del minimo aggiungo 'numero_punti' punti e determino la rho down
       # se invece è dopo, ne tolgo 'numero_punti'
       # in ambedue i casi tolgo poi un piccolo margine
       if pos_pluto < pos[ind]:
           rhodown = diagm[id_pluto+numero_punti,indice2] - 0.001
       else:
           rhodown = diagm[id_pluto-numero_punti,indice2] - 0.001
       ranges[i,:] = [rhodown,rhoup]   
       # if id_pluto < indT:
       #     rhodown = diagm[id_pluto+numero_punti,indice2] - 0.001
       # else:
       #     rhodown = diagm[id_pluto-numero_punti,indice2] - 0.001
       # ranges[i,:] = [rhodown,rhoup]       
       
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

   ratio = temp_eceM_rho/temp_tsM_rho     # Ratio between temperaure values (averaged)
   distance = temp_eceM_rho-temp_tsM_rho  # Difference between temperaure values (averaged)
   err_dist = np.sqrt((err_eceM_rho/2)**2+(err_tsM_rho/2)**2)
   err_dist_perc = (1/temp_eceM_rho)*(np.sqrt(temp_tsM_rho/temp_eceM_rho*err_eceM_rho)**2+(err_tsM_rho)**2)
   
   print('Dimensioni temp_eceM_rho = ',temp_eceM_rho.shape)
   print('Dimensioni temp_TsM_rho = ',temp_tsM_rho.shape)    
   
   vars['temp_tsM_rho'] = temp_tsM_rho
   vars['err_tsM_rho'] = err_tsM_rho
   vars['temp_eceM_rho'] = temp_eceM_rho
   vars['err_eceM_rho'] = err_eceM_rho
   vars['timeTs2'] = timeTs2
   vars['timeEce22'] = timeEce22
   vars['ranges'] = ranges
   vars['ratio'] = ratio
   vars['distance'] = distance
   vars['err_dist'] = err_dist
   vars['err_dist_perc'] = err_dist_perc
   # Si trova la matrice 'ranges' con tutti i valoori di rhoup/down a tutti gli istanti del TS
   # Si trovano qiundi le due matrici: temp_eceM_rho e temp_tsM_rho 
   # che contengono i valori di Te a tutte le rho e a tutte gli istanti tra t1 e t2
   # assi temporali: timeTs2 e timeEce2
   # le posizioni: rTs e rEce
   
   # Print dei valori degli intervalli sui quali si effetttua la media
   tlim = vars['tlim'] # tempo al quale si ha la max Te-max
   