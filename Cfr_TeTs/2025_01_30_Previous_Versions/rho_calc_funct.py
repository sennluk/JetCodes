#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:38:51 2024

@author: lsenni
"""

import numpy as np
import my_flush
from ppfeg import ppfs
import matplotlib.pyplot as plt
from scipy import signal,interpolate

def rhocalc(shot, w, tlim, psiKk1, psi1, psi2):
    # Prove di calcolo e plots

    tTs = w.hrts.te       # Chan Te HRTS  
    psiTs = w.hrts.psi    # Channel psi hrts
    tEce = w.ecm1.prfl 
    
    timeEce = tEce.t
    timeTs = tTs.t
    
    rho_ece_ = np.zeros(tEce.v.shape)
    rho_ts_ = np.zeros(tTs.v.shape)
    
    for i,time in enumerate(timeEce):
        psi = psiKk1[:,i]
        fl_int = interpolate.make_interp_spline(w.efit.ftor.r, w.efit.ftor.v[:,i]/w.efit.ftor.v[-1,i]) # function to intertpolate the flux
        fl_int_ece = fl_int(psi)
        rho_ece_[:,i] = np.sqrt(fl_int_ece)
        
    
    for i,time in enumerate(timeTs):
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
    
    idx= np.argmin(abs(timeEce - tlim))
    rho_ece = rho_ece_[:,idx]
    
    idx = np.argmin(abs(timeTs-tlim))
    rho_ts = rho_ts_[:,idx]              
    
    return