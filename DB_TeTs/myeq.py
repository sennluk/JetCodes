#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:21:09 2024
@author: lsenni
Function to calculate equilibrium at Jet.
From Edmondo Giovannozzi
Inputs:
    shot, self (diagnostica/canale), dda (in base a quale equilibrio voglio la psi?)
"""
import numpy as np
import ppfeg 
from scipy.interpolate import RegularGridInterpolator

class Equilibrium:
    ftor: ppfeg.V2d
    psi2d: ppfeg.V2d
    fbnd: ppfeg.V2d
    faxs: ppfeg.V2d
    psir: ppfeg.V2d
    psiz: ppfeg.V2d

    def get_rho(self, t):
        psi, phi = self.ftor.get(t=t)
        rho = np.sqrt(phi/phi[-1])
        return psi, rho

    def get_psi(self, t, r, z):
        _,rr = self.psir.get(t=t)
        _,zz = self.psiz.get(t=t)
        p2d = self.psi2d.get(t=t)[1].reshape((33,33))
        _,psi_bnd = self.fbnd.get(t=t)
        _,psi_ax = self.faxs.get(t=t)
        psi_norm = (p2d - psi_ax)/(psi_bnd - psi_ax)
        rg_interpolator = RegularGridInterpolator((zz,rr), psi_norm)
       
        return rg_interpolator((z,r))

def load_equilibrium(shot):
    dda = 'efit'          # egio aveva EFTF
    ftor = ppfeg.jetdata(shot, dda, 'ftor')
    psi2d = ppfeg.jetdata(shot, dda, 'psi')
    psir = ppfeg.jetdata(shot, dda, 'psir')
    psiz = ppfeg.jetdata(shot, dda, 'psiz')
    fbnd = ppfeg.jetdata(shot, dda, 'fbnd')
    faxs = ppfeg.jetdata(shot, dda, 'faxs')
    return Equilibrium(ftor, psi2d, fbnd, faxs, psir, psiz)