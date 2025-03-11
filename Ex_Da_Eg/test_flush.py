#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 15:14:48 2024

@author: egio
"""
import my_flush
import numpy as np
import ppfeg


shot = 95986
w = ppfeg.ppfs(shot)

zece = w.ecm1.antp.v[1,0]
prfl = w.ecm1.prfl
time = 47.0
r, te = prfl.get(t=time) 

#r = np.linspace(2.5, 3.8, 1000)
z = np.full_like(r, zece)


ts, ier = my_flush.flushinit(15, shot, time)
psi, ier = my_flush.Flush_getFlux(r*100, z*100)
