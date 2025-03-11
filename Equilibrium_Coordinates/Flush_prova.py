#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 13:42:29 2024

@author: lsenni
"""

import numpy as np
import py_flush as Flush

def main():
    #****Initialisation****
    shot = 73569
    time = 50.0
    time, ier = Flush.flushquickinit(shot,time)
    
    nup = 2
    r = np.zeros(nup)
    z = np.zeros(nup)
    r[0] = 300.0
    z[0] = 0.0
    z[1] = -10.0
    r[1] = 380.0
    z[1] = -10.0
    
    #****Call the routine****
    f,ier = Flush.Flush_getFlux(nup,r,z)
    
    #****Print output****
    print("Flux at first point", f[0])
    print("Flux at second point", f[1])

main()
