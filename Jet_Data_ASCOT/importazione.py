# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd

file_t = '/home/lsenni/Python_LS/Jet_Data_ASCOT/T_E_104522.csv'
file_n = '/home/lsenni/Python_LS/Jet_Data_ASCOT/N_E_104522.csv'

dft = pd.read_csv(file_t, delimiter=',', decimal='.')
dfn = pd.read_csv(file_n, delmiter=',', decimal = '.')

dfpsi = dft[1:101]
dft = dfpsi
dfn = dft[1:101]




