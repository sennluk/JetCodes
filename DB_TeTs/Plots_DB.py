#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:46:08 2024
@author: lsenni
Plot ed elaborazioni partendo dal file salvato con la creazione del Database

Ritrasformare in DataFrame ?
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


DB_Name = 'DB_Fish'

# Load the database
db = np.load(f'{DB_Name}.npy',allow_pickle='TRUE').item()

shot_list = db.keys()
print(shot_list)

shot = 96994    # selszionare lo shot per numero o con l'indie della lista

time = db[shot]['time']
Thrts = db[shot]['Tts']
Tece = db[shot]['Tkk1']

plt.figure()
plt.plot(time, Thrts)
plt.plot(time, Tece)

plt.figure()
plt.scatter(Tece, Thrts)