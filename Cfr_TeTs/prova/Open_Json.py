#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 16:21:18 2025
@author: lsenni
Prova riapertura DB json 
"""


import json

## Per riaprire il dizionario
with open("/home/lsenni/Python_LS/Cfr_TeTs/prova/DB_23_Shots.json", "r") as f:
     DB = json.load(f)