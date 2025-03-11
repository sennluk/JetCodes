#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 12:38:41 2024

@author: egio
"""
import numpy as np
from io import StringIO
import pandas as pd

SHOTS = """
shot    t_start t_end   kg10_seq
102068	48.5    49.5    188
102069	48.5    50.0    182
102071	48.6    49.8    184
102072	48.0    50.5    168
102074	48.785  49.5    178
102076	49.5    50.1    153
102079	49.7    50.3    170
103481	45.0    46.5    0
103055	44.5    46.5    0
103057	46.5    47.5    180
102081	46.5    47.5    177
102081	48.5    49.5    177
102082	46.0    47.5    171
101861	47.0    48.0    142
101861	48.5    49.5    142
103894	47.5    49.0    146
102068	48.5    49.5    188
102069	48.5    50.0    182
103705	48.0    48.7    208
102238	47.2    48.0    224
102752	47.5    49.0    113
102754	47.5    50.0    119
103705	47.5    48.8    208
103708	47.5    49.3    208
103800	48.0    49.0    157
103800	49.5    50.5    157
"""

# 102071 PKK3 283

def _separate(line):
    shot = int(line[0])
    v1,v2 = line[1].split('-')
    t_start = float(v1)
    t_end = float(v2)
    return shot, t_start, t_end

#_shot_info = [_separate(line.split('\t')) for line in SHOTS.split('\n') if line]

#shot_info = np.array(_shot_info, dtype=[('shot',int),('t_start', float),('t_end', float)]).view(np.recarray)

shot_io = StringIO(SHOTS)
shot_info_pd = pd.read_csv(shot_io, delim_whitespace=True)
shot_info = shot_info_pd.to_records(index=False)

