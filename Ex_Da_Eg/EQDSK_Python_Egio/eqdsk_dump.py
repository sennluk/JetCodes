#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 12:24:14 2021

@author: egio
"""
import numpy as np
import argparse
import ppfeg
from dataclasses import dataclass
import v2d
from typing import Any

_EFIT_VARIABLES = ('dfdp', 'dpdp','f','p', 'faxs','fbnd',
                   'psi', 'psir', 'psiz','q','rbnd', 'zbnd',
                   'rmag','zmag','rlim','zlim', 'bvac','xip')

def float_format(v):
    return [f"{vv:16.9e}" for vv in v]

def join_float(*v):
    return "".join(float_format(v))

def five_float(v):
    v = float_format(v)
    line = []
    while v:
        line.append("".join(v[:5]))
        v = v[5:] 
    return line


def first_line(shot, dda, t, nw, nh):
    info = f"JPN: {shot} DDA: {dda}  t:{t:.6f}"
    line = f"{info:48s}{0:4d}{nw:4d}{nh:4d}"
    return line
    
def pos_info(v):
    return v[-1] - v[0], v[0], 0.5*(v[0] + v[-1])


def print_lines(fid, lines):
    for line in lines:
        print(line, file=fid)

@dataclass
class EqlSlice:
    shot: int
    dda: str
    user: str
    t: float
    
    nw: int
    nh: int
    
    dfdp: Any
    dpdp: Any
    f: Any
    p: Any
    
    faxs: Any
    fbnd: Any
    
    psi: Any
    psir: Any
    psiz: Any
    
    q: Any
    
    rbnd: Any
    zbnd: Any
    
    rmag: Any
    zmag: Any
    
    rlim: Any
    zlim: Any

    bvac: Any
    xip: Any

    def __post_init__(self):
        for nm in ['rmag', 'zmag', 'faxs', 'fbnd', 'bvac','xip']:
            setattr(self, nm, getattr(self, nm)[0])
        

    def dump_file(self, suffix, case, sign=1):
        filename = f"{suffix}_pulse_{self.shot}_time_{self.t:.3f}_{case}_{self.dda}_{self.user}.eqdsk"
        print(self.bvac)
        bvac = self.bvac*sign
        xip = self.xip*sign
        f = self.f*sign
        with open(filename,'wt') as fid:
            line = first_line(self.shot, self.dda, self.t, self.nw, self.nh)
            print(line, file=fid)
            rdim, rleft, rmid = pos_info(self.psir)
            zdim, zbottom, zmid = pos_info(self.psiz)
            rcentr = 2.96 # as stated in EFIT
            print(join_float(rdim, zdim, rcentr, rleft, zmid), file=fid)
            print(join_float(self.rmag, self.zmag, self.faxs, self.fbnd, bvac), file=fid)
            print(join_float(xip, self.faxs, 0.0, self.rmag, 0.0), file=fid)
            print(join_float(self.zmag, 0.0, self.fbnd, 0.0, 0.0), file=fid)
            
            print_lines(fid, five_float(f))
            print_lines(fid, five_float(self.p))
            print_lines(fid, five_float(self.dfdp))
            print_lines(fid, five_float(self.dpdp))
            print_lines(fid, five_float(self.psi.ravel()))
            print_lines(fid, five_float(self.q))
            nb = self.rbnd.size
            nl = self.rlim.size
            print(f"{nb:5d}{nl:5d}", file=fid)
            rzbnd = np.vstack((self.rbnd, self.zbnd)).T.ravel()
            print_lines(fid, five_float(rzbnd))
            rzlim = np.vstack((self.rlim, self.zlim)).T.ravel()
            print_lines(fid, five_float(rzlim))
            
            

            
@dataclass
class EqlInfo:
    
    shot: int
    dda: str
    user: str
    
    dfdp: Any
    dpdp: Any
    f: Any
    p: Any
    
    faxs: Any
    fbnd: Any
    
    psi: Any
    psir: Any
    psiz: Any
    
    q: Any
    
    rbnd: Any
    zbnd: Any
    
    rmag: Any
    zmag: Any
    
    rlim: Any
    zlim: Any
    
    bvac: Any
    xip: Any
    
    def __init__(self, shot, dda='EFIT', user='JETPPF', seq=0):
        self.shot = shot
        self.dda = dda
        self.user = user
        self.seq = seq
        for nm in _EFIT_VARIABLES:
            setattr(self, nm, ppfeg.jetdata(shot, dda, nm, userid=user))

    def get_time(self):
        return self.dfdp.t
    
    def time_slice(self, t):
        v = {}
        for nm in _EFIT_VARIABLES:
            psi_lin, v[nm] = getattr(self, nm).get(t=t)
        nw = self.psir.r.size
        nh = self.psiz.r.size
        v['psi'] = v['psi'].reshape((nh,nw))
        v['shot'] = self.shot
        v['t'] = t
        v['dda'] = self.dda
        v['user'] = self.user
        v['nw'] = nw
        v['nh'] = nh
        v['dfdp'] = v['dfdp']* 4e-7 * np.pi
        return EqlSlice(**v)
        
def is_running_from_ipython():
    from IPython import get_ipython
    return get_ipython() is not None

if __name__ == "__main__":
   if not is_running_from_ipython():
      parser = argparse.ArgumentParser('Writing EQDSK file form equilibrium PPF')
      parser.add_argument('shot', type=int)
      parser.add_argument('time', nargs='*', type=float)
      parser.add_argument('-dda', default='EFIT')
      parser.add_argument('-exp', default='bae')
      parser.add_argument('-case', default='00')
      parser.add_argument('-sign', default=-1, type=int)
      parser.add_argument('-tlim', default='40,70')
      parser.add_argument('-user', default='JETPPF')
      args = parser.parse_args()
      eqlinfo = EqlInfo(args.shot, dda=args.dda, user=args.user)
      
      if args.time:
          time = np.array(args.time)
      else:
          time = eqlinfo.get_time()
      tlim = [float(t) for t in args.tlim.split(',')]
      idt = (time >= tlim[0]) & (time<= tlim[1])
      if np.any(idt):
          time = time[idt]

      for t in time:
          print("Saving time:", t)
          eqlinfo.time_slice(t).dump_file(args.exp,args.case,sign=args.sign)
          