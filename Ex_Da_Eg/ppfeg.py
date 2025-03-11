from __future__ import print_function
import numpy as np
import ppf
import getdat
try:
    import kc1m
except ImportError:
    from . import kc1m

#Disable warning when run under pytest
import warnings
warnings.filterwarnings("ignore")

import xarray as xr

from netCDF4 import Dataset
import os
from operator import itemgetter
try:
    from v2d import V2d
except ImportError:
    from .v2d import V2d


# global config
JETDATAEXCEPTION = False
_GLOBAL_SET = {'raise_exception':JETDATAEXCEPTION,
                'return_xarray':False}
                
def raise_exception_if_error(v):
    """Set to True if you want an exception when no data is found.
    Normally jetdata will return an empty object."""
    global JETDATAEXCEPTION
    JETDATAEXCEPTION = v
    _GLOBAL_SET['raise_exception'] = v
    
def return_xarray(v):
    """Set to True if you want jetdata return an xarray instead of a V2d object"""
    _GLOBAL_SET['return_xarray'] = v

def reset_flags():
    """Reset the global flag to their default. Don't raise any exception
    and return a V2d Object"""
    global JETDATAEXCEPTION
    JETDATAEXCEPTION = False
    _GLOBAL_SET['raise_exception'] = False
    _GLOBAL_SET['return_xarray'] = False

class JetDataNotFound(RuntimeError):
    pass



class jetppf:
    """Class to create JET PPFs:
    It can use a context manager:
    
    with jetppf(shot, userid) as w:
        w.add(dda, dtype_1, data_1, x_1, t_1, description1, [dunit=..., xunit=..., tunit=...])
        w.add(dda, dtype_2, data_2, x_2, t_2, description2, [dunit=..., xunit=..., tunit=...])
        ...
        
    The close is automatically called at the end by the context manager.
    w.seq after the end of the context manager will give you the sequence number.
    
    In the open a lot of default value (basically all set to zero) have been chosen. 
    Look at the init function.
    Write in order all the DTYPE of a given DDA, then switch to a new DDA. 
    No check is perfomed but remember that the PPF Library doesn't permit to reopen
    a closed DDA, so you will get an error. A DDA is automatically closed
    when you start writing to a new DDA. 
    """
    def __init__(self, shot, userid, com="PPF",status=0):
        """Open the PPF. The userid should be yours, you should have the permission to wite
        PPFs with that userid.
        """
        date = 0 #np.zeros(6,dtype=int)
        time = 0 #np.zeros(6,dtype=int)
        self.shot = shot
        self.userid = userid
        self.date = date
        self.time = time
        self.com = com
        self.status = status
        self.seq = 0
        self.channels = []
        self.completed = False
        ppf.ppfuid(self.userid,"W")
        ier=ppf.ppfopn(self.shot,self.date,self.time,self.com,self.status)
        if ier == 0:
            self.ok = True
        else:
            self.ok = False
            
    def add(self, dda, dtyp, data, x, t, desc, dunit=" ", xunit="m", tunit="s"):
        """Add a DDA/DTYPE to the open PPF file.
        It uses 'ppfwritedouble'. 
        According to the PPF library convention the fastest data coordinate 
        is the X one.
        That is the latest coordinate in the C ordering (like in the standard
        numpy ordering) or the first coordinate if the F ordering is used. 
        
        No actual check is made at the moment.
        
        """
        if self.ok and not self.completed:
            xx = np.asarray(x).flatten()
            tt = np.asarray(t).flatten()
            ddata = np.asarray(data).flatten()
            nx = len(xx)
            nt = len(tt)
            #ddata = data.flatten()
        
            irdat=ppf.ppfwritedouble_irdat(nx,nt)
            ihdat=ppf.ppfwritedouble_ihdat(dunit,xunit,tunit,"F","F","F",desc)
            iwdat,ier=ppf.ppfwritedouble(self.shot,dda,dtyp,irdat,ihdat,ddata,xx,tt)
            if ier != 0:
                print("there is an error:{0} {1}/{2}".format(ier, dda, dtyp))
            else:
                self.channels += ((dda, dtyp),)
        elif self.completed:
            print("The ppf is already closed")
            
    def close(self,program="JetPPF Python",vers=1):
        """Close the PPF and get the sequence number (stored in the filed 'seq')"""
        if self.ok:
            seq,ier=ppf.ppfclo(self.shot,program,vers)
            self.seq = seq
            if ier != 0:
                self.ok = False
            
            self.completed = True
            
    def __repr__(self):
        rstr = []
        rstr += ["Pulse:     %d" % self.shot]
        rstr += ["Ok:        %d" % self.ok]
        rstr += ["Completed: %d" % self.completed]
        rstr += ["Sequence:  %d" % self.seq]
        for ch in self.channels:
            rstr += ["  -> %4s / %4s" % ch] 
        return "\n".join(rstr)

    def __enter__(self):
        """Entering the contex manager (basically returning self)"""
        return self
    
    def __exit__(self, *args, **kwargs):
        """Exiting the contex manager (basically calling 'close')"""
        self.close()
        
    
class Transp:
    def __init__(self, shot, run):
        self.shot = shot
        self.run = run
        base_dir = "/home/pshare/transp-data/transp/result/JET/"
        full_name = os.path.join(base_dir, str(shot), run,run+'.CDF')
        self.td = Dataset(full_name)
    def printinlong(self, name):
        for k, v in self.td.variables.items():
            if name.upper() in v.long_name:
                print("  {:5} {:10} {}".format(k,str(v.shape),v.long_name.strip()))
    
    
    def __getattr__(self, name):
        try:
            v = getattr(self.td, name)
        except AttributeError:
            v = self.td[name.upper()]
        return v
    
    def __getitem__(self,name):
        return self.td[name]

# Single JET data used by the Fishbones module
class Kc1mJetData(V2d):
    def __init__(self, shot=None, channel='', tstart=40, tend=60, dt=1.0):
        if shot is None:
            return
        tm = kc1m.Kc1mTimingCollection(shot,[channel], tstart, tend)
        tb = tm.getdata(dt)
        ier = tm.ier
        super().__init__(r=np.array([0.0]), t=tb.gettimes(), v=tb.vdata[0,:], 
                shot=shot, userid=' ', units='T/s', xunits=' ', tunits='s', 
                type=' ', xtype=' ', ttype=' ', desc=channel,  ier = ier)
        if self.ier == 0:
            self.fs = 1.0/np.mean(np.diff(self.t))
            self.Fn = self.fs/2.0 # Nyquist frequency
        else:
            self.fs = 0.0
            self.Fn = 0.0
    
    def _getacopy(self, copydata=False):
        result = super()._getacopy(copydata=copydata)
        result.fs = self.fs
        result.Fn = self.Fn
        return result
        
def _jetdata_dispatcher(r=None, t=None, v=None, shot=0, userid='', units='', xunits='', tunits='', 
                type=None, xtype=None, ttype=None, desc='',  ier = 0, additional_arguments=None):
    """Internal function that given the values of the flags create a V2d or a xarray DataArray"""
    
    if _GLOBAL_SET['return_xarray']:
        result = _jetdata_xarray(r,t,v,shot,userid,units,xunits,tunits,type,xtype,ttype,desc,ier, additional_arguments)
    else:
        result = V2d(r,t,v,shot,userid,units,xunits,tunits,type,xtype,ttype,desc,ier, additional_arguments)
    
    return result 
    
def _jetdata_xarray(r=None, t=None, v=None, shot=0, userid='', units='', xunits='', tunits='', 
                type=None, xtype=None, ttype=None, desc='',  ier = 0, additional_arguments=None):
    """Internal function that that create an xarray DataArray, the input arguments are the same 
    as those required by a V2d. Small change the 'type' has been changed to 'dtype' to avoid
    clash with the standard python function."""
    result = xr.DataArray(data=v, coords=(('x',r),('t',t)), 
                attrs={'shot':shot,'userid':userid,'units':units,'xunits':xunits,'tunits':tunits,
                'dtype':type, 'xtype':xtype,'ttype':ttype,'desc':desc, 'ier':ier})
    return result
    

def transp_dataset(shot, filename):
    base_dir = "/home/pshare/transp-data/transp/result/JET/"
    full_name = os.path.join(base_dir, str(shot), filename,filename+'.CDF')
    return Dataset(full_name)

def jetdata(shot, dda='', dtype='',seq=0,userid="JETPPF", copydata=True, jpf=None, transp=None):
    """Retrieve a jet data. Can be use to get PPFs, JPFs and Transp data."""
    if (isinstance(shot, V2d)): # Copy constructor it is a relic, but it may be used some where
        other = shot # a better name
        if (copydata):
            v = other.v.copy()
        else:
            v = None
        try:
            additional_arguments = other.additional_arguments
        except AttributeError:
            additional_arguments = None
        
        result = _jetdata_dispatcher(other.r.copy(), other.t.copy(), v=v, shot=other, userid=other.userid, 
                units = other.units, xunits=other.xunits, tunits=other.tunits, 
                type=other.type, xtype=other.xtype, ttype=other.ttype, 
                desc=other.desc,  ier=other.ier, additional_arguments=additional_arguments)
    else:
        if jpf is not None:
            data,t,nwds,title,units,ier=getdat.getdat8(jpf,shot)
            userid = None
            
            xunits = ' '
            tunits = 's'
            type = ' '
            xtype = ' '
            ttype = ' '
            desc = title
            r = np.array([0.0])
            v = np.reshape(data,(len(r),len(t)),order='F')
        elif transp is not None:
            base_dir = "/home/pshare/transp-data/transp/result/JET/"
            filename, channel = transp.split('/')
            full_name = os.path.join(base_dir, str(shot), filename,filename+'.CDF')
            r = np.array([])
            t = np.array([])
            v = np.array([])
            xunits = ' '
            tunits = 's'
            type = ' '
            xtype = ' '
            ttype = ' '
            try:
                with Dataset(full_name) as dataset:
                    var_list = [key for key in dataset.variables.keys() if dataset.variables[key].shape]
                    if channel.upper() in var_list:
                        variable = dataset.variables[channel.upper()]
                        units = variable.units.strip();
                        v = variable[:].T
                        desc = variable.long_name.strip()
                        if len(variable.shape)==1:
                            time_var = variable.dimensions[0]
                            tunits = dataset.variables[time_var].units.strip()
                            t = dataset.variables[time_var][:]
                            r = np.array([])
                        elif len(variable.shape)==2:
                            time_var = variable.dimensions[0]
                            x_var = variable.dimensions[1]
                            t = dataset.variables[time_var][:]
                            tunits = dataset.variables[time_var].units.strip()
                            r = dataset.variables[x_var][:].T
                            xunits =  dataset.variables[x_var].units.strip()
                            if np.all(np.diff(r,axis=-1)==0,axis=-1).all():
                                r = r[:,0]
                        ier = 0
                    else:
                        ier = -1
            except FileNotFoundError:
                ier = -1
        else:
            if userid != 'JETPPF':
                ppf.ppfuid(userid, rw='R')
            if seq>0:
                ppf.ppfgo(shot, seq)
            else:
                ppf.ppfgo(0,0)
            ihdat,iwdat,v,r,t,ier = ppf.ppfreaddouble(shot, dda, dtype) #, seq, userid)
            if userid != 'JETPPF':
                ppf.ppfuid('JETPPF', rw='R')
            units, xunits, tunits, type, xtype, ttype, desc=ppf.ppfget_ihdat(ihdat)
            v = np.reshape(v,(len(r),len(t)),order='F')
            
        result = _jetdata_dispatcher(r, t, v=v, shot = shot, userid = userid, units=units,
                xunits=xunits, tunits=tunits, 
                type = type, xtype = xtype, ttype = ttype, 
                desc = desc, ier = ier, additional_arguments=None)
        if JETDATAEXCEPTION and not result.ok():
            if jpf is not None:
                chname = jpf
            else:
                chname = dda + '/' + dtype
            raise JetDataNotFound("Error reading Shot:{} data:{}".format(shot,chname))
    return result


#The next classes are used to navigate the ppfs tree
class ddas:
    """Get a list of a dtype given a DDA.
    the single dtype can be accessed as it was an attribute.
    """
    def __init__(self,pulse,dda, userid=None, sequence=0, actual_sequence=0):
        self.pulse = pulse
        self.dda = dda
        self.userid = userid
        self.sequence = sequence
        self.actual_sequence = actual_sequence

        if userid is not None:
            ppf.ppfuid(userid)

        if sequence>0:
            ppf.ppfgo(pulse, sequence)
        else:
            ppf.ppfgo(0,0)
        
        ndt,dtnames,lxtv,dtcomms,ier = ppf.pdtinf(pulse,dda)
        self.error = ier
        self.ndt = ndt
        if ier == 0:
            self.dtnames = [dtnm.strip()  for dtnm in dtnames]
            self.dtcomms = dtcomms
            self.lxtv = lxtv.reshape(ndt,2)
        else:
            self.dtnames = []
            self.dtcomms = []
            self.lxtv = np.zeros((0,2))

    def __repr__(self):
        """Show a list of dtype inside the DDA.
        It shows:
           name   nx  nt   comments
        """
        rstr =[]
        llist = sorted(zip(self.dtnames, self.lxtv[:,0], self.lxtv[:,1], self.dtcomms), key=lambda x: x[0])
        rstr += [" %4s:  %3d : %3d -- %s" % ll for ll in llist]
        return "\n".join(rstr)

    def __getattr__(self, name):
        if name.upper() in self.dtnames:
            return jetdata(self.pulse, dda=self.dda, dtype=name,seq=self.sequence,userid=self.userid)
        else:
            raise AttributeError("The requested DTYPE <%s> is not present" % name.upper())
    
    def __getitem__(self, name):
        return getattr(self, name)
        
    def __dir__(self):
        return sorted(self.dtnames)

class _ppfs:
    def __init__(self,pulse,userid='JETPPF', sequence=0):
        
        ppf.ppfuid(userid,"R")
        if sequence != 0:
            dup = 1
        else:
            dup = 0
        nseq,lseq,ldda,ndda,cdda,ier=ppf.pdlppf(pulse,uid=-1,dup=dup)
        self.error = ier
        self.pulse = pulse
        self.userid = userid
        self.sequence = sequence

        if ier == 0:
            assert ndda == len(cdda)
            assert nseq == lseq.size == ldda.size
            assert ldda.sum() == ndda
            dda_seq = []

            for i in range(nseq):
                dda_seq += ldda[i]*[lseq[i]]
            self.cdda = []
            self.dda_seq = []
            for dda,seq in zip(cdda, dda_seq):
                if sequence <= 0 or sequence == seq:
                    self.dda_seq.append(seq)
                    self.cdda.append(dda.strip())
        else:
            self.cdda = []
            self.dda_seq = []

        dda_seq_sorted = sorted(zip(self.cdda, self.dda_seq),key=itemgetter(1))
        dda_seq_sorted = sorted(dda_seq_sorted,key=itemgetter(0))
        if dda_seq_sorted:
            self.cdda, self.dda_seq = zip(*dda_seq_sorted)
        self.ddas_dict = dict()
    def __repr__(self):
        rstr =  ["pulse: %d" % self.pulse]
        rstr += ["userid: %s" % self.userid]
        if self.error != 0:
            rstr += ["Error no.: %d" % self.error]
        ndda = len(self.cdda)
        ncol = 4
        
        rstr_dda = ["   {:4}-{:3d}".format(dda,seq) for dda,seq in zip(self.cdda, self.dda_seq)]
        rstr_3dda = []
        for i in range(0,ndda,ncol):
            rstr_3dda += [" ".join(rstr_dda[i:i+ncol])]
        rstr += rstr_3dda
        #rstr += [" %s" % dda.lower() for dda in sorted(self.cdda)]
        return "\n".join(rstr)
        
    def __dir__(self):
        return sorted([nm.lower() for nm in self.cdda])
        
    def __getattr__(self,name):
        if self.error != 0:
            raise AttributeError("No ppfs for the requested user <%s>" % self.userid)
        if name.upper() in self.cdda:
            if name.upper() not in self.ddas_dict.keys():
                idnm = self.cdda.index(name.upper())
                actual_sequence = self.dda_seq[idnm]
                
                self.ddas_dict[name.upper()] =  ddas(self.pulse, name, 
                                                     userid=self.userid, 
                                                     sequence=self.sequence, 
                                                     actual_sequence=actual_sequence)
                
            return self.ddas_dict[name.upper()]

        else:
            raise AttributeError("The requested DDA <%s> is not present" % name.upper())

        
class ppfs:
    def __init__(self, pulse, userid='JETPPF', sequence=0):
        self.pulse = pulse
        self.base_userid = userid
        self.base_sequence = sequence
        self.ppfuser = {(userid,sequence):_ppfs(pulse, userid=userid,sequence=sequence)}

    def __repr__(self):
        return self.ppfuser[(self.base_userid,self.base_sequence)].__repr__()

    def __getattr__(self, name):
        return self.ppfuser[(self.base_userid,self.base_sequence)].__getattr__(name)

    def __dir__(self):
        return self.ppfuser[(self.base_userid,self.base_sequence)].__dir__()
        
    def __getitem__(self,sliced_userid):
        if isinstance(sliced_userid,slice):
            userid = sliced_userid.start
            sequence = sliced_userid.stop
        else:
            userid = sliced_userid
            sequence = None
        if userid is None:
            userid = self.base_userid
        if sequence is None:
            sequence = 0
            
        if (userid,sequence) not in self.ppfuser:
            value = _ppfs(self.pulse,userid=userid,sequence=sequence)
            self.ppfuser[(userid,sequence)] = value
        
        return self.ppfuser[(userid,sequence)]
    
    def __call__(self, jpf=''):
        return jetdata(self.pulse,jpf=jpf)
