from __future__ import print_function

import numpy as np
from scipy.signal import hilbert, resample, medfilt, lfilter, butter, buttord, filtfilt
from scipy.integrate import cumtrapz
import inspect
import matplotlib.pyplot as plt

def _v(obj):
    if isinstance(obj, V2d):
        return obj.v
    else:
        return obj

class V2d:
    def __init__(self, r=None, t=None, v=None, shot=0, userid='', units='', xunits='', tunits='', 
                type=None, xtype=None, ttype=None, desc='',  ier = 0, additional_arguments=None):
        self.shot = shot
        self.userid = userid
        self.units = units
        self.xunits = xunits
        self.tunits = tunits
        self.type = type
        self.xtype = xtype
        self.ttype = ttype
        self.desc = desc
        self.r = r
        self.t = t
        self.v = v
        self.ier = ier
        self.additional_arguments = additional_arguments
     
    def copy(self):
        """return a copy of the given object
        """
        res = self._getacopy()
        return res        

    def _getacopy(self, copydata=True, inplace=False):
        """Return a copy of the given object.
        copydata: if True also the field v is copied
        inplace:  if True the actual object is reeturned isntead of a copy
        """
        if inplace:
            return self
            
        r = self.r.copy() if self.r is not None else None
        t = self.t.copy() if self.t is not None else None
        v = self.v.copy() if copydata and self.v is not None else None
        other = self.__class__()
        other.r = r
        other.t = t
        other.v = v
        other.shot = self.shot
        other.userid = self.userid
        other.units = self.units
        other.xunits = self.xunits
        other.tunits = self.tunits
        other.type = self.type
        other.xtype = self.xtype
        other.ttype = self.ttype
        other.desc =  self.desc
        other.ier =  self.ier
        
        if hasattr(self, 'additional_arguments'):
            other.additional_arguments = self.additional_arguments
        else:
            other.additional_arguments = None
        return other

    def __repr__(self):
        nr = len(self.r)
        nt = len(self.t)
        rstr = []
        rstr += ["    v: Float %dx%d" % (nr,nt)]
        rstr += ["    t: Float  %d" % nt]
        rstr += ["    r: Float  %d" % nr]
        rstr += ["  ier:        {}".format(self.ier)]
        return "\n".join(rstr)


    def ok(self):
        """Return True if the object contains correct data
        """
        return self.ier == 0 and self.v.size > 0

    def get(self, r=np.nan, t=np.nan, return_tr_found=False):
        if (np.isnan(t) and np.isnan(r)):
            return (self.r, self.t, self.v)
        elif (np.isnan(t)):
            ir = np.argmin(abs(self.r - r))
            if return_tr_found:
                return (self.t, self.v[ir,:], self.r[ir])
            else:
                return (self.t, self.v[ir,:])
        elif (np.isnan(r)):
            it = np.argmin(abs(self.t - t))
            if return_tr_found:
                return (self.r, np.atleast_2d(self.v)[:,it], self.t[it])
            else:
                return (self.r, np.atleast_2d(self.v)[:,it])
        else:
            ir = np.argmin(abs(self.r - r))
            it = np.argmin(abs(self.t - t))
            if return_tr_found:
                return np.atleast_2d(self.v)[ir,it], self.r[ir], self.t[it]
            else:
                return np.atleast_2d(self.v)[ir,it]

    def slice(self,  r=np.nan, t=np.nan):
        result = self._getacopy(copydata=False)
        if not self.ok():
            return result
        if (np.isnan(t) and np.isnan(r)):
            return result
        elif (np.isnan(t)):
            ir = np.argmin(abs(self.r - r))
            result.r = self.r[ir][np.newaxis]
            result.v = self.v[ir,:][np.newaxis,:]
        elif (np.isnan(r)):
            it = np.argmin(abs(self.t - t)) 
            result.t = self.t[it][np.newaxis]
            result.v = self.v[:,it][:,np.newaxis]
        else:
            ir = np.argmin(abs(self.r - r))
            it = np.argmin(abs(self.t - t))
            result.r = self.r[ir][np.newaxis]
            result.t = self.t[it][np.newaxis]
            result.v = self.v[ir,it][np.newaxis,np.newaxis]
        return result

    def get_r_profile(self, t = np.nan):
        it = np.argmin(abs(self.t - t))
        return (self.r, self.v[:,it])

    def get_t_profile(self, r = np.nan):
        ir = np.argmin(abs(self.r - r))
        return (self.t, self.v[ir,:])

    def trange(self, tlim):
        """Restrict the time axis in the specified range
        """
        tlimm = np.asarray(tlim,dtype=float)
        if tlimm.size == 2:
            idt = (self.t >= np.min(tlimm)) & (self.t <= np.max(tlimm))
            self.t = self.t[idt]
            self.v = self.v[:,idt]
        return self
         
    def tresample(self,n):
        """Resample the data along the time axis.
        Calling:
            obj.tresample(n,inplace=[True,False])
            n specify the number of time points it will get 
        """
        v,t = resample(self.v, n, t=self.t, axis=1)
        self.v = v
        self.t = t
        return self

    def moveto(self, t):
        result = self._getacopy(copydata=True)
        result.t = self.t - t
        return result

    def clip(self, a_min, a_max):
        result = self._getacopy(copydata=False)
        result.v = np.clip(self.v, a_min, a_max)
        return result

    def max(self):
        if self.ok():
            return self.v.max()
        else:
            return np.finfo('d').min

    def min(self):
        if self.ok():
            return self.v.min()
        else:
            return np.finfo('d').max

    def mean(self, trange):
        if not self.ok():
            return None
        trange = np.asarray(trange, dtype=float)
        idt = (self.t >= np.min(trange)) & (self.t <= np.max(trange))
        if len(self.v.shape) == 2 and self.v.shape[0] == self.t.size:
            return np.mean(self.v[idt,...])
        elif len(self.v.shape) == 2 and self.v.shape[1] == self.t.size:
            return np.mean(self.v[:,idt])
        else:
            return np.mean(self.v[idt])

    def medfilt(self, nm=5):
        if self.v is None:
            return
        result = self._getacopy(copydata=False)
        if len(self.v.shape) == 1:
            result.v = medfilt(self.v, kernel_size=nm)
        elif len(self.v.shape) == 2:
            vv = np.zeros(self.v.T.shape)
            for v,sv in zip(vv,self.v.T):
                v[:] = medfilt(sv, kernel_size=nm)
            result.v = vv.T
        return result 
        
    # TODO it may have to be deleted
    # def _not_working_decimate(self, q):
    #     result = self._getacopy(copydata=False)
    #     result.v = decimate(self.v, q)
    #     result.t = decimate(self.t, q)
    #     return result

    def filter(self, Wn, ftype='low', symmetric=True):
        if self.v is None:
            return
        result = self._getacopy(copydata=False)
        
        dt = np.mean(np.diff(self.t))
        
        FNy2 = 1.0 / 2.0 / dt
        Wno = Wn / FNy2
        if ftype == 'low':
            bord =  buttord(Wno,Wno*1.7,1,10)
        elif ftype == 'high':
            bord =  buttord(Wno*0.6,Wno,1,10)
        else:
            raise  Exception('What are you doing')
            
        bf,af = butter(*bord,btype=ftype)
        if len(self.v.shape) == 1:
            if symmetric:
                result.v = filtfilt(bf,af, self.v)
            else:
                result.v = lfilter(bf, af, self.v)
        elif len(self.v.shape) == 2:
            result.v = lfilter(bf, af, self.v, axis=0)
        return result


    def dt(self):
        """Derivative along the t direction
        """
        result = self._getacopy(copydata=False)
        vdiff = self.v[:,1:] - self.v[:,:-1]
        tdiff = self.t[1:] - self.t[:-1]
        tmea = 0.5*(self.t[1:] + self.t[:-1])
        result.t = tmea
        result.v = vdiff/tdiff
        return result

    def integrate(self, Wn=None, symmetric=True):
        if self.v is None:
            return self
        if not self.ok():
            return self
        result = self._getacopy(copydata=False)
        
        if len(self.v.shape) == 1:
            ii = cumtrapz(self.v, x=self.t)
            result.v = np.hstack((np.zeros(1),ii))
        else:
            assert(False)
        if Wn is not None:
            result.filter(Wn,ftype='high', symmetric=symmetric)
        return result
        
    def decimate(self, q):
        result = self._getacopy(copydata=False)
        result.v = np.atleast_2d(self.v)[:,::q].copy()
        result.t = self.t[::q].copy()
        return result
        
    def angle(self, inplace=False):
        result = self._getacopy(copydata=False)
        result.v = np.angle(self.v)
        return result
        
    def unwrap(self, axis=0, inplace=False):
        """Convenience method for np.unwrap.
        Contrary to the unwrap convention, it is applied on the first
        dimension (that generally correspond to the time coordinate)
        Calling:
            obj.unwrap(axis=<axis>)
        """
        result = self._getacopy(copydata=False)
        result.v = np.unwrap(self.v, axis=axis)
        return result
    
            
    def hilbert(self, inplace=False):
        """return the hilbert transform of the signal.
        The original signal is left unchanged.
        """
        if self.v is None:
            return
        result = self._getacopy(copydata=False)
            
        if 'axis' in inspect.getargspec(hilbert)[0]:
            result.v = hilbert(self.v, axis=-1)
        else:
            if len(self.v.shape) == 1:
                result.v = hilbert(self.v)
            elif len(self.v.shape) == 2:
                vv = np.zeros(self.v.T.shape,dtype=np.complex)
                for v,sv in zip(vv,self.v.T):
                    v[:] = hilbert(sv)
                result.v = vv.T
        return result


    def vtline(self,*arg, **args):
        """Plot a vertical line for each time which is present 
        """
        for t in self.t:
            plt.axvline(t,*arg, **args)
            
    def plot(self, t = np.nan, r = np.nan, ax=None, func=None, **kargs):
        """Plot the data, various keyword are defined:
            t:  if not is NAN the time slice closest to the given time is plotted 
        """
        if func is None:
            func = lambda x:x
        if not self.ok():
            return
        if ax is None:
            ax = plt
        if (np.isnan(t) and np.isnan(r)):
            temp = ax.plot(self.t,func(self.v.T), **kargs)
        elif (np.isnan(t)):
            ir = np.argmin(abs(self.r - r))
            temp = ax.plot(self.t,func(self.v[ir,:]), **kargs)    
        elif (np.isnan(r)):
            it = np.argmin(abs(self.t - t))
            temp = ax.plot(self.r,func(self.v[:,it]), **kargs)
        return temp

    def specgram(self, ax=None, **kargs):
        if not self.ok():
            return
        if ax is None:
            ax=plt
        dt = np.mean(np.diff(self.t))
        tmin = self.t[0]
        tmax = self.t[-1]
        ax.specgram(self.v.ravel(), Fs=1.0/dt, xextent=(tmin, tmax), **kargs);
        
    def psd(self, ax=None, **kargs):
        if not self.ok():
            return
        if ax is None:
            ax=plt
        dt = np.mean(self.t)
        ax.psd(self.v.ravel(), Fs=1.0/dt,  **kargs);


    def pcolormesh(self, ax=None, **kwargs):
        if not self.ok():
            return
        if ax is None:
            ax=plt
        return ax.pcolormesh(self.t, self.r, self.v, **kwargs)
    
    def contour(self, *args, ax=None, **kwargs):
        if not self.ok():
            return
        if ax is None:
            ax=plt
        return ax.contour(self.t, self.r, self.v, *args, **kwargs)

    def contourf(self, *args, ax=None, **kwargs):
        if not self.ok():
            return
        if ax is None:
            ax=plt
        return ax.contourf(self.t, self.r, self.v, *args, **kwargs)
    
    def show(self):
         plt.show()
         
        
    # TODO it may be deleted
    # def basecopy(self):
    #
    #     result.units = self.units
    #     result.xunits = self.xunits,self.tunits,self.type,self.xtype,self.ttype,self.desc

    def __add__(self, other):
        result = self._getacopy(copydata=False)
        if isinstance(other, V2d):
            result.v = self.v + other.v
        else:
            result.v = self.v + other
        return result

    def __radd__(self, other):
        result = self._getacopy(copydata=False)
        result.v = self.v + other
        return result

    def __iadd__(self, other):
        self.v = self.v + other
        return self
        
    def __sub__(self, other):
        result = self._getacopy(copydata=False)
        if isinstance(other, V2d):
            result.v = self.v - other.v
        else:
            result.v = self.v - other
        return result

    def __rsub__(self, other):
        result = self._getacopy(copydata=False)
        result.v = other - self.v
        return result

    def __isub__(self, other):
        self.v = self.v - _v(other)
        return self
        
    def __mul__(self, other):
        result = self._getacopy(copydata=False)
        if isinstance(other, V2d):
            result.v = self.v * other.v
        else:
            result.v = self.v * other
        return result
    
    def __rmul__(self, other):
        result = self._getacopy(copydata=False)
        result.v = other * self.v
        return result
        
    def __imul__(self, other):
        self.v = self.v * _v(other)
        return self

    def __div__(self, other):
        result = self._getacopy(copydata=False)
        if isinstance(other, V2d):
            result.v = self.v / other.v
        else:
            result.v = self.v / other
        return result
    
    __truediv__ = __div__ # python 3 corrispondence
    
    def __rdiv__(self, other):
        result = self._getacopy(copydata=False)
        result.v = other / self.v
        return result
    
    __rtruediv__ = __rdiv__ # python 3 corrispondence
    
    def __idiv__(self, other):
        self.v = self.v / _v(other)
        return self
        
    __itruediv__ = __idiv__ # python 3 corrispondence

    def __pow__(self, other):
        result = self._getacopy(copydata=False)
        result.v = self.v**other
        return result

    def __ipow__(self, other):
        self.v = self.v ** other
        return self
        
    def __neg__(self, other):
        result = self._getacopy(copydata=False)
        result.v = -self.v
        return result

    def __pos__(self, other):
        return self

    def abs(self):
        result = self._getacopy(copydata=False)
        result.v = np.abs(self.v)
        return result
    
    def __call__(self, *args, **kwargs):
        self.additional_arguments = (args, kwargs)
        return self

    def sqrt(self):
        result = self._getacopy(copydata=False)
        result.v = np.sqrt(self.v)
        return result

    def exp(self):
        result = self._getacopy(copydata=False)
        result.v = np.exp(self.v)
        return result

    def log(self):
        result = self._getacopy(copydata=False)
        result.v = np.log(self.v)
        return result

