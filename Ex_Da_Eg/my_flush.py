"""
my_flush.py

My interface to flush routines


"""
import numpy as np
from numpy.ctypeslib import ndpointer
from ctypes import CDLL, c_int, c_uint, c_char_p, c_double, POINTER, create_string_buffer, byref, c_void_p

lib_flush = CDLL('libflush.so')
lib_flush.flushinit_.argtypes = [
                    POINTER(c_int),         # igo
                    POINTER(c_int),         # ishot
                    POINTER(c_double),      # time
                    POINTER(c_int),         # lunget
                    POINTER(c_int),         # iseq
                    c_char_p,               # uid
                    c_char_p,               # dda
                    POINTER(c_int),         # lunmsg
                    POINTER(c_int),         # ier
                    c_uint,                  # uid_len
                    c_uint]                  # dda_len
lib_flush.flushinit_.restype = None

# time, ier = flush.flushinit(15,npulse,time,0,0,uid,dda,0);
# flushinit (IGO,ISHOT,TIME,LUNGET,ISEQ,UID,DDA,LUNMSG,IER)
def flushinit(igo, ishot, time, lunget=0, iseq=0, uid='JETPPF', dda='EFIT', lunmsg=0):
    igo = c_int(igo)
    ishot = c_int(ishot)
    time = c_double(time)
    lunget = c_int(lunget)
    iseq = c_int(iseq)
    uid = create_string_buffer(uid.encode())
    dda = create_string_buffer(dda.encode())
    lunmsg = c_int(lunmsg)
    ier = c_int(0)
    uid_len = c_uint(len(uid))
    dda_len = c_uint(len(dda))
    
    lib_flush.flushinit_(
                byref(igo), 
                byref(ishot), 
                byref(time), 
                byref(lunget), 
                byref(iseq), 
                uid, 
                dda, 
                byref(lunmsg), 
                byref(ier), 
                uid_len, dda_len)
    
    return time.value, ier.value

# nbar, psibar, xpout, ypout, ier = flush.fluqax(1, np.atleast_1d(np.abs(float(M)/float(N))))
# void fluqax(int32_t *nQ, double *q, int32_t *nOut, double *psi, double *rInner, double *rOuter, int32_t *ierr)
    
lib_flush.fluqax_.argtypes = [
                    POINTER(c_int),                 # nQ
                    ndpointer(dtype=np.float64),    # q, 
                    POINTER(c_int),                 # nOut, 
                    ndpointer(dtype=np.float64),    # psi, 
                    ndpointer(dtype=np.float64),    # rInner, 
                    ndpointer(dtype=np.float64),    # rOuter, ierr)
                    POINTER(c_int)]                 # ierr
lib_flush.fluqax_.restype = None

def fluqax(nQ, q):
    q = np.ascontiguousarray(np.atleast_1d(q),dtype=np.float64)
    nn = q.size
    nQ = c_int(nn)
    nOut = c_int(nn)
    psi = np.empty(nn,dtype=np.float64)
    rInner = np.empty(nn,dtype=np.float64)
    rOuter = np.empty(nn,dtype=np.float64)
    ierr = c_int(0)
    lib_flush.fluqax_(
                byref(nQ),
                q,
                byref(nOut),
                psi,
                rInner,
                rOuter,
                byref(ierr)
                )
    return nOut.value, psi, rInner, rOuter, ierr.value

# Flush_getPsiGivenQ (np, q, psi, ier)
lib_flush.flush_getpsigivenq_.argtypes = [
                POINTER(c_int),                 # np
                ndpointer(dtype=np.float64),    # q
                ndpointer(dtype=np.float64),    # psi
                POINTER(c_int)]                 # ier

def Flush_getPsiGivenQ(q):
    q = np.ascontiguousarray(np.atleast_1d(q),dtype=np.float64)
    nn = q.size
    print("Flush q size", nn)
    
    psi = np.empty(nn, dtype=np.float64)
    
    nn = c_int(nn)
    ier = c_int(0)
    lib_flush.flush_getpsigivenq_(byref(nn), q, psi, byref(ier))
    return psi, ier.value


# xp, yp, br, bz, ier = flush.flusu2(1, np.atleast_1d(psibar), npoint, npoint, np.zeros(npoint), np.zeros(npoint,dtype=int), 2)
# void flusu2(int32_t *nPsi,double *psi,int32_t *np0,int32_t *np,double *r,double *z,double *br,double *bz,double *work,int32_t *jwork,int32_t *lopt,int32_t *ier)

lib_flush.flusu2_.argtypes = [
                    POINTER(c_int), #nPsi,
                    ndpointer(dtype=np.float64),  # psi,
                    POINTER(c_int),         # np0,
                    POINTER(c_int),         # np,
                    ndpointer(dtype=np.float64), # r,
                    ndpointer(dtype=np.float64), # z,
                    ndpointer(dtype=np.float64), # br,
                    ndpointer(dtype=np.float64), # bz,
                    ndpointer(dtype=np.float64), # work,
                    ndpointer(dtype=np.int32),   # jwork,
                    POINTER(c_int),              # lopt,
                    POINTER(c_int)]              # ier)
lib_flush.flusu2_.restype = None

def flusu2(nPsi, psi, npoint, npdim=0, work=None, jwork=None, lopt=2):
    psi = np.ascontiguousarray(np.atleast_1d(psi),dtype=np.float64)
    nn = psi.size
    
    nPsi = c_int(nn)
    npoint = c_int(npoint)
    npdim = c_int(npdim)
    
    nnpoint = max(npoint.value, npdim.value)
    r = np.zeros((nn,nnpoint))
    z = np.zeros((nn,nnpoint))
    br = np.zeros((nn,nnpoint))
    bz = np.zeros((nn,nnpoint))
    
    work = np.zeros(nnpoint)
    jwork = np.zeros(nnpoint, dtype=np.int32)
    
    lopt = c_int(lopt)
    ier = c_int(0)
    lib_flush.flusu2_(
                byref(nPsi),
                psi,
                byref(npoint),
                byref(npdim),
                r, z, br, bz, work, jwork,
                byref(lopt),
                byref(ier) )
    return r, z, br, bz, ier.value
                
# Flush_getClosedFluxSurfaces (nflux, flux, nPoints, r, z, ier)

lib_flush.Flush_getClosedFluxSurfaces.argtypes = [
                POINTER(c_int),                 # nflux
                ndpointer(dtype=np.float64),    # flux
                POINTER(c_int),                 # nPoints
                ndpointer(dtype=np.float64),    # r
                ndpointer(dtype=np.float64),    # z
                POINTER(c_int)]                 # ier)
def Flush_getClosedFluxSurface(flux,nPoints=360):
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    r = np.empty((nflux,nPoints))
    z = np.empty((nflux,nPoints))
    nflux = c_int(nflux)
    nPoints = c_int(nPoints)
    ier = c_int(0)
    lib_flush.Flush_getClosedFluxSurfaces(byref(nflux),flux, byref(nPoints), r, z, byref(ier))
    return r, z, ier.value

# No sane way to implement it (a memory leak is the minimum one can get)
# Flush_getFluxSurfaces(nflux, flux, accuracy, nPieces, nPoints, rSurf, zSurf, ier)
#lib_flush.Flush_getFluxSurfaces.argtypes = [
                    #POINTER(c_int),                 # nflux, 
                    #ndpointer(dtype=np.float64),    # flux, 
                    #POINTER(c_double),              # accuracy, 
                    #POINTER(c_int),                 # nPieces, 
                    #POINTER(c_int),                 # nPoints, 
                    #ndpointer(dtype=np.float64),    # rSurf, 
                    #ndpointer(dtype=np.float64),    # zSurf, 
                    #POINTER(c_int)]                 # ier

#-----------------------------------------------------------------------
# xp0, yp0, fp0, ier = flush.Flush_getMagAxis()
# Flush_getMagAxis (rAxis, zAxis, fAxis, ier)

# Flush_getLCFSboundary (nLCFS, rLCFS, zLCFS, accuracy, ier)
lib_flush.Flush_getLCFSboundary.argtypes = [
                    POINTER(c_int),                 # nLCFS
                    POINTER(c_void_p),              # rLCFS
                    POINTER(c_void_p),              # zLCFS
                    POINTER(c_double),              # accuracy
                    POINTER(c_int)                  # ier
                    ]
def Flush_getLCFSboundary(accuracy=1e-4):
    nLCFS = c_int(0)
    accuracy = c_double(accuracy)
    rLCFS = c_void_p(0)
    zLCFS = c_void_p(0)
    ier = c_int(0)
    print("Prima della chiamata")
    lib_flush.Flush_getLCFSboundary(byref(nLCFS), byref(rLCFS), byref(zLCFS), byref(accuracy), byref(ier))
    print("dopo la chiamata")
    return nLCFS.value, rLCFS, zLCFS, ier.value


lib_flush.flush_getmagaxis_.argtypes = [
                        POINTER(c_double),   # rAxis
                        POINTER(c_double),   # zAxis
                        POINTER(c_double),   # fAxis
                        POINTER(c_int) ]     # ier
def Flush_getMagAxis():
    rAxis = c_double(0.0)
    zAxis = c_double(0.0)
    fAxis = c_double(0.0)
    ier = c_int(0)
    lib_flush.flush_getmagaxis_(byref(rAxis), byref(zAxis), byref(fAxis), byref(ier))
    return rAxis.value, zAxis.value, fAxis.value, ier.value,

# fp, br, bz, bt, xpa, elong, ier = flush.flulax(2, npoint, xp, yp, np.zeros(npoint*5))

lib_flush.flush_geterror_.argtypes = [POINTER(c_int)]

def Flush_getError(errNum):
    errNum = c_int(errNum)
    lib_flush.flush_geterror_(errNum)



_some_routines = """
lib_flush.{1}.argtypes = [
                    POINTER(c_int),                 # np
                    ndpointer(dtype=np.float64),    # r
                    ndpointer(dtype=np.float64),    # z
                    ndpointer(dtype=np.float64),    # V
                    POINTER(c_int)]                 # ier

def {0}(r,z):
    r = np.ascontiguousarray(np.atleast_1d(r),dtype=np.float64)
    z = np.ascontiguousarray(np.atleast_1d(z),dtype=np.float64)
    assert r.shape == z.shape
    n = r.size
    V = np.empty(n)
    n = c_int(n)
    ier = c_int(0)
    lib_flush.{1}(byref(n), r, z, V, byref(ier))
    return V.reshape(r.shape), ier.value
"""
_globs = globals()
_locs = locals()

def deffun(name):
    name = "Flush_"+name
    lib_name = name.lower()+"_"
    exec(_some_routines.format(name,lib_name), _globs, _locs)

deffun("getBr")
deffun("getBz")
deffun("getBt")
deffun("getJr")
deffun("getJz")
deffun("getJt")
deffun("getFlux")
deffun("getAbsoluteFlux")
deffun("getdpsidr")
deffun("getdpsidz")
deffun("getdpsidrdr")
deffun("getdpsidzdz")
deffun("getdpsidrdz")
deffun("getdpsidzdr")
deffun("getdFdPsi")
deffun("getFdFdPsi")
deffun("getdPdPsi")


_routines_vs_psi = """
lib_flush.{1}.argtypes = [
                    POINTER(c_int),                 # np
                    ndpointer(dtype=np.float64),    # psi
                    ndpointer(dtype=np.float64),    # V
                    POINTER(c_int)]                 # ier

def {0}(psi):
    psi = np.ascontiguousarray(np.atleast_1d(psi),dtype=np.float64)
    n = psi.size
    V = np.empty(n)
    n = c_int(n)
    ier = c_int(0)
    lib_flush.{1}(byref(n), psi, V, byref(ier))
    return V.reshape(psi.shape), ier.value
"""
_globs = globals()
_locs = locals()

def deffunvspsi(name):
    name = "Flush_"+name
    lib_name = name.lower()+"_"
    exec(_routines_vs_psi.format(name,lib_name), _globs, _locs)

deffunvspsi("getPProfile")
deffunvspsi("getFtorProfile")
deffunvspsi("getqProfile")

# flush_getMidPlanProjRight(np, r, z, Rm, ier)
lib_flush.flush_getmidplaneprojright_.argtypes =[
                            POINTER(c_int),         # np
                            ndpointer(dtype=float), # r
                            ndpointer(dtype=float), # z
                            ndpointer(dtype=float), # rm
                            POINTER(c_int),         # ier
                            ]
lib_flush.flush_getmidplaneprojright_.restype = None
def Flush_getMidPlaneProjRight(r, z):
    r = np.ascontiguousarray(r, dtype=float)
    z = np.ascontiguousarray(z, dtype=float)
    n = r.size
    if z.size != n:
        raise RuntimeError(f"Incompatible dimension r[{n}], z[{z.size}]")
    Rm = np.zeros(n, dtype=float)
    n = c_int(n)
    ier = c_int(0)
    lib_flush.flush_getmidplaneprojright_(n, r, z, Rm, ier)
    return Rm
    

# flush_getminorradius(np, flux, minRad, ier)
lib_flush.flush_getminorradius_.argtypes = [
                            POINTER(c_int),                 # np
                            ndpointer(dtype=np.float64),    # flux
                            ndpointer(dtype=np.float64),    # minRad
                            POINTER(c_int)]                 # ier
def Flush_getMinorRadius(flux):
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    minRad = np.empty((nflux))
    ier = c_int(0)
    nflux = c_int(nflux)
    lib_flush.flush_getminorradius_(byref(nflux),flux, minRad, byref(ier))
    return minRad, ier.value
# flush_getminorradius(np, flux, minRad, ier)
lib_flush.flush_getnormalisedminorradius_.argtypes = [
                            POINTER(c_int),                 # np
                            ndpointer(dtype=np.float64),    # flux
                            ndpointer(dtype=np.float64),    # minRad
                            POINTER(c_int)]                 # ier

def Flush_getNormalisedMinorRadius(flux):
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    sh = flux.shape
    minRad = np.empty((nflux))
    ier = c_int(0)
    nflux = c_int(nflux)
    lib_flush.flush_getnormalisedminorradius_(byref(nflux),flux, minRad, byref(ier))
    return minRad.reshape(sh), ier.value

# flush_getfluxlabelrho(nflux, flux, rho, ier)
lib_flush.flush_getfluxlabelrho_.argtypes = [
                            POINTER(c_int),                 # np
                            ndpointer(dtype=np.float64),    # flux
                            ndpointer(dtype=np.float64),    # rho
                            POINTER(c_int)]                 # ier
def Flush_getRho(flux):
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    rho = np.empty((nflux))
    ier = c_int(0)
    nflux = c_int(nflux)
    lib_flush.flush_getfluxlabelrho_(byref(nflux),flux, rho, byref(ier))
    return rho, ier.value

# flush_readfirstwall(np, rlim, zlim, ier)

lib_flush.flush_readfirstwall_.argtypes = [
                            POINTER(c_int),                 # np
                            ndpointer(dtype=np.float64),    # rlim
                            ndpointer(dtype=np.float64),    # zlim
                            POINTER(c_int)]                 # ier
def Flush_readFirstWall():
    rLim = np.zeros(300)
    zLim = np.zeros(300)
    nn = c_int(0)
    ier = c_int(0)
    lib_flush.flush_readfirstwall_(
                byref(nn),
                rLim,
                zLim,
                byref(ier))
    if ier.value == 0:
        rLim = rLim[:nn.value]
        zLim = zLim[:nn.value]
    else:
        rLim = np.zeros(0)
        zLim = np.zeros(0)
    return rLim, zLim, ier.value
        
# Flush_getVolume (nflux, flux, volume, ier)
lib_flush.flush_getvolume_.argtypes = [
                            POINTER(c_int),      # nflux
                            ndpointer(dtype=np.float64),    # flux
                            ndpointer(dtype=np.float64),    # volume
                            POINTER(c_int)]                 # ier
def Flush_getVolume(flux):
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    volume_cm3 = np.empty((nflux))
    ier = c_int(0)
    nflux = c_int(nflux)
    lib_flush.flush_getvolume_(byref(nflux),flux, volume_cm3, byref(ier))
    volume_m3 = volume_cm3/1e6 
    return volume_m3, ier.value

#Flush_getIntersections
# (r, z, tht, acc, np, flux, nfound, r1,z1, r2,z2, r3,z3, r4,z4, ier)
lib_flush.flush_getintersections_.argtypes = [
                            POINTER(c_double), # r
                            POINTER(c_double), # z
                            POINTER(c_double), # tht
                            POINTER(c_double), # acc
                            POINTER(c_int),    # np
                            ndpointer(dtype=np.float64),  # flux
                            ndpointer(dtype=c_int),    # nfound
                            ndpointer(dtype=np.float64),  # r1
                            ndpointer(dtype=np.float64),  # z1
                            ndpointer(dtype=np.float64),  # r2
                            ndpointer(dtype=np.float64),  # z2
                            ndpointer(dtype=np.float64),  # r3
                            ndpointer(dtype=np.float64),  # z3
                            ndpointer(dtype=np.float64),  # r4
                            ndpointer(dtype=np.float64),  # z4
                            POINTER(c_int)       # ier
                            ]
def Flush_getIntersections(r, z, tht, acc, flux):
    r = c_double(r)
    z = c_double(z)
    tht = c_double(tht)
    acc = c_double(acc)
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    nfound = np.empty(nflux, dtype=c_int)
    r1 = np.empty(nflux)
    z1 = np.empty(nflux)
    r2 = np.empty(nflux)
    z2 = np.empty(nflux)
    r3 = np.empty(nflux)
    z3 = np.empty(nflux)
    r4 = np.empty(nflux)
    z4 = np.empty(nflux)
    nflux = c_int(nflux)
    ier = c_int(0)
    lib_flush.flush_getintersections_(
                byref(r), byref(z), byref(tht), byref(acc), 
                byref(nflux), flux, nfound, r1, z1, r2, z2, r3, z3, r4, z4, byref(ier))
    return nfound, r1, z1, r2, z2, r3, z3, r4, z4, ier.value
    
lib_flush.flush_writeeqdskfile_.argtypes = [
                        c_char_p,        #  filename
                        POINTER(c_int),  #  nR
                        POINTER(c_int),  #  nZ
                        POINTER(c_int),  #  ier
                        c_uint
                        ]

def Flush_writeEQDSKfile(filename, nR, nZ):
    filename = create_string_buffer(filename.encode())
    nR = c_int(nR)
    nZ = c_int(nZ)
    ier = c_int(0)
    filename_len = c_uint(len(filename))
    lib_flush.flush_writeeqdskfile_(
        filename,
        byref(nR),
        byref(nZ),
        byref(ier),
        filename_len)
    return ier.value
    
lib_flush.flush_getsplinedomainlimits_.argtypes = [ 
                        POINTER(c_double),    # rMin, 
                        POINTER(c_double),    # zMin, 
                        POINTER(c_double),    # rMax, 
                        POINTER(c_double),    # zMax, 
                        POINTER(c_int)]       # ier)
def Flush_getSplineDomainLimits():
    rMin = c_double(0.0)
    zMin = c_double(0.0)
    rMax = c_double(0.0)
    zMax = c_double(0.0)
    ier = c_int(0)
    lib_flush.flush_getsplinedomainlimits_(byref(rMin), byref(zMin), 
                                           byref(rMax), byref(zMax), byref(ier))
    return rMin.value, zMin.value, rMax.value, zMax.value, ier.value 

# Flush_getBref (rRef, bRef, ier)
lib_flush.flush_getbref_.argtypes = [
                        POINTER(c_double),    # rRef, 
                        POINTER(c_double),    # bRef, 
                        POINTER(c_int)]       # ier)

def Flush_getBref():
    rRef = c_double(0.0)
    bRef = c_double(0.0)
    ier = c_int(0)
    lib_flush.flush_getbref_(byref(rRef), byref(bRef), byref(ier))
    return rRef.value, bRef.value, ier.value
 
# Flush_getFluxAveragedQuantities
# (nflux, flux, nQuantities, quantitiesNames, quantities, ier)

lib_flush.flush_getfluxaveragedquantities_.argtypes = [
                        POINTER(c_int),                 # nflux
                        ndpointer(dtype=np.float64),    # flux
                        POINTER(c_int),                 # nQuantities
                        ndpointer(dtype=np.byte),       # quantitiesNames
                        ndpointer(dtype=np.float64),    # quantities
                        POINTER(c_int),                 # ier
                        c_uint]                         # N quantitiesNames
def Flush_getFluxAveragedQuantities(flux, quantitiesNames):
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    names = np.ascontiguousarray(np.atleast_1d(quantitiesNames),dtype=np.string_)
    nQuantities = names.size
    n_char = names.dtype.itemsize
    names = names.view(dtype=np.byte)
    quantities = np.empty((nQuantities,nflux),dtype=np.float64)
    nflux = c_int(nflux)
    nQuantities = c_int(nQuantities)
    ier = c_int(0)
    n_char = c_uint(n_char)
    lib_flush.flush_getfluxaveragedquantities_(
                        nflux, flux, nQuantities, names, quantities, ier, n_char)
    return quantities, ier.value

#Flush_getlcfsFlux (lcfsFlux, ier)
lib_flush.flush_getlcfsflux_.argtypes = [
                        POINTER(c_double),    # lcfsFlux, 
                        POINTER(c_int)]       # ier)
def Flush_getlcfsFlux():
    lcfsFlux = c_double(0.0)
    ier = c_int(0)
    lib_flush.flush_getlcfsflux_(lcfsFlux, ier)
    return lcfsFlux.value, ier.value
    
    
#Flush_getTangentsToSurfaces (r, z, angle, nflux, flux, iSide, rTan, zTan, accuracy, nBeams, ier)
lib_flush.flush_gettangentstosurfaces_.argtypes = [
                       POINTER(c_double), # r
                       POINTER(c_double), # z
                       ndpointer(dtype=np.float64), # angle(nflux)
                       POINTER(c_int),    # nflux
                       ndpointer(dtype=np.float64), # flux(nflux)
                       POINTER(c_int),    # iside
                       ndpointer(dtype=np.float64), # rTan
                       ndpointer(dtype=np.float64), # zTan
                       POINTER(c_double), # accuracy
                       POINTER(c_int),    # nBeams
                       POINTER(c_int)]    # ier

def Flush_getTangentsToSurfaces(r, z, flux, iside, accuracy, nBeams):
    r = c_double(r)
    z = c_double(z)
    flux = np.ascontiguousarray(np.atleast_1d(flux),dtype=np.float64)
    nflux = flux.size
    angle = np.empty(nflux,dtype=np.float64)
    rTan = np.empty(nflux,dtype=np.float64)
    zTan = np.empty(nflux,dtype=np.float64)
    nflux = c_int(nflux)
    accuracy = c_double(accuracy)
    nBeams = c_int(nBeams)
    iside = c_int(iside)
    ier = c_int(0)
    lib_flush.flush_gettangentstosurfaces_(
                r, z, angle, nflux, flux, 
                iside, rTan, zTan, accuracy, nBeams, ier)
    return angle, rTan, zTan, ier.value
    

