import numpy as np
import inspect
from collections import namedtuple
import matplotlib.pyplot as plt
from textwrap import dedent
import argparse
import colonsubplots
from scipy.integrate import cumtrapz

class Eqdsk:
    pass
   
def DD(variables):
    fr_actual = inspect.currentframe()
    fr = fr_actual.f_back # The calling frame
    
    try:
        wout = {}
        for var in variables.split(','):
            var  = var.strip()

            if var in fr.f_locals:
                wout[var] = fr.f_locals[var]
            else:
                raise RuntimeError(f'Not existent variable: {var}') 
    finally:
        del fr
        del fr_actual
    NTout = namedtuple('D', list(wout.keys()))
    ntout = NTout(**wout)
    return  ntout


def read_case(line):
    #case = line[:48]
    #i3 = int(line[48:48+5])
    #nrbox = int(line[48+5:48+10])
    #nzbox = int(line[48+10:])
    s_line = line.split()
    nzbox = int(s_line[-1])
    nrbox = int(s_line[-2])
    case = " ".join(s_line[:-2])
    return case, nrbox, nzbox


def read_float_line(line):
    values = []
    n = len(line)//16
    for i in range(n):
        str_value = line[16*i: 16*(i+1)]
        if len(str_value) > 1:
            values.append(float(str_value))
    return values


def read_1d(fid, n):
    v = []
    for i in range(0, n, 5):
        v += read_float_line(next(fid))
    return np.array(v)


def read_2_int_values(fid):
    line = next(fid)
    s_line = line.split()
    return int(s_line[0]), int(s_line[1])


def eqdsk_read(filename):
    with open(filename) as fid:
        case, nrbox, nzbox = read_case(next(fid))
        rboxlen, zboxlen, r0exp, rboxleft, zboxmid = read_float_line(next(fid))
        raxis, zaxis, psiaxis, psiedge, b0exp = read_float_line(next(fid))
        current, psiaxis, zero, raxis, zero = read_float_line(next(fid))
        zaxis, *zero = read_float_line(next(fid))
        f_dia = read_1d(fid, nrbox)
        pressure = read_1d(fid, nrbox)
        ffprime = read_1d(fid, nrbox)
        pprime = read_1d(fid, nrbox)
        psi = read_1d(fid, nrbox*nzbox).reshape(nzbox, nrbox)
        q = read_1d(fid, nrbox)

        nedge, nlimiter = read_2_int_values(fid)
        rz_edge = read_1d(fid, 2*nedge)
        rz_limiter = read_1d(fid, 2*nlimiter)
        r_edge = rz_edge[::2]
        z_edge = rz_edge[1::2]
        r_limiter = rz_limiter[::2]
        z_limiter = rz_limiter[1::2]


    r = np.linspace(rboxleft, rboxleft + rboxlen, nrbox)
    z = np.linspace(zboxmid - zboxlen/2, zboxmid + zboxlen/2, nzbox)
    psil = np.linspace(psiaxis, psiedge, nrbox)
    phi = cumtrapz(q, psil, initial=0)
    rho_tor_norm = np.sqrt(phi/phi[-1]) 
    sigma_ip = np.sign(current)
    sigma_b0 = np.sign(b0exp)
    sigma_psi = np.sign(psiedge - psiaxis)
    sigma_bp = sigma_psi/sigma_ip
    sigma_rhothetaphi = np.sign(q[0])/sigma_ip/sigma_b0

    return DD('''case, nrbox, nzbox,
      rboxlen, zboxlen, r0exp, rboxleft, zboxmid,
      raxis, zaxis, psiaxis, psiedge, b0exp,
      current,
      f_dia, pressure, ffprime, pprime, psi, q,
      r, z, psil,
      sigma_ip, sigma_b0, sigma_psi, sigma_bp, sigma_rhothetaphi,
      r_edge, z_edge, r_limiter, z_limiter,
      phi, rho_tor_norm''')


def create_eqdsk_fig():
    fig, ax = colonsubplots.subplots(
         {'dims': (1, 1)},
         {'dims': (3, 2)},
         sharex=True, figsize=(9, 6)
      )
    fig.subplots_adjust(hspace=0, wspace=0.07, left=0.1)
    return fig, ax

def ploteqdsk(w):
    fig, ax = create_eqdsk_fig()
    first = True
    for ww in w:
        fig.suptitle(ww.case)

        ploteqdsk_psi(ww, ax[0], first=first)
        ploteqdsk_lin(ww, ax[1], first=first)
        first = False
    return fig, ax


def ploteqdsk_psi(w, ax, first=True):
    if first:
        linestyles = 'solid'
        marker = 'o'
        ls = '-'
    else:
        linestyles = 'dashed'
        marker = '+'
        ls = '--'
    rr,zz = np.meshgrid(w.r, w.z)
    ax.contour(w.r, w.z, w.psi, 40, linestyles=linestyles)
    ax.plot(w.r_edge, w.z_edge, 'k', ls=ls)
    ax.plot(w.r_limiter, w.z_limiter, ls=ls)
    ax.plot(w.raxis, w.zaxis, marker)
    ax.plot(rr,zz,',', color=[0.7,0.7,0.7], zorder=-10)
    ax.axis('scaled')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')


def ploteqdsk_lin(w, ax, first=True):
    if first:
        ls = '-'
    else:
        ls = '--'
    psin = np.linspace(0, 1, w.psil.size)
    ax[0,0].plot(psin, w.f_dia, ls)
    ax[0,0].set_ylabel('F dia')
    ax[1,0].plot(psin, w.pressure/1e3, ls)
    ax[1,0].set_ylabel('P (kPa)')
    ax[2,0].plot(psin, w.q, ls)
    ax[2,0].set_ylabel('q')
    
    ax[0,1].plot(psin, w.ffprime, ls)
    ax[0,1].set_ylabel("F F'")
    ax[1,1].plot(psin, w.pprime/1e3, ls)
    ax[1,1].set_ylabel("P' (kPa/m)")
    ax[2,1].axis('off')
    br0 = w.r0exp*w.b0exp
    ax[0,0].axhline(br0, color='k', ls='--')
    psi_approx_weber = 4e-7*np.pi * w.current / 2 * w.raxis
    psi_approx = psi_approx_weber/2/np.pi
    info = f"""
    nr, nz:   {w.nrbox} {w.nzbox}
    psi_axis: {w.psiaxis}
    psi_edge: {w.psiedge}
    delta psi: {w.psiedge-w.psiaxis}
    R0:       {w.r0exp}
    B0:       {w.b0exp}
    R0*B0:    {br0}
    Ip:       {w.current/1e3} kA
    Rax, Zax: {w.raxis:.3f}, {w.zaxis:.3f}
    psi appr: {psi_approx_weber:.3f} Wb
    psi appr: {psi_approx:.3f} Wb/rad
    """
    if first:
        ax[2,1].text(0, 0.98,dedent(info), va='top', ha='left', 
                     transform=ax[2,1].transAxes)
    ax[-1,0].set_xlabel('psi norm')

    for axx in ax[:,1]:
        axx.yaxis.set_label_position("right")
        axx.yaxis.tick_right()

    for axx in ax.ravel():
        axx.grid(True)
        axx.set_xlim([0,1])


def is_running_from_ipython():
    from IPython import get_ipython
    return get_ipython() is not None


if __name__ == "__main__":
    if not is_running_from_ipython():
        parser = argparse.ArgumentParser('Displaying EQDSK file')
        parser.add_argument('file', nargs='+')
        args = parser.parse_args()
        w = [eqdsk_read(filename) for filename in args.file]
        ploteqdsk(w)
        plt.show()
