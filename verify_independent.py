#!/usr/bin/env python3
"""
verify_independent.py - independent numerical verification (Table 4).

Written separately from the MATLAB program and using a different method at
every step:
  * core mode      : ratio form  u J1(u)/J0(u) = w K1(w)/K0(w), Brent's method,
                     plus the Rudolph-Neumann approximation [21]
  * cladding modes : 4x4 boundary-condition determinant, Brent's method
  * coupling       : adaptive quadrature (scipy.integrate.quad), evanescent
                     tail outside the cladding INCLUDED in the normalisation
  * orthogonality  : normalised overlap of different cladding modes, 0..inf
  * band depth     : Tmin = cos^2(kappa L) with kappa recomputed at lambda_r
  * resonances     : direct root of the phase-matching condition
  * SRI response   : chain rule  d(lambda_r)/dn3 = lambda_r*gamma*(-dn_cl/dn3)/dn_eff,
                     integrated over each SRI interval (Gauss-Legendre) and
                     compared with the shift computed by the MATLAB program

Reads the CSV files written by run_all.m in ./results and writes
results/table4_verification.csv. Requires numpy and scipy.
Usage:  python verify_independent.py  [results_dir]
"""
import sys, os, csv, time
import numpy as np
from scipy.special import jv, yv, kve
from scipy.optimize import brentq
from scipy.integrate import quad

RES = sys.argv[1] if len(sys.argv) > 1 else "results"

def read_csv(name):
    with open(os.path.join(RES, name)) as f:
        rows = list(csv.reader(f))
    return rows[0], np.array([[float(x) for x in r] for r in rows[1:]])

_, prm = read_csv("params.csv")
a1, a2, n1, n2, n3, dn, L, Lam, M, lam0, disp = prm[0]
M = int(M)
if disp:
    sys.exit("params.csv has disp = 1; verification is written for the constant-index model.")

# ---------------------------------------------------------------- core mode
def core_neff(lam):
    k = 2*np.pi/lam
    def f(n):
        u = k*a1*np.sqrt(n1**2 - n**2); w = k*a1*np.sqrt(n**2 - n2**2)
        return u*jv(1, u)/jv(0, u) - w*kve(1, w)/kve(0, w)
    lo = max(n2, np.sqrt(max(n1**2 - (2.4048/(k*a1))**2, 0))) + 1e-12
    return brentq(f, lo, n1 - 1e-12, xtol=1e-15, rtol=1e-15)

def rudolph_neumann(lam):
    V = 2*np.pi/lam*a1*np.sqrt(n1**2 - n2**2)
    b = (1.1428 - 0.996/V)**2
    return np.sqrt(n2**2 + b*(n1**2 - n2**2))

# ------------------------------------------------------- cladding modes (4x4)
def _mats(X, lam, s3):
    k = 2*np.pi/lam; u2 = X/a2
    n = np.sqrt(n2**2 - (u2/k)**2)
    u1 = k*np.sqrt(n1**2 - n**2); w = k*np.sqrt(n**2 - s3**2)
    z = np.zeros_like(X)
    # unknowns: core J0 | cladding J0, Y0 | outside K0(wr)/K0(wa2)
    Mx = np.array([
        [jv(0, u1*a1),       -jv(0, u2*a1),      -yv(0, u2*a1),      z],
        [-u1*jv(1, u1*a1),    u2*jv(1, u2*a1),    u2*yv(1, u2*a1),   z],
        [z,                   jv(0, u2*a2),       yv(0, u2*a2),     -np.ones_like(X)],
        [z,                  -u2*jv(1, u2*a2),   -u2*yv(1, u2*a2),   w*kve(1, w*a2)/kve(0, w*a2)]])
    return Mx, n, u1, u2, w

def det4(X, lam, s3):
    Mx = _mats(np.atleast_1d(X), lam, s3)[0]
    return np.linalg.det(np.moveaxis(Mx, (0, 1), (-2, -1)))

def clad_modes(lam, s3, Mm):
    k = 2*np.pi/lam
    if s3 >= n2:
        return np.full(Mm, np.nan), [None]*Mm
    Xmax = k*a2*np.sqrt(n2**2 - s3**2) - 1e-7
    Xs = np.linspace(1e-4, min(Xmax, np.pi*(Mm + 3)), 60*(Mm + 3))
    v = det4(Xs, lam, s3)
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    ns, fields = np.full(Mm, np.nan), [None]*Mm
    for j in range(min(Mm, len(idx))):
        X = brentq(lambda x: det4(x, lam, s3)[0], Xs[idx[j]], Xs[idx[j]+1], xtol=1e-14, rtol=1e-15)
        Mx, n, u1, u2, w = _mats(np.array([X]), lam, s3)
        _, _, vh = np.linalg.svd(Mx[:, :, 0]); c = vh[-1]
        ns[j] = n[0]
        fields[j] = (c, u1[0], u2[0], w[0])
    return ns, fields

def clad_field(fd, r):
    c, u1, u2, w = fd
    if r <= a1: return c[0]*jv(0, u1*r)
    if r <= a2: return c[1]*jv(0, u2*r) + c[2]*yv(0, u2*r)
    return c[3]*kve(0, w*r)/kve(0, w*a2)*np.exp(-w*(r - a2))

def core_field(lam, nco, r):
    k = 2*np.pi/lam; u = k*np.sqrt(n1**2 - nco**2); w = k*np.sqrt(nco**2 - n2**2)
    if r <= a1: return jv(0, u*r)/jv(0, u*a1)
    return kve(0, w*r)/kve(0, w*a1)*np.exp(-w*(r - a1))

def qint(f, a, b, pieces=1):
    edges = np.linspace(a, b, pieces + 1)
    return sum(quad(f, edges[i], edges[i+1], limit=400, epsabs=0, epsrel=1e-11)[0] for i in range(pieces))

def norm_clad(fd):
    w = fd[3]
    return (qint(lambda r: clad_field(fd, r)**2*r, 0, a1)
            + qint(lambda r: clad_field(fd, r)**2*r, a1, a2, 60)
            + qint(lambda r: clad_field(fd, r)**2*r, a2, a2 + 60/w))

def kappa_over_dn(lam, nco, fd):
    Nco = qint(lambda r: core_field(lam, nco, r)**2*r, 0, a1) + \
          qint(lambda r: core_field(lam, nco, r)**2*r, a1, a2, 20)
    ov = qint(lambda r: core_field(lam, nco, r)*clad_field(fd, r)*r, 0, a1, 4)
    return abs(np.pi/lam*ov/np.sqrt(Nco*norm_clad(fd)))

# --------------------------------------------------------- phase matching
def pm_res(lam, m, s3):
    return (core_neff(lam) - clad_modes(lam, s3, m)[0][m-1])*Lam - lam

def resonance(m, s3, lo=0.78, hi=1.22, n=23):
    g = np.linspace(lo, hi, n); v = np.array([pm_res(x, m, s3) for x in g])
    i = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(i) == 0: return np.nan
    return brentq(lambda x: pm_res(x, m, s3), g[i[0]], g[i[0]+1], xtol=1e-12)

def chain_rule_S(m, s3, lr):
    """d(lambda_r)/dn3 = lambda_r*gamma*(-dn_cl/dn3)/dn_eff  (nm/RIU)"""
    h, hn = 1e-5, 1e-6
    ncl = lambda lam, s: clad_modes(lam, s, m)[0][m-1]
    dneff = lambda lam: core_neff(lam) - ncl(lam, s3)
    D = dneff(lr)
    dD_dlam = (dneff(lr + h) - dneff(lr - h))/(2*h)
    gamma = 1.0/(1.0 - lr/D*dD_dlam)
    dncl_dn3 = (ncl(lr, s3 + hn) - ncl(lr, s3 - hn))/(2*hn)
    return lr*gamma*(-dncl_dn3)/D*1e3

# ===================================================================== checks
t0 = time.time(); rows = []
def report(item, method, value, crit, unit=""):
    ok = value < crit
    rows.append([item, method, value, crit, "PASS" if ok else "FAIL"])
    print(f"{item:<38s} {value:11.3e} {unit:<6s} (criterion < {crit:g})  {'PASS' if ok else 'FAIL'}")

# 1) core neff
_, cm = read_csv("core_mode.csv")
lam_c, nco_mat = cm[0][0], cm[0][1]
nco_py = core_neff(lam_c)
report("neff LP01 vs ratio form", "Brent, ratio form", abs(nco_py - nco_mat), 1e-9)
report("neff LP01 vs Rudolph-Neumann", "approximation [21]", abs(rudolph_neumann(lam_c) - nco_mat), 1e-5)

# 2) cladding neff, 3) kappa/dn, 4) orthogonality at lam0
_, cp = read_csv("coupling_at_lam0.csv")
ncl_py, fds = clad_modes(lam0, n3, M)
ok = ~np.isnan(cp[:, 1]) & ~np.isnan(ncl_py)
report(f"neff cladding nu=1-{2*M-1}", "4x4 determinant", np.max(np.abs(ncl_py[ok] - cp[ok, 1])), 1e-6)
K = 8   # nu = 1..15 (adaptive quadrature is slow; the modes used in the paper)
kd = np.array([kappa_over_dn(lam0, nco_py, fds[j]) for j in range(K)])
report("kappa/dn nu=1-15", "adaptive quadrature", 100*np.max(np.abs(kd - cp[:K, 3])/kd), 1.0, "%")
N = [norm_clad(fds[j]) for j in range(K)]
ortho = 0.0
for i in range(K):
    for j in range(i+1, K):
        o = (qint(lambda r: clad_field(fds[i], r)*clad_field(fds[j], r)*r, 0, a1)
             + qint(lambda r: clad_field(fds[i], r)*clad_field(fds[j], r)*r, a1, a2, 60)
             + qint(lambda r: clad_field(fds[i], r)*clad_field(fds[j], r)*r, a2, a2 + 60/fds[i][3]))
        ortho = max(ortho, abs(o)/np.sqrt(N[i]*N[j]))
report("orthogonality nu=1-15", "normalised overlap", ortho, 1e-4)

# 5) band depth, 6) resonant wavelengths
_, rs = read_csv("resonances.csv")
# depth from an independently computed kappa at the resonant wavelength
# (adaptive quadrature, evanescent tail included): Tmin = cos^2(kappa L)
dd = []
for nu, lr, kap, kL, tmin in rs[:, :5]:
    m = int((nu + 1)//2); lam_r = lr/1e3
    ncl_r, fd_r = clad_modes(lam_r, n3, m)
    kL_py = kappa_over_dn(lam_r, core_neff(lam_r), fd_r[m-1])*dn*L
    d_py = 10*np.log10(max(np.cos(kL_py)**2, 1e-6))
    dd.append((int(nu), abs(d_py - tmin), abs(kL_py - kL)))
    print(f"   nu={int(nu):2d}: kappa*L MATLAB {kL:.5f}  Python {kL_py:.5f}   Tmin MATLAB {tmin:8.2f} dB  Python {d_py:8.2f} dB")
report("kappa*L at resonance nu=1-15", "adaptive quadrature at lambda_r", max(d[2] for d in dd), 1e-3, "rad")
report("band depth nu=1-13", "cos^2(kappa L), independent kappa", max(d[1] for d in dd if d[0] <= 13), 0.05, "dB")
deep = [d for d in dd if d[0] > 13]
if deep:
    print(f"   (info) nu=15 depth difference {deep[0][1]:.2f} dB: a band deeper than 40 dB depends on kappa*L to 1e-4 rad")
dl = []
for nu, lr in rs[:, :2]:
    m = int((nu + 1)//2); lp = resonance(m, n3)
    dl.append(abs(lp*1e3 - lr))
    print(f"   nu={int(nu):2d}: MATLAB {lr:9.3f} nm   Python {lp*1e3:9.3f} nm")
report("resonant wavelength nu=1-15", "direct phase-matching root", max(dl), 0.1, "nm")

# 7) SRI response via chain rule
hdr, sw = read_csv("sri_sweep.csv")
modes = [5, 6, 7, 8]; gl_x, gl_w = np.polynomial.legendre.leggauss(4)
worst_abs, worst_rel = 0.0, 0.0
print("   SRI interval   nu   MATLAB d(lambda) [nm]   chain rule [nm]   diff [%]")
for s in range(len(sw) - 1):
    sa, sb = sw[s, 0], sw[s+1, 0]
    for j, m in enumerate(modes):
        if np.isnan(sw[s, j+1]) or np.isnan(sw[s+1, j+1]): continue
        d_mat = sw[s+1, j+1] - sw[s, j+1]
        nodes = 0.5*(sb - sa)*gl_x + 0.5*(sa + sb)
        Ss = [chain_rule_S(m, x, resonance(m, x, lo=(sw[s+1, j+1]-5)/1e3, hi=(sw[s, j+1]+5)/1e3, n=5))
              for x in nodes]
        d_chain = 0.5*(sb - sa)*np.dot(gl_w, Ss)          # nm
        rel = 100*abs(d_chain - d_mat)/max(abs(d_mat), 1e-9)
        worst_rel = max(worst_rel, rel) if abs(d_mat) > 0.05 else worst_rel
        worst_abs = max(worst_abs, abs(d_chain - d_mat))
        print(f"   {sa:.2f}-{sb:.2f}   {2*m-1:3d}   {d_mat:12.4f}          {d_chain:12.4f}     {rel:6.2f}")
report("SRI shift per interval (chain rule)", "chain rule, Gauss-Legendre", worst_rel, 3.0, "%")

with open(os.path.join(RES, "table4_verification.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["item", "independent_method", "max_difference", "criterion", "result"])
    w.writerows(rows)
print(f"\nVerification done in {time.time()-t0:.0f} s -> {RES}/table4_verification.csv")
