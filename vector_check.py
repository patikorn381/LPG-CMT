#!/usr/bin/env python3
"""
vector_check.py - how accurate is the scalar (LP) approximation for this fibre?

The MATLAB program uses the scalar three-layer model. This script solves the
EXACT (full-vector) three-layer problem for the hybrid modes of azimuthal
order l = 1 (HE1m / EH1m, Erdogan's cladding-mode numbers nu = 1, 2, 3, ...)
and compares it with the MATLAB results:

  0. self-tests of the vector solver against the textbook exact
     characteristic equation of a two-layer step-index fibre
       (a) weakly guiding core mode HE11,
       (b) a bare silica rod in air (the high-contrast case), 16 modes
  1. effective indices at the reference wavelength (odd nu <-> scalar modes)
  2. resonant wavelengths of the eight bands (Table 2)
  3. shift and sensitivity versus surrounding index (Table 5)
  4. period of the turning point of nu = 15

Method. Fields ~ exp(i(l*phi + beta*z - w*t)), Ez = e(r), Hz = (i/Z0) g(r).
In every layer e and g are combinations of J_l, Y_l (n > neff) or I_l, K_l
(n < neff). Continuity of e, g, E_phi and H_phi at r = a1 and r = a2 gives a
real 8x8 homogeneous system; its determinant vanishes at the modes.

Reads the CSV files written by run_all.m in ./results and writes
results/table_vector_modes.csv, table_vector_sri.csv and vector_summary.csv.
Requires numpy and scipy.   Usage:  python vector_check.py [results_dir]
"""
import sys, os, csv, time
import numpy as np
from scipy.special import jv, yv, ive, kve, jvp, yvp
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline

RES = sys.argv[1] if len(sys.argv) > 1 else "results"
L_AZ = 1


def _fun(kind, q, r, r_ref):
    """value and d/dr of radial function; q = k*sqrt(|n^2-neff^2|). Scaled relative to r_ref for I/K."""
    x = q*r
    if kind == 'J':
        return jv(L_AZ, x), q*jvp(L_AZ, x)
    if kind == 'Y':
        return yv(L_AZ, x), q*yvp(L_AZ, x)
    if kind == 'I':   # I_l(x) scaled by exp(-q*r_ref)
        s = np.exp(x - q*r_ref)
        v = ive(L_AZ, x)*s
        dv = q*0.5*(ive(L_AZ-1, x) + ive(L_AZ+1, x))*s
        return v, dv
    if kind == 'K':   # K_l(x) scaled by exp(+q*r_ref)
        s = np.exp(-(x - q*r_ref))
        v = kve(L_AZ, x)*s
        dv = -q*0.5*(kve(L_AZ-1, x) + kve(L_AZ+1, x))*s
        return v, dv
    raise ValueError

def _rows(neff, k, n, kinds, r, r_ref):
    """For a layer with index n and list of radial-function kinds, return the 4 x (2*len(kinds)) block
    of [e, g, Ephi-like, Hphi-like] at radius r. Columns: e-coefficients then g-coefficients."""
    q2 = n*n - neff*neff
    q = k*np.sqrt(abs(q2)); Q2 = k*k*q2
    m = len(kinds); B = np.zeros((4, 2*m))
    for j, kd in enumerate(kinds):
        v, dv = _fun(kd, q, r, r_ref)
        # e column
        B[0, j] = v
        B[2, j] = k*neff*L_AZ*v/r/Q2
        B[3, j] = k*n*n*dv/Q2
        # g column
        B[1, m+j] = v
        B[2, m+j] = -k*dv/Q2
        B[3, m+j] = -k*neff*L_AZ*v/r/Q2
    return B

def matrix(neff, lam, a1, a2, n1, n2, n3):
    k = 2*np.pi/lam
    k1 = ['J'] if neff < n1 else ['I']
    k2 = ['J', 'Y'] if neff < n2 else ['I', 'K']
    k3 = ['K']
    M = np.zeros((8, 8))
    # unknown order: core e,g (2) | clad e(2), g(2) | out e,g (2)
    C1 = _rows(neff, k, n1, k1, a1, a1); D1 = _rows(neff, k, n2, k2, a1, a1 if neff >= n2 else a1)
    M[0:4, 0:2] = C1; M[0:4, 2:6] = -D1
    D2 = _rows(neff, k, n2, k2, a2, a1 if neff >= n2 else a1); C3 = _rows(neff, k, n3, k3, a2, a2)
    M[4:8, 2:6] = D2; M[4:8, 6:8] = -C3
    return M

def det(neff, *p):
    M = matrix(neff, *p)
    # row/column equilibration so the sign is robust
    M = M/np.abs(M).max(axis=1, keepdims=True)
    M = M/np.abs(M).max(axis=0, keepdims=True)
    return np.linalg.det(M)

def v_clad_scan(lam, a1, a2, n1, n2, n3, nmodes, pts_per_pi=120):
    """First nmodes l=1 cladding modes (n3 < neff < n2), ordered by decreasing neff: nu = 1,2,3,..."""
    k = 2*np.pi/lam; p = (lam, a1, a2, n1, n2, n3)
    Xmax = min(k*a2*np.sqrt(n2*n2 - n3*n3) - 1e-6, np.pi*(nmodes/2 + 4))
    Xs = np.linspace(1e-3, Xmax, int(pts_per_pi*Xmax/np.pi))
    ns = np.sqrt(n2*n2 - (Xs/(k*a2))**2)
    v = np.array([det(x, *p) for x in ns])
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    out = []
    for i in idx[:nmodes]:
        out.append(brentq(lambda x: det(x, *p), ns[i+1], ns[i], xtol=1e-15, rtol=1e-15))
    return np.array(out)

def v_core_scan(lam, a1, a2, n1, n2, n3):
    p = (lam, a1, a2, n1, n2, n3)
    ns = np.linspace(n1 - 1e-9, n2 + 1e-9, 2000)
    v = np.array([det(x, *p) for x in ns])
    i = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0][0]
    return brentq(lambda x: det(x, *p), ns[i+1], ns[i], xtol=1e-15, rtol=1e-15)

# ---------- classical exact two-layer equation (independent of the 8x8 matrix) ----------
def two_layer_exact(lam, a, nco, ncl, nmodes, npts=200000):
    """Roots of the textbook HE/EH (l=1) characteristic equation of a step-index rod:
       (J'/(uJ) + K'/(wK)) (nco^2 J'/(uJ) + ncl^2 K'/(wK)) = (l*neff)^2 (1/u^2 + 1/w^2)^2 """
    k = 2*np.pi/lam
    def F(neff):
        u = k*a*np.sqrt(nco**2 - neff**2); w = k*a*np.sqrt(neff**2 - ncl**2)
        J = jv(1, u); Jp = jvp(1, u)
        Kp_over_K = -0.5*(kve(0, w) + kve(2, w))/kve(1, w)
        # multiply through by (u J)^2 to remove poles
        t1 = (Jp/u + J*Kp_over_K/w)
        t2 = (nco**2*Jp/u + ncl**2*J*Kp_over_K/w)
        return t1*t2 - (L_AZ*neff)**2*(1/u**2 + 1/w**2)**2*J*J
    ns = np.linspace(nco - 1e-10, ncl + 1e-10, npts)
    v = F(ns); idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    return np.array([brentq(F, ns[i+1], ns[i], xtol=1e-15, rtol=1e-15) for i in idx[:nmodes]])

# ---------- scalar (LP) three-layer model, independent Python implementation ----------
def s_core(lam, a1, n1, n2):
    k = 2*np.pi/lam
    def f(n):
        u = k*a1*np.sqrt(n1**2 - n**2); w = k*a1*np.sqrt(n**2 - n2**2)
        return u*jv(1, u)*kve(0, w) - w*kve(1, w)*jv(0, u)
    ns = np.linspace(n1-1e-9, n2+1e-9, 4000); v = f(ns)
    i = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0][0]
    return brentq(f, ns[i+1], ns[i], xtol=1e-15, rtol=1e-15)
def _sf(X, lam, a1, a2, n1, n2, n3):
    k = 2*np.pi/lam; u2 = X/a2; n = np.sqrt(n2**2 - (u2/k)**2)
    u1 = k*np.sqrt(n1**2 - n**2); w = k*np.sqrt(n**2 - n3**2)
    J0c = jv(0,u1*a1); J1c = jv(1,u1*a1); J0a = jv(0,u2*a1); J1a = jv(1,u2*a1); Y0a = yv(0,u2*a1); Y1a = yv(1,u2*a1)
    d = u2*(J0a*Y1a - Y0a*J1a)
    A = (J0c*u2*Y1a - Y0a*u1*J1c)/d; B = (J0a*u1*J1c - u2*J1a*J0c)/d
    psi = A*jv(0,u2*a2) + B*yv(0,u2*a2); dpsi = -u2*(A*jv(1,u2*a2) + B*yv(1,u2*a2))
    return dpsi*kve(0,w*a2) + w*kve(1,w*a2)*psi
def s_clad(lam, a1, a2, n1, n2, n3, M):
    k = 2*np.pi/lam; Xmax = k*a2*np.sqrt(n2**2 - n3**2) - 1e-6
    Xs = np.linspace(1e-3, min(Xmax, np.pi*(M+3)), 80*(M+3)); v = _sf(Xs, lam, a1, a2, n1, n2, n3)
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]; out = np.full(M, np.nan)
    for j in range(min(M, len(idx))):
        X = brentq(lambda x: _sf(x, lam, a1, a2, n1, n2, n3), Xs[idx[j]], Xs[idx[j]+1], xtol=1e-14, rtol=1e-15)
        out[j] = np.sqrt(n2**2 - (X/(k*a2))**2)
    return out

# ======================================================================
def read_csv(name):
    with open(os.path.join(RES, name)) as f:
        rows = list(csv.reader(f))
    return rows[0], np.array([[float(x) for x in r] for r in rows[1:]])

def write_csv(name, header, rows):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)

_, prm = read_csv("params.csv")
a1, a2, n1, n2, n3, dn, L, Lam, M, lam0, disp = prm[0]
if disp:
    sys.exit("params.csv has disp = 1; this check is written for the constant-index model.")
t0 = time.time(); summary = []
def report(item, value, crit, unit=""):
    ok = value < crit
    summary.append([item, value, crit, "PASS" if ok else "FAIL"])
    print(f"{item:<58s} {value:10.3e} {unit:<4s} (criterion < {crit:g})  {'PASS' if ok else 'FAIL'}")
def info(item, value, unit=""):
    summary.append([item, value, "", "info"])
    print(f"{item:<58s} {value:10.4g} {unit}")

# vector modes found in a narrow bracket around the scalar value (fast);
# the full scan (v_clad_scan / v_core_scan) is used to validate the bracket.
def v_clad(lam, s3, nsc):
    p = (lam, a1, a2, n1, n2, s3)
    return brentq(lambda x: det(x, *p), nsc - 2e-5, nsc + 5e-6, xtol=1e-15, rtol=1e-15)
def v_core(lam, s3, nsc):
    p = (lam, a1, a2, n1, n2, s3)
    return brentq(lambda x: det(x, *p), nsc - 3e-5, nsc + 3e-5, xtol=1e-15, rtol=1e-15)
def residual(lam, m, s3, vector):
    nc = s_core(lam, a1, n1, n2); ncl = s_clad(lam, a1, a2, n1, n2, s3, m)[m-1]
    if vector:
        nc = v_core(lam, s3, nc); ncl = v_clad(lam, s3, ncl)
    return (nc - ncl)*Lam - lam
def root_near(f, guess_um, half=0.012, n=9):
    g = np.linspace(guess_um - half, guess_um + half, n); v = np.array([f(x) for x in g])
    i = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(i) == 0:
        return np.nan
    return brentq(f, g[i[0]], g[i[0]+1], xtol=1e-11)*1e3

print("== 0. self-tests of the vector solver ==")
hc = two_layer_exact(lam0, a1, n1, n2, 1)[0]
h8 = v_core_scan(lam0, a1, a2, n1, n2, n3)
report("HE11 core mode: 8x8 system vs textbook equation", abs(hc - h8), 1e-9)
rod = two_layer_exact(lam0, a2, n2, 1.0, 16, npts=2000000)
m8 = v_clad_scan(lam0, a1, a2, n2 + 1e-9, n2, 1.0, 16)
report("silica rod in air, 16 modes: 8x8 vs textbook equation", np.max(np.abs(rod - m8)), 1e-9)

print("\n== 1. effective indices at %.2f um, n3 = %.2f ==" % (lam0, n3))
_, cp = read_csv("coupling_at_lam0.csv")
K = 8
vs = v_clad_scan(lam0, a1, a2, n1, n2, n3, 2*K)
sc = s_clad(lam0, a1, a2, n1, n2, n3, K)
vb = np.array([v_clad(lam0, n3, sc[j]) for j in range(K)])
report("bracketed vector roots equal the scanned odd-nu roots", np.max(np.abs(vb - vs[0::2])), 1e-12)
report("MATLAB scalar neff vs this script's scalar neff (nu 1-15)", np.max(np.abs(sc - cp[:K, 1])), 1e-8)
nco_s = s_core(lam0, a1, n1, n2); nco_v = v_core(lam0, n3, nco_s)
print("   core mode: scalar LP01 %.8f   vector HE11 %.8f   difference %+.2e" % (nco_s, nco_v, nco_v - nco_s))
print("   nu | scalar (MATLAB) | vector HE (odd nu) | difference | vector EH (nu+1)")
rows = []
for j in range(K):
    print("   %2d | %.8f | %.8f | %+.2e | %.8f" % (2*j+1, cp[j, 1], vs[2*j], vs[2*j] - cp[j, 1], vs[2*j+1]))
info("max |neff(vector) - neff(scalar)|, cladding nu 1-15", np.max(np.abs(vs[0::2] - cp[:K, 1])))
info("|neff(vector) - neff(scalar)|, core mode", abs(nco_v - nco_s))

print("\n== 2. resonant wavelengths at Lambda = %g um, n3 = %.2f ==" % (Lam, n3))
_, rs = read_csv("resonances.csv")
print("   nu | scalar MATLAB (nm) | scalar Python (nm) | vector (nm) | vector - scalar (nm) | (%)")
for nu, lr in rs[:, :2]:
    m = int((nu + 1)//2)
    ls = root_near(lambda x: residual(x, m, n3, False), lr/1e3)
    lv = root_near(lambda x: residual(x, m, n3, True), lr/1e3)
    rows.append([int(nu), lr, ls, lv, lv - lr, 100*(lv - lr)/lr])
    print("   %2d | %9.2f | %9.2f | %9.2f | %+6.2f | %+.3f" % tuple(rows[-1]))
rows = np.array(rows)
write_csv("table_vector_modes.csv", ["nu", "lr_scalar_matlab_nm", "lr_scalar_python_nm", "lr_vector_nm", "diff_nm", "diff_percent"], rows.tolist())
report("scalar resonances: MATLAB vs Python", np.nanmax(np.abs(rows[:, 2] - rows[:, 1])), 0.1, "nm")
info("max |lambda_r(vector) - lambda_r(scalar)|", np.nanmax(np.abs(rows[:, 4])), "nm")
info("max relative difference of lambda_r", np.nanmax(np.abs(rows[:, 5])), "%")

print("\n== 3. response to the surrounding index ==")
hdr, sw = read_csv("sri_sweep.csv")
modes = [5, 6, 7, 8]; V = np.full((len(sw), 4), np.nan)
for i in range(len(sw)):
    for j, m in enumerate(modes):
        if not np.isnan(sw[i, j+1]):
            V[i, j] = root_near(lambda x: residual(x, m, sw[i, 0], True), sw[i, j+1]/1e3, half=0.02, n=13)
Ss = sw[:, 1:] - sw[0, 1:]; Vs = V - V[0]
print("   SRI  | shift scalar (nm): nu 9, 11, 13, 15      | shift vector (nm)")
for i in range(len(sw)):
    print("   %.2f | %s | %s" % (sw[i, 0], " ".join("%8.2f" % x for x in Ss[i]), " ".join("%8.2f" % x for x in Vs[i])))
dS = np.diff(sw[:, 1:], axis=0)/np.diff(sw[:, 0])[:, None]; dV = np.diff(V, axis=0)/np.diff(sw[:, 0])[:, None]
mid = 0.5*(sw[1:, 0] + sw[:-1, 0])
print("   sensitivity (nm/RIU), scalar | vector | vector/scalar")
for i in range(len(mid)):
    print("   %.3f | %s | %s | %s" % (mid[i], " ".join("%9.2f" % x for x in dS[i]), " ".join("%9.2f" % x for x in dV[i]),
                                      " ".join("%5.2f" % x for x in dV[i]/dS[i])))
write_csv("table_vector_sri.csv",
          ["SRI", "lr_s_nu9", "lr_s_nu11", "lr_s_nu13", "lr_s_nu15", "lr_v_nu9", "lr_v_nu11", "lr_v_nu13", "lr_v_nu15",
           "shift_s_nu9", "shift_s_nu11", "shift_s_nu13", "shift_s_nu15", "shift_v_nu9", "shift_v_nu11", "shift_v_nu13", "shift_v_nu15"],
          np.hstack([sw[:, :1], sw[:, 1:], V, Ss, Vs]).tolist())
write_csv("table_vector_sensitivity.csv",
          ["SRI_mid", "S_s_nu9", "S_s_nu11", "S_s_nu13", "S_s_nu15", "S_v_nu9", "S_v_nu11", "S_v_nu13", "S_v_nu15"],
          np.hstack([mid[:, None], dS, dV]).tolist())
info("sensitivity ratio vector/scalar, first SRI interval (max)", np.nanmax(dV[0]/dS[0]))
info("sensitivity ratio vector/scalar, last SRI interval (max dev. from 1)", np.nanmax(np.abs(dV[-1]/dS[-1] - 1)))
info("shift at the last SRI value: max relative difference", 100*np.nanmax(np.abs(Vs[-1]/Ss[-1] - 1)), "%")
# operating range (|S| >= 10 nm/RIU) with the vector model
for j, m in enumerate(modes):
    a = np.abs(dV[:, j]); i1 = np.where(a >= 10)[0]
    if len(i1) == 0: lo = np.nan
    elif i1[0] == 0: lo = sw[0, 0]
    else:
        i = i1[0]; lo = np.interp(1.0, np.log10(a[i-1:i+1]), mid[i-1:i+1])
    info("vector model: |S| >= 10 nm/RIU from SRI (nu = %d)" % (2*m-1), lo)

print("\n== 4. turning point of nu = 15 ==")
ls = np.linspace(1.10, 1.50, 41)
def period(l, vector):
    nc = s_core(l, a1, n1, n2); ncl = s_clad(l, a1, a2, n1, n2, n3, 8)[7]
    if vector: nc = v_core(l, n3, nc); ncl = v_clad(l, n3, ncl)
    return l/(nc - ncl)
lf = np.linspace(1.10, 1.50, 40001)
for name, vec in (("scalar", False), ("vector", True)):
    cs = CubicSpline(ls, [period(l, vec) for l in ls])(lf); i = np.argmax(cs)
    info("turning point of nu = 15, %s: period" % name, cs[i], "um (at %.0f nm)" % (lf[i]*1e3))

write_csv("vector_summary.csv", ["item", "value", "criterion", "result"], summary)
print("\nVector check done in %.0f s -> %s/table_vector_modes.csv, table_vector_sri.csv, table_vector_sensitivity.csv, vector_summary.csv" % (time.time() - t0, RES))
