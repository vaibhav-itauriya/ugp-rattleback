"""Equations of motion for a semi-ellipsoidal rattleback, rolling or slipping.

Independent Python implementation of the model in rattleback_simulator_123H.m.
State y = [w1 w2 w3 phi theta psi v1 v2 v3 xG yG]; v in body components.
Checked against the MATLAB code: reversal times 2.831 s (rolling) and 4.223 s
(Coulomb, mu = 0.15, eps = 3e-5 m/s), against 2.832 s and 4.224 s.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
from scipy.integrate import solve_ivp

m, g = 1.0, 9.81
a1, a2, a3, b = 0.2, 0.03, 0.02, 0.005
I = np.array([[2e-4, -5e-5, 0.0], [-5e-5, 1.6e-3, 0.0], [0.0, 0.0, 1.75e-3]])
Iinv = np.linalg.inv(I)
A2 = np.array([a1**2, a2**2, a3**2])

def rotmat(ang):
    p, t, s = ang
    R1 = np.array([[1,0,0],[0,np.cos(p),-np.sin(p)],[0,np.sin(p),np.cos(p)]])
    R2 = np.array([[np.cos(t),0,np.sin(t)],[0,1,0],[-np.sin(t),0,np.cos(t)]])
    R3 = np.array([[np.cos(s),-np.sin(s),0],[np.sin(s),np.cos(s),0],[0,0,1]])
    return R1 @ R2 @ R3

def geom(ang, w):
    phi, th, ps = ang
    n = np.array([np.sin(phi)*np.sin(ps) - np.cos(ps)*np.sin(th)*np.cos(phi),
                  np.cos(ps)*np.sin(phi) + np.sin(ps)*np.sin(th)*np.cos(phi),
                  np.cos(th)*np.cos(phi)])
    nd = -np.cross(w, n)
    S = A2 @ (n*n); sq = np.sqrt(S)
    r = -A2*n/sq; r[2] += b
    Sd = A2 @ (n*nd)
    rd = -A2*nd/sq + A2*n*Sd/S**1.5
    return n, nd, r, rd, sq - b*n[2]

def Phi(ang):
    _, th, ps = ang
    return np.array([[np.cos(th)*np.cos(ps), np.sin(ps), 0.0],
                     [-np.cos(th)*np.sin(ps), np.cos(ps), 0.0],
                     [np.sin(th), 0.0, 1.0]])

def rolling_wd(w, ang):
    n, nd, r, rd, h = geom(ang, w)
    X = np.array([[0,-r[2],r[1]],[r[2],0,-r[0]],[-r[1],r[0],0]])
    W = np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    Ir = I - m*(X @ X)
    rhs = -W @ (I @ w) - m*(X @ (W @ (W @ r))) - m*(X @ (W @ rd)) + m*g*(X @ n)
    return np.linalg.solve(Ir, rhs)

def contact(w, ang, v, P):
    n, nd, r, rd, h = geom(ang, w)
    vp = v + np.cross(w, r); vp = vp - (vp @ n)*n
    d = np.cross(n, r); c = -(rd @ nd) + r @ np.cross(w, nd)
    A0 = Iinv @ (-np.cross(w, I @ w)); B = Iinv @ np.cross(r, n)
    if P['law'] == 'viscous':
        Fv = -P['cv']*vp
        A0 = Iinv @ (np.cross(r, Fv) - np.cross(w, I @ w)); A1 = np.zeros(3)
    else:
        u = vp/np.sqrt(vp @ vp + P['eps']**2)
        tau = np.cross(r, -P['mu']*u)
        if P['rho'] > 0:
            sc = P['rho']*(w @ n)
            sg = sc/np.sqrt(sc**2 + vp @ vp + P['eps']**2)
            tau = tau - P['mu']*P['rho']*sg*n
        A1 = Iinv @ tau
    if P.get('fixN'):                      # comparison only: N held at m*g
        N = m*g
    else:
        N = m*(g + c + d @ A0)/(1 - m*(d @ B) - m*(d @ A1))
    F = Fv if P['law'] == 'viscous' else -P['mu']*N*u
    wd = A0 + N*(B + A1)
    vd = (N*n + F)/m - g*n - np.cross(w, v)
    return wd, vd, N, F, vp

def rhs(t, y, P):
    w, ang = y[0:3], y[3:6]
    if P['law'] == 'rolling':
        wd = rolling_wd(w, ang)
        _, _, r, _, _ = geom(ang, w); v = -np.cross(w, r); vd = np.zeros(3)
    else:
        v = y[6:9]; wd, vd, _, _, _ = contact(w, ang, v, P)
    vs = rotmat(ang) @ v
    return np.concatenate([wd, np.linalg.solve(Phi(ang), w), vd, vs[:2]])

def ic(w3):
    ang = np.array([1e-6]*3); w = np.array([0, 0, w3])
    _, _, r, _, _ = geom(ang, w)
    return np.concatenate([w, ang, -np.cross(w, r), [0.0, 0.0]])

def simulate(P, w3, T=8.0, hz=1000):
    tev = np.linspace(0, T, int(T*hz)+1)
    s = solve_ivp(rhs, [0, T], ic(w3), t_eval=tev, args=(P,),
                  rtol=1e-8, atol=1e-10, method='LSODA')
    assert s.success, s.message
    return s.t, s.y.T
