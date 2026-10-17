"""Run every simulation used in the paper and save it to runs/<name>.npz.

Usage:  python3 run_simulations.py            (all runs, two worker processes)
        python3 run_simulations.py --list     (show the run names)
About 25 runs; roughly 10-15 minutes on two cores.
"""
import os, sys, time
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from model import simulate

OUT = os.path.join(HERE, 'runs'); os.makedirs(OUT, exist_ok=True)
EPS = 3e-5                                    # converged regularisation (paper, Fig. 2)
C = lambda mu=0.15, rho=0.0, eps=EPS, **k: dict(law='coulomb', mu=mu, eps=eps, rho=rho, **k)
V = lambda cv: dict(law='viscous', cv=cv)
R = dict(law='rolling')

JOBS = {                                      # name: (parameters, initial spin, output rate in Hz)
    'rolling_nonpref': (R, -5, 1000), 'rolling_pref': (R, +5, 1000),
    'coulomb_non': (C(), -5, 1000), 'coulomb_pref': (C(), +5, 500),
    'spin05_non': (C(rho=0.5e-3), -5, 1000), 'viscous_nonpref': (V(500.0), -5, 1000),
}
for e, tag in [(1e-2,'0.01'), (3e-3,'0.003'), (1e-3,'0.001'), (3e-4,'0.0003'), (1e-4,'0.0001'), (1e-5,'1e-05')]:
    JOBS[f'eps_{tag}'] = (C(eps=e), -5, 500)
for mu in (0.05, 0.2, 0.5, 1.0, 2.0): JOBS[f'mu_{mu}'] = (C(mu=mu), -5, 500)
for mu in (0.15, 0.2, 0.5, 1.0, 2.0): JOBS[f'fixN_mu_{mu}'] = (C(mu=mu, fixN=True), -5, 500)
for r in (0.1, 0.2, 0.5, 1.0, 2.0):
    if r != 0.5: JOBS[f'rho_{r}_non'] = (C(rho=r*1e-3), -5, 500)
    JOBS[f'rho_{r}_pref'] = (C(rho=r*1e-3), +5, 500)
for c in (1.5, 10.0, 50.0, 200.0, 2000.0): JOBS[f'cv_{c:g}'] = (V(c), -5, 500)

def run(item):
    name, (P, w3, hz) = item
    path = os.path.join(OUT, name + '.npz')
    if os.path.exists(path): return name, 'already done'
    t0 = time.time(); t, Y = simulate(P, w3, T=8.0, hz=hz)
    np.savez_compressed(path, t=t, Y=Y); return name, f'{time.time()-t0:.0f} s'

if __name__ == '__main__':
    if '--list' in sys.argv:
        for k in JOBS: print(k)
        sys.exit()
    with Pool(2) as p:
        for name, msg in p.imap_unordered(run, JOBS.items()): print(f'{name:22s} {msg}', flush=True)
