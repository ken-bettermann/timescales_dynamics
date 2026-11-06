# -*- coding: utf-8 -*-
"""
@author: kgbet
"""

from dataclasses import dataclass
import numpy as np
from scipy.differentiate import jacobian
from scipy.optimize import root
from timescale import extended_times
from system import AsynchSystem
import matplotlib.pyplot as plt

@dataclass
class StabilityResult:
    fixed_point:tuple
    period:float
    jacobian: np.ndarray
    eigenvalues: np.ndarray
    spectral_radius: float
    result: str
    

def determine(rho, tol=1e-9):
    if rho<1-tol:
        return 'asymptotically stable'
    elif rho>1+tol:
        return 'unstable'
    return 'inconclusive'

def fixed_points(system, guesses=None, bounds=((-10, 10), (-10, 10)), n=10, tol=1e-9):
    "Finds fixed points of the system where F(x,y)=G(x,y)=0 by searching with Newton's method"
    
    if guesses is None:
        guesses = [(a,b) for a in np.linspace(*bounds[0], n) for b in np.linspace(*bounds[1], n)]
    results = []
    FG = lambda z: system.FG(z[0],z[1])
    for guess in guesses:
        sol = root(FG, guess, method='hybr')
        if not sol.success:
            continue
        if max(abs(np.array(FG(sol.x)))) > tol:
            continue
        if not any(np.allclose(sol.x, p, atol=1e-6) for p in results):
            results.append(sol.x)
    return [tuple(p) for p in results]

def siegmund_matrix(t1, t2, s1, s2, mu, nu, fx, fy, gx, gy):
    "gives A_{ijkl}"
    both = t1 and t2
    A = np.zeros((4,4))
    # u
    if not s1:
        A[0,0] = 1
    else:
        A[0,0] = 1 + mu * fx
        A[0, 2 if both else 3] = mu * fy
    # ubar
    if both or ((not t2) and s2):
        A[1,0] = 1
    else:
        A[1,1] = 1
    # v
    if not s2:
        A[2,2] = 1
    else:
        A[2,2] = 1 + nu * gy
        A[2, 0 if both else 1] = nu * gx
    #vbar
    if both or ((not t1) and s1):
        A[3,2] = 1
    else:
        A[3,3] = 1
    return A

def psiT(T1, T2, fx, fy, gx, gy):
    "gives the 2x2 submatrix we care about"
    E1, E2 = extended_times(T1, T2)
    set1,set2 = set(E1), set(E2)
    times = sorted(set1 | set2)
    M = np.eye(4)
    for i, t in enumerate(times[:-1]):
        s = times[i+1]
        mu = 0
        nu = 0
        if s in set1:
            mu = s - max(p for p in E1 if p <= t)
        if s in set2:
            nu = s - max(p for p in E2 if p <= t)
        A = siegmund_matrix(t in set1, t in set2, s in set1, s in set2, mu, nu, fx, fy, gx, gy)
        M = A @ M
    Psi = np.array([[M[0,0], M[0,2]], [M[2,0], M[2,2]]])
    return Psi
        

def stability(system, fixed_point=(0.0,0.0), tol=1e-9):
    "local stability of fixed point of the asynchronous system"
    
    x0, y0 = fixed_point[0], fixed_point[1]
    check1, check2 = system.FG(x0,y0)
    if abs(check1)>tol or abs(check2)>tol:
        raise Exception("That's not a fixed point!")
    f = lambda z: system.F(z[0],z[1])
    g = lambda z: system.G(z[0],z[1])
    res1 = jacobian(f,[x0,y0])
    df = res1.df
    res2 = jacobian(g,[x0,y0])
    dg = res2.df
    Psi = psiT(system.T1, system.T2, df[0], df[1], dg[0], dg[1])
    eigs = np.linalg.eigvals(Psi)
    rho = max(abs(eigs))
    return determine(rho), rho

def bifurcation_1param(system, par_name, par_vals, bounds=((-10, 10), (-10, 10)), n=5):
    "records (parameter, fixed point, spectral radius) for every fixed point at each parameter value"
    records = []
    for par in par_vals:
        s = system.replace(**{par_name: par})
        for fp in fixed_points(s, bounds=bounds, n=n):
            res, rho = stability(s, fixed_point=fp, tol=1e-6)
            records.append((par, fp, rho))
    return records

def plot_bif1param(records, par_name, coordinate):
    "Plots the result of the bifurcation_1param function with x-axis parameter and y-axis specified coordinate of fixed point"
    pars = np.array([r[0] for r in records])
    xs   = np.array([r[1][coordinate-1] for r in records])
    rhos = np.array([r[2] for r in records])

    stable = rhos < 1

    plt.scatter(pars[stable],  xs[stable],  c='red',  s=15, label='stable')
    plt.scatter(pars[~stable], xs[~stable], c='blue', s=15, label='unstable')
    plt.xlabel(par_name)
    plt.ylabel("coordinate %s of fixed point" % coordinate)
    plt.legend()
    plt.show()
    return
    
F = lambda x, y, r: r - x - x**3 - 3*y
G = lambda x, y: 0.5 * (x - y)
S = AsynchSystem([0, 1], [0, 1], F, G, params={'r': 1.0})
recs = bifurcation_1param(S, 'r', np.linspace(0.5, 4.5, 500), bounds=((-5, 5), (-5, 5)), n=5)
plot_bif1param(recs, 'r', 1)

F = lambda x, y, r: r - x**2
G = lambda x, y: x - y

S = AsynchSystem([0, 1], [0, 1], F, G, params={'r': 1.0})
recs = bifurcation_1param(S, 'r', np.linspace(-0.5, 1.5, 500), bounds=((-10, 10), (-10, 10)), n=5)
plot_bif1param(recs, 'r', 2)

from fractions import Fraction

linearF = lambda x, y, alpha=0, beta=0: alpha*x + beta*y
linearG = lambda x, y, gamma=0, delta=0: gamma*x + delta*y

S = AsynchSystem(lambda mu: [0, mu], lambda nu: [0, nu], linearF, linearG,
                 params=dict(mu=1, nu=1, alpha=-0.5, beta=0.3, gamma=0.2, delta=-0.3))

nus = [Fraction(k, 50) for k in range(4, 300)] # note that we need commensurability so this is weird.
recs = bifurcation_1param(S, 'nu', nus)
plot_bif1param(recs, 'nu', 1)