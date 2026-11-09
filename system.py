# -*- coding: utf-8 -*-
"""
@author: kgbet
"""

import inspect
from dataclasses import dataclass
from timescale import common_period,extended_times,entire_period_timing, build_timing
import matplotlib.pyplot as plt

def filter_params(fun, params):
    "the entries of params that fun accepts"
    sig = inspect.signature(fun)
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
        return dict(params)
    return {k: v for k, v in params.items() if k in sig.parameters}

def fun_wrapper(fun, params):
    kw = filter_params(fun, params)
    return lambda x, y: fun(x, y, **kw)


def partials(f, x, y, jac=None):
    "returns partials in x and y directions of function f at specified x and y"
    if jac is not None:
        fx, fy = jac(x,y)
        return fx,fy
    hx = 1e-6* max(1,abs(x))
    hy = 1e-6 * max(1,abs(y))
    fx = (f(x+hx,y) - f(x-hx,y))/(2*hx)
    fy = (f(x,y+hy) - f(x,y-hy))/(2*hy)
    return fx,fy

def run_timing(events, F, G, x0, y0):
    "runs all events and keeps all states of x and y"
    lastx = x0
    lasty = y0
    x = [x0]
    y= [y0]
    for t, ax, ay, mu, nu in events:
        if ax:
            lastx = x[-1]
        if ay:
            lasty = y[-1]
        if mu is not None:
            x.append(lastx + F(lastx,lasty)*mu)
        if nu is not None:
            y.append(lasty + G(lastx,lasty)*nu)
    return x, y

def linearF(x, y, alpha=0, beta=0):
    return alpha*x + beta*y

def linearG(x, y, gamma=0, delta=0):
    return gamma*x + delta*y

def plot(x, E, xname=None, Ename='Time'):
    fig, ax = plt.subplots()
    ax.plot(E,x, '--', marker='o')
    ax.set_xlabel(Ename)
    ax.set_ylabel(xname or '')
    ax.set_title('Simulated Trajectory')
    plt.show()
    return fig, ax

def plot2D(x, y, E1, E2):
    E = sorted(set(E1) | set(E2))
    xnew = [x[0]]
    ynew = [y[0]]
    for t in E[1:]:
        if t in E1:
            i = E1.index(t)
            xnew.append(x[i])
        else:
            xnew.append(xnew[-1])
        if t in E2:
            i = E2.index(t)
            ynew.append(y[i])
        else:
            ynew.append(ynew[-1])
    fig, ax = plt.subplots()
    ax.plot(xnew,ynew, '--', marker='o')
    plt.show()
    return fig, ax

@dataclass

class AsynchSystem:
    
    def __init__(self, T1, T2, F, G, params=None, name="", jac_F=None, jac_G=None):
        self.params = params or {}
        self.T1_0, self.T2_0 = T1, T2
        self.T1 = T1(**filter_params(T1, self.params)) if callable(T1) else T1
        self.T2 = T2(**filter_params(T2, self.params)) if callable(T2) else T2
        self.F0 = F
        self.G0 = G
        self.F = fun_wrapper(F, self.params)
        self.G = fun_wrapper(G, self.params)
        self.name = name
        self.timing = entire_period_timing(self.T1,self.T2)
        self.period = common_period(self.T1,self.T2)
    
    def simulate(self,x0,y0, n_periods=1):
        E1, E2 = extended_times(self.T1, self.T2, n_periods=n_periods)
        timings = build_timing(E1,E2)
        events = timings.events
        x, y = run_timing(events, self.F, self.G, x0, y0)
        return x, y, E1, E2
    
    @classmethod
    def linear(cls, T1, T2, alpha, beta, gamma, delta):
        "creates a linear system asynchsystem object"
        params = dict(alpha=alpha, beta = beta, gamma=gamma,delta=delta)
        return cls(T1, T2, linearF, linearG, params=params, name='linear system')
    def uniform_time(cls, mu, nu, F, G, **kw):
        "creates a system with underlying time 0, mu, 2mu... and 0, nu, 2nu..."
        return cls([0,mu], [0,nu], F, G, **kw)
    def replace(self, **kw):
        "Copy a system with specified change in parameter or time scale"
        T1 = kw.pop("T1", self.T1_0)
        T2 = kw.pop("T2", self.T2_0)
        pars = dict(self.params)
        pars.update(kw)
        return AsynchSystem(T1, T2, self.F0, self.G0, params=pars, name=self.name)
    def __repr__(self):
        return "AsyncSystem(%s, T1=%s, T2=%s, period=%s, params=%s)" % (self.name, self.T1, self.T2, self.period, self.params)
     
    def periodic_map(self, z, numsteps=1):
        "steps the system forward numsteps common periods"
        try:
            E1, E2 = extended_times(self.T1, self.T2, n_periods=numsteps)
            timings = build_timing(E1,E2)
            events = timings.events
            x, y = run_timing(events, self.F, self.G, z[0], z[1])
            return x[-1], y[-1]
        except:
            raise Exception("Aw shucks! Something's wrong here (probably Nans)")
    def FG(self, x, y):
        return self.F(x,y), self.G(x,y)