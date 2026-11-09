# -*- coding: utf-8 -*-
"""
@author: kgbet
"""

from dataclasses import dataclass
from fractions import Fraction

def validate_scale(T, name="T"):
    """
    Validates one period of an isolated time scale
    and converts entries to exact fractions

    Parameters
    ----------
    T : time scale
    name : optional
        used to specify time scale in error message

    Returns
    -------
    T : validated time scale

    """
    if len(T) < 2:
        raise Exception("%s must include at least two pointt (endpoints of period)" % name)
    if any(T[i+1] <= T[i] for i in range(len(T)-1)):
        raise Exception("%s must be strictly increasing" % name)
    
    return T

def common_period(T1, T2, tol=1e-9):
    """
    Tests whether periods are commensurable and returns the common period

    Parameters
    ----------
    T1 : first time scale
    T2 : second time scale

    Returns
    -------
    commonperiod: the common period they share

    """
    T1 = validate_scale(T1,"T1")
    T2 = validate_scale(T2, "T2")
    period1 = T1[-1] - T1[0]
    period2 = T2[-1] - T2[0]
    r = Fraction(period1/period2).limit_denominator(1000)
    if abs(float(r) - period1/period2) > tol*max(1,period1/period2):
        raise Exception("Periods are not commensurable")
    
    commonperiod = period1*r.denominator
    return commonperiod

def extendtime(T, n):
    """Extends the first period T of a time scale n times (n+1 periods in total)."""
    i = 0
    timescale = T.copy()
    while i < n:
        for j in range(1, len(T)):
            timescale.append(T[j] + (T[-1] - T[0]) * (i + 1))
        i += 1
    return timescale

def extended_times(T1, T2, n_periods=1):
    
    T1 = validate_scale(T1,"T1")
    T2 = validate_scale(T2,"T2")
    if T1[0] != T2[0]:
        raise Exception("Time scales must start at the same time")
    P = common_period(T1,T2)
    k1 = int(P * n_periods / (T1[-1]-T1[0]))
    k2 = int(P * n_periods / (T2[-1]-T2[0]))
    
    extendedT1 = extendtime(T1,k1-1)
    extendedT2 = extendtime(T2,k2-1)
    return extendedT1,extendedT2


@dataclass
class Timing:
    events: list
    all_time: list
    end: float
    
def build_timing(T1,T2):
    T1 = validate_scale(T1,"T1")
    T2 = validate_scale(T2,"T2")
    if T1[0] != T2[0]:
        raise Exception("T1 and T2 must start together")
    map1 = {T1[i]:T1[i+1] for i in range(len(T1)-1)}
    map2 = {T2[i]:T2[i+1] for i in range(len(T2)-1)}
    set1 = set(T1)
    set2 = set(T2)
    events = []
    all_time = sorted(set1 | set2)
    for t in all_time:
        ax = (t in set1) and t != T1[0]
        ay = (t in set2) and t != T2[0]
        mu = float(map1[t] - t) if t in map1 else None
        nu = float(map2[t] - t) if t in map2 else None
        events.append((t, ax, ay, mu, nu))
    return Timing(events,all_time, all_time[-1])

def entire_period_timing(T1,T2):
    E1, E2 = extended_times(T1,T2)
    return build_timing(E1,E2)