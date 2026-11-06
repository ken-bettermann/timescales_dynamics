# Asynchronous Discrete Dynamical Systems

Python tools for simulating and analyzing the stability of **asynchronously updated discrete dynamical systems**: two coupled variables, `x` and `y`, each updated on its own discrete time scale, where each update uses the most recently available value of the other variable.
Python tools for simulating and analyzing the stability of 2D asynchronous discrete dynamical systems. Each variable updates to the next point on its time based on its current value and the value of the most recently available value of the other variable.

This code is based on the linear system posed by Dr. Stefan Siegmund and Dr. Petr Stehlik in `Time-scale Induced Asynchronous Discrete Dynamical Systems` and is a result of research supervised by Dr. Sabrina Streipert and funded by Pitt Momentum.

---

## Background
Many models contain components that do not update at the same moment: a continuously reproducing host and a seasonally reproducing vector, processors in a distributed algorithm, or components in a control loop with delayed data transfer. This project studies systems of the form

```
x<sup>Δ</sup>(t) = F( x(t), y(t*) )      t ∈ T1
y<sup>Δ</sup>(t) = G( x(t*), y(t) )      t ∈ T2
```

where `T1` and `T2` are periodic isolated time scales and `t*` denotes the most recent update time of the *other* variable (sample-and-hold). Concretely, at each time in `T1`, `x` takes a forward-Euler step of length `μ(t)` (the graininess function in `T1`) using the last known value of `y`; `y` is handled the same way on `T2` using graininess function `ν(t)`.

### Stability method

Following Siegmund, the delayed system is rewritten as a delay-free system in four variables `(u, ū, v, v̄)`, where `ū` and `v̄` store the "last seen" values of `x` and `y`. At each event time we multiply by the 4x4 linearization of dynamics. Stability can only be computed for time scales `T1` and `T2` which have commensurable period, i.e. period(T1)/period(T2) = q∈ℚ. In this circumstance, there exists some `T` such that the combined time scale `T1 ∪ T2` is T-periodic. To find stability under this condition:

1. Multiply the 4×4 matrices to get the evolution operator `Φ(T, 0)`.
2. Take the 2×2 submatrix on the `u`/`v` coordinates (rows/columns 1 and 3) to get the solution operator `Ψ(T, 0)`.
3. Its spectral radius `ρ(Ψ(T, 0))` decides stability.

The stability condition at fixed point `(x*,y*)` is determined as follows:

- `ρ < 1` → locally asymptotically stable
- `ρ > 1` → unstable
- `ρ = 1` → inconclusive

For linear `F`, `G` and commensurable time scales this recovers Siegmund's result at the origin.

---

## Repository contents

| File | Purpose |
|---|---|
| `timescale.py` | Time-scale utilities: validation, common period, periodic extension, and event/timing construction |
| `system.py` | The `AsynchSystem` class: parameters, simulation, periodic map, and plotting |
| `stability.py` | Fixed-point search, Siegmund matrix construction, local stability, and 1-parameter bifurcation diagrams |
| `Pitt_Momentum_S2025.pdf` | Research notes: theory, proofs, and examples behind the code |

---

## Installation

Requires Python 3.9+ and:

```
numpy
scipy      # >= 1.15
matplotlib
```

```bash
pip install numpy "scipy>=1.15" matplotlib
```

No package install is needed. Keep the three `.py` files in the same folder.

---

## Quick start

### 1. Defineing a system

A time scale is given as **one period** of strictly increasing points, and both scales must start at the same time. `F` and `G` take `(x, y)` plus optional keyword parameters, which are picked up by name from `params`.

```python
from system import AsynchSystem

F = lambda x, y, r: r - x - x**3 - 3*y
G = lambda x, y: 0.5 * (x - y)

# x updates at t = 0, 1, 2, ...   y updates at t = 0, 1, 2, ...
S = AsynchSystem([0, 1], [0, 1], F, G, params={'r': 1.0}, name='example')

print(S.period)   # common period of T1 and T2
```

Time scales may also depend on parameters by passing a function, allowing for adjustable underlying time.

```python
S = AsynchSystem(lambda mu: [0, mu, 3*mu], lambda nu: [0, nu], F, G,
                 params={'mu': 1, 'nu': 2, 'r': 1.0})
```

### 2. Simulating

```python
from system import plot, plot2D

x, y, E1, E2 = S.simulate(x0=0.5, y0=0.5, n_periods=3)
plot(x, E1, xname='x')
plot(y, E2, xname='y')
plot2D(x, y, E1, E2)       # phase portrait (sample-and-hold)
```

`S.periodic_map([x0, y0], numsteps=1)` returns the state after one (or more if numsteps!=1) common periods.

### 3. Find fixed points and test stability

```python
from stability import fixed_points, stability

fps = fixed_points(S, bounds=((-5, 5), (-5, 5)), n=10)
for fp in fps:
    verdict, rho = stability(S, fixed_point=fp)
    print(fp, verdict, rho)
```

### 4. Bifurcation diagram in one parameter

```python
import numpy as np
from stability import bifurcation_1param, plot_bif1param

recs = bifurcation_1param(S, 'r', np.linspace(0.5, 4.5, 200),
                          bounds=((-5, 5), (-5, 5)), n=5)
plot_bif1param(recs, 'r', coordinate=1)   # red = stable, blue = unstable
```

### 5. Linear systems

```python
S = AsynchSystem.linear([0, 2], [0, 3], alpha=-0.5, beta=0.3, gamma=0.2, delta=-0.3)
```

For linear systems the origin is the fixed point, and `stability(S)` reproduces Siegmund's criterion.

### Commensurable time scales

The periods of `T1` and `T2` must be commensurable (their ratio rational). When sweeping a step size, use `fractions.Fraction` to ensure successful sweeping:

```python
from fractions import Fraction
nus = [Fraction(k, 50) for k in range(4, 300)]
recs = bifurcation_1param(S, 'nu', nus)
```

---

## API summary

**`timescale.py`**
- `validate_scale(T)` checks that `T` has at least two points and is strictly increasing.
- `common_period(T1, T2)` returns the common period, or raises if the periods are not commensurable.
- `extendtime(T, n)` repeats one period `n` more times.
- `extended_times(T1, T2, n_periods=1)` extends both scales over `n_periods` common periods.
- `build_timing(T1, T2)` builds the ordered event list `(t, update_x, update_y, mu, nu)`.
- `entire_period_timing(T1, T2)` builds the timing for one common period.

**`system.py`**
- `AsynchSystem(T1, T2, F, G, params=None, name="")` is the main class. Methods: `simulate`, `periodic_map`, `FG`, `replace` (copy with changed parameters or time scales), and `linear(...)` (classmethod).
- `run_timing`, `plot`, `plot2D`, `partials`, `fun_wrapper`, `filter_params` are helpers.

**`stability.py`**
- `fixed_points(system, ...)` finds solutions of `F = G = 0`.
- `siegmund_matrix(...)` builds the 4×4 linearization for one step.
- `psiT(T1, T2, fx, fy, gx, gy)` composes the matrices over one period and returns the 2×2 evolution operator `Ψ(T, 0)`.
- `stability(system, fixed_point, tol)` returns `(verdict, spectral_radius)`.
- `bifurcation_1param(...)` and `plot_bif1param(...)` compute and plot stability against a parameter.

---

## Known limitations

- Only 2-variable systems are supported.
- Only **local** stability is tested for nonlinear systems; `ρ = 1` is reported as inconclusive.
- Time scales must be commensurable and share the same starting point.
- Fixed-point search is grid-based, so some fixed points may be missed if the bounds or grid are too coarse.

---

## Reference

- Siegmund, Stefan, and Petr Stehlik. "Time scale-induced asynchronous discrete dynamical systems." Discrete & Continuous Dynamical Systems-Series B 26.2 (2021): 1011.

## How to cite

If you use this software, please cite it as follows:
```bibtex
@software{timescales_dynamics2026,
  author        = {Kenneth Bettermann},
  title         = {timescales_dynamics},
  month         = November,
  year          = 2026,
  publisher     = {GitHub},
  journal       = {GitHub repository},
  howpublished  = {\url{[https://github.com/ken-bettermann/timescales_dynamics](https://github.com/ken-bettermann/timescales_dynamics)}}
}
