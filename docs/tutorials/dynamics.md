# Dynamics: differential equations from data

Curves often come from a system that obeys a *differential equation*: a rule
that links a curve to its own derivatives. A spring, for example, follows

D²x(t) + β₁ Dx(t) + β₀ x(t) = 0,

where `Dx` is the velocity and `D²x` the acceleration. β₀ sets how stiff the
spring is and β₁ how quickly the motion dies out (the *damping*).

*Principal differential analysis* (PDA) turns this around. Given many curves,
it estimates the weights β₀(t), β₁(t) that make the left-hand side as close to
zero as possible for all of them. The weights may change over time.

This tutorial studies the movement of the lower lip while a speaker says
"bob" 20 times (downloaded once, then cached), and then looks at the growth
data in the *phase plane*.

## 1. The lip data

```python
import matplotlib.pyplot as plt
import numpy as np
import fdatools as fdt
from fdatools.dynamics import PDA, phase_plane
```

We smooth the 20 lip curves with an order-6 B-spline and a small penalty on
the fourth derivative, so that the first two derivatives are smooth.

```python
# requires-data: lip
lip = fdt.datasets.load_lip()
t = lip.t                                      # 51 time points, 0 to 0.35 s
lip_basis = fdt.BSpline(domain=(0.0, 0.35), n_basis=31, order=6)
lip_fd = fdt.smooth(lip.value, t, basis=lip_basis, penalty=4, lam=1e-8).fd

fig, ax = plt.subplots()
lip_fd.plot(ax=ax, color="grey", alpha=0.6)
ax.set(xlabel="time (s)", ylabel="lip position (mm)", title="20 repetitions of 'bob'")
```

## 2. A constant-coefficient equation

The simplest PDA uses constant weights. `PDA` is a scikit-learn estimator:

```python
# requires-data: lip
constant = PDA(order=2).fit(lip_fd)
beta0 = float(constant.weights_[0].coefs[0, 0])
beta1 = float(constant.weights_[1].coefs[0, 0])
round(beta0, 2), round(beta1, 2)               # about 54.4 and 2.1
```

For a constant-coefficient spring the period of the oscillation is about
2π / sqrt(β₀):

```python
# requires-data: lip
round(2 * np.pi / np.sqrt(beta0), 2)           # about 0.85 seconds
```

`weights_[j]` is β_j as a one-curve `FData`. `operator_` is the fitted
operator L = β₀ + β₁ D + D² as an `LDO`, and `residuals_` holds `Lx` for every
curve: the part of each curve the equation does not explain.

## 3. Weights that change over time

Speech is not a simple spring: the lip opens, closes and opens again. Let the
weights vary with a B-spline basis, with a small roughness penalty:

```python
# requires-data: lip
weight_basis = fdt.BSpline(domain=(0.0, 0.35), n_basis=21)
pda = PDA(order=2, weight_basis=weight_basis, lam=1e-8).fit(lip_fd)

grid = np.linspace(0.0, 0.35, 200)
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(grid, pda.weights_[0](grid))
axes[0].set(xlabel="time (s)", title="beta 0 (stiffness)")
axes[1].plot(grid, pda.weights_[1](grid))
axes[1].axhline(0.0, color="grey", linewidth=0.5)
axes[1].set(xlabel="time (s)", title="beta 1 (damping)")
fig.tight_layout()
```

### Does the equation reproduce the data?

`solve()` integrates the fitted equation Lx = 0 from a starting position and
velocity. Start it where the mean lip curve starts and compare:

```python
# requires-data: lip
mean_lip = lip_fd.mean()
start = np.array([0.0])
initial = [float(mean_lip(start)[0, 0]), float(mean_lip(start, 1)[0, 0])]

path_varying = pda.solve(t, initial)
path_constant = constant.solve(t, initial)

fig, ax = plt.subplots()
ax.plot(t, mean_lip(t), color="black", linewidth=2, label="mean curve")
ax.plot(t, path_varying, "--", label="time-varying weights")
ax.plot(t, path_constant, ":", label="constant weights")
ax.legend()
ax.set(xlabel="time (s)", ylabel="lip position (mm)")

gap = float(np.max(np.abs(path_varying - mean_lip(t)[:, 0])))
round(gap, 2)                                  # under 1 mm on a 22 mm movement
```

The time-varying equation follows the mean curve closely. The constant one
cannot: one fixed stiffness is not enough to describe opening and closing.

### The stability diagram

`plot_overlay()` draws the path (β₁(t), β₀(t)) over the parabola
β₀ = β₁² / 4. Above the parabola the system *oscillates*; below it, it only
moves towards rest. Left of zero (β₁ < 0) the motion grows instead of dying
out.

```python
# requires-data: lip
ax = pda.plot_overlay(labels={0.0: "start", 0.175: "middle", 0.35: "end"})
ax.set_title("Lip dynamics on the stability diagram")
```

## 4. The phase plane

A *phase-plane plot* draws one derivative against another. The default of
`phase_plane` is acceleration (`D²x`) against velocity (`Dx`). An oscillation
then traces a loop around the origin. Energy is a useful way to read it: the
horizontal distance from the origin shows kinetic energy (speed), the
vertical one potential energy (force).

Here are the height curves of three girls from the growth data (it ships with
fdatools). The pubertal growth spurt is a large loop: speed rises, peaks when the
acceleration crosses zero, and falls again.

```python
growth = fdt.datasets.load_growth()
girls = fdt.smooth(
    growth.hgtf[:, :3], growth.age,
    basis=fdt.BSpline(domain=(1.0, 18.0), n_basis=35, order=6), penalty=4, lam=1e-2,
).fd

ages = np.linspace(4.0, 17.0, 300)
ax = phase_plane(girls, ages, labels={6.0: "6", 11.0: "11", 14.0: "14"})
ax.axhline(0.0, color="grey", linewidth=0.5)
ax.set(xlabel="velocity (cm / year)", ylabel="acceleration (cm / year²)")
```

`deriv=(0, 1)` plots position against velocity instead.

## R equivalent

```r
library(fda)
lipbasis <- create.bspline.basis(range(liptime), nbasis = 31, norder = 6)
lipfd    <- smooth.basis(liptime, lip, fdPar(lipbasis, 4, 1e-8))$fd

# constant weights
cbasis   <- create.constant.basis(range(liptime))
pdaconst <- pda.fd(list(lipfd), list(fdPar(cbasis), fdPar(cbasis)))

# time-varying weights
wbasis   <- create.bspline.basis(range(liptime), 21)
pdavary  <- pda.fd(list(lipfd), list(fdPar(wbasis, 2, 1e-8), fdPar(wbasis, 2, 1e-8)))
pda.overlay(pdavary)

phaseplanePlot(seq(4, 17, length.out = 300), hgtfd[1:3])
```

Differences to know:

- R's `pda.fd` integrates with the trapezoid rule on 501 points. fdatools does the
  same by default (`n_grid=501`), so the weights match R. `n_grid=None` uses
  exact integrals instead; the results differ from R by about 1e-5.
- R has no `solve()`; you would call an ODE solver yourself.
- Weight functions of a system (several variables) are nested like R's
  `bwtlist`: `weights_[i][k][j]`. Forcing functions (`awtlist`, `ufdlist`)
  and fixed weights are not available.
