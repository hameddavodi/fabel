# Registration: aligning the growth spurt

Every girl has a *pubertal growth spurt*: a few years in which she grows fast,
then slows down. But the spurt comes at different ages, from about 9 to 14.
If you average the curves as they are, the spurts smear out and the mean curve
shows a weak, wide spurt that no girl really had.

Curves can differ in two ways:

- **amplitude**: how strong a feature is (a big spurt or a small one);
- **phase**: when the feature happens (an early spurt or a late one).

*Registration* removes the phase differences. It bends the time axis of each
curve with a *warping function* h(t), so that the features line up. A warping
function always increases (time never runs backwards) and keeps the start and
end of the interval fixed.

This tutorial uses the growth data, which ships with fdatools.

## 1. Growth acceleration curves

We smooth 20 girls with an order-6 B-spline and a penalty on the fourth
derivative, so that the second derivative (the acceleration) is smooth too.

```python
import matplotlib.pyplot as plt
import numpy as np
import fdatools as fdt
from fdatools.registration import landmark_register, register

growth = fdt.datasets.load_growth()
age = growth.age
basis = fdt.BSpline(domain=(1.0, 18.0), n_basis=35, order=6)
fit = fdt.smooth(growth.hgtf[:, :20], age, basis=basis, penalty=4, lam=1e-2)

velocity = fit.fd.derivative(1)       # cm / year
accel = fit.fd.derivative(2)          # cm / year^2, exact: a lower-order spline
t = np.linspace(3.0, 17.0, 300)

fig, ax = plt.subplots()
ax.plot(t, accel(t), color="grey", alpha=0.6)
ax.plot(t, accel.mean()(t), color="black", linewidth=2, label="mean")
ax.axhline(0.0, linestyle=":")
ax.legend()
ax.set(xlabel="age (years)", ylabel="acceleration", title="Before registration")
```

We plot from age 3: near age 1 the acceleration of an infant is very large
and hides everything else.

## 2. Landmark registration

A *landmark* is a feature you can find in every curve. Here we use the age of
peak growth speed: the maximum of the velocity between age 8 and 16. We
search only inside that interval, because a spline's derivatives can be large
and unreliable right at the ends of the data.

```python
search = np.linspace(8.0, 16.0, 801)
peak_age = search[np.argmax(velocity(search), axis=0)]   # one age per girl
peak_age.round(1)
```

`register(fd, landmarks=...)` warps every curve so that its landmark moves to
the mean landmark age (or to `target_landmarks=` if you give them).
`landmark_register(fd, landmarks)` is the same call with a shorter name.

```python
landmark = register(accel, landmarks=peak_age)

target = peak_age.mean()
landmark.warp_values(np.array([target])).round(2)   # h_i(target) = each girl's own peak age

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(t, landmark.registered(t), color="grey", alpha=0.6)
axes[0].plot(t, landmark.registered.mean()(t), color="black", linewidth=2)
axes[0].axvline(target, linestyle=":")
axes[0].set(xlabel="registered age", title="Landmark registered")
axes[1].plot(t, landmark.warp_values(t))
axes[1].plot(t, t, color="black", linestyle="--")
axes[1].set(xlabel="registered age", ylabel="clock age h(t)", title="Warping functions")
fig.tight_layout()
```

Read a warping function like this: at registered age `t`, girl `i` was really
`h_i(t)` years old. At the mean peak age, a girl whose warp lies above the
diagonal had her spurt later than average; below it, earlier.

The result is a `RegistrationResult`:

- `registered`: the aligned curves x_i(h_i(t));
- `warp`: the warping functions as curves; `warp_values(t)` evaluates them
  exactly;
- `warp_inverse`: the inverse warps (landmark registration only);
- `unregistered`: the input curves.

## 3. Continuous registration

Landmarks need you to find a feature by hand. *Continuous registration* does
not: it warps each curve to be as close as possible to a target curve (the
mean, by default). fdatools' default criterion is the *minimum eigenvalue*
criterion of Ramsay and Silverman: it asks each registered curve to be
proportional to the target, so amplitude differences do not disturb the warp.

The warping functions are built from a latent curve W in a B-spline basis
(`warp_basis`). More basis functions allow more complex warps; `lam`
penalises rough warps.

```python
warp_basis = fdt.BSpline(domain=(1.0, 18.0), n_basis=6)
continuous = register(accel, warp_basis=warp_basis, lam=1.0)

continuous.n_iter                 # Newton iterations per curve
continuous.criterion.round(3)     # final criterion value per curve

fig, ax = plt.subplots()
ax.plot(t, continuous.registered(t), color="grey", alpha=0.6)
ax.plot(t, continuous.registered.mean()(t), color="black", linewidth=2)
ax.set(xlabel="registered age", title="Continuous registration")
```

For periodic data such as the weather, pass `periodic=True`: then a time
shift per curve is estimated as well.

## 4. How much variation was phase?

Kneip and Ramsay (2008) split the total variation of the curves into an
amplitude part and a phase part. `.decompose()` returns the two mean squared
errors, their ratio `rsq` (the share due to phase) and a constant `c`. We
leave out the first two years, where the infant acceleration dominates.

```python
amp_mse, phase_mse, rsq, c = landmark.decompose(domain=(3.0, 17.0))
round(rsq, 2)                     # about 0.6: most variation near the spurt is timing

round(continuous.decompose(domain=(3.0, 17.0)).rsq, 2)   # about 0.55
```

About 60% of the variation in growth acceleration between ages 3 and 17 is
phase: the girls differ mostly in *when* they have their spurt, and less in
how strong it is.

`phase_mse` can be negative. That happens when a registration makes the
curves *less* aligned than before, and it is a useful warning sign.

## 5. In a scikit-learn pipeline

`Registrator` wraps continuous `register()` as a scikit-learn transformer.
`fit` registers the training curves to their mean and keeps that mean as the
target; `transform` registers curves to the same target. It accepts an
`FData` or, inside a pipeline, a coefficient matrix with one row per curve
plus the curves' basis:

```python
from fdatools.registration import Registrator

reg = Registrator(warp_basis, lam=1.0, basis=accel.basis)
registered_coefs = reg.fit_transform(np.asarray(accel.coefs).T)
registered_coefs.shape            # (20, 33): 20 girls, 33 basis functions
```

## R equivalent

```r
library(fda)
hgtbasis <- create.bspline.basis(c(1, 18), 35, norder = 6)
hgtfd    <- smooth.basis(growth$age, growth$hgtf[, 1:20],
                         fdPar(hgtbasis, 4, 1e-2))$fd
accelfd  <- deriv.fd(hgtfd, 2)
agefine  <- seq(8, 16, length.out = 801)
pgs      <- agefine[apply(eval.fd(agefine, hgtfd, 1), 2, which.max)]

# landmark registration
wbasis   <- create.bspline.basis(c(1, 18), norder = 4,
                                 breaks = c(1, mean(pgs), 18))
lmreg    <- landmarkreg(accelfd, ximarks = pgs, x0marks = mean(pgs),
                        WfdPar = fdPar(fd(matrix(0, 5, 1), wbasis), 2, 1e-4))
lmreg$regfd; lmreg$warpfd

# continuous registration
Wfd0     <- fdPar(fd(matrix(0, 6, 20), create.bspline.basis(c(1, 18), 6)), 2, 1)
contreg  <- register.fd(mean.fd(accelfd), accelfd, Wfd0)

AmpPhaseDecomp(accelfd, contreg$regfd, contreg$warpfd, c(3, 17))
```

Differences to know:

- One `register()` covers `register.fd` and `landmarkreg`; the argument
  `landmarks=` switches between them.
- `deriv.fd` re-expands the derivative in the original basis and is about 1%
  off; `fd.derivative()` is exact.
- fdatools' continuous criterion uses the same discretisation as R, so `lam`
  means the same thing. Its Newton method uses the exact Hessian and iterates
  until the gradient is zero; R often stops earlier, so the warps can differ.
- The default warp basis for landmarks is a cubic spline with knots at the
  target landmarks and `lam=1e-4`; for continuous registration it is R's
  default, a two-function linear basis with `lam=0`.
