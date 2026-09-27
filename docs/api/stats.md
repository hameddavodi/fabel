# Statistics

`fabel.stats`: covariance and correlation of functional data, functional
depth, the functional boxplot, and permutation t and F tests.

`f_test` takes either raw regression inputs, `f_test(y, x, basis=, lam=,
penalty=)` like R's `Fperm.fd`, or a model fitted by
[`fregress`](regression.md), `f_test(model, n_perm=1000)`, which refits that
exact model (its terms, intercept included, coefficient settings and weights)
under every permutation.

::: fabel.stats
