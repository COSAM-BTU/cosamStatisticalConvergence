# Changelog

## 1.0.0

First public release, accompanying the article listed in `README.md`.

- `statisticalConvergence` function object for OpenFOAM v2406: marginal-standard-error
  (MSER-5) transient detection, autocorrelation integral time scale, effective sample
  size, small-sample confidence intervals for the mean (Student-*t*) and for the root
  mean square (chi-squared), and a Welch periodogram, all evaluated in situ.
- Signals may be single-cell probes or volume-weighted spatial averages of a field
  component or magnitude.
- Optional plateau criterion that writes and stops the run through the solver time
  control, with a per-signal gate, a settling time and a configurable comparison window.
- `statConvCheck`, a standalone utility sharing one compiled kernel with the function
  object, so that the verified code path is the one the solver executes.
- Tutorials reproducing the synthetic verification case of the article and showing the
  function object attached to a small transient case.
