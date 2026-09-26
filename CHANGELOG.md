# Changelog

## 1.1.0

The statistics computed by the function object and `statConvCheck` are unchanged; results
computed with 1.0.0 are reproduced exactly.

- Function object: the sample times are now reset together with the signal histories when the
  dictionary is re-read during a run (previously only the histories were cleared), and the stop
  message states the tolerance relative to `max(|mean|, refScale)`, which is what is tested.
- New `tools/coverage.py`: coverage of the intervals of the kernel on synthetic stationary
  AR(1) records, with and without MSER truncation, together with the median effective sample
  size that the tool reports and the frequency of large MSER truncations of stationary records
  (seeded, reproducible; `--T0 1` gives an uncorrelated control).
- `tools/gci.py` exposes the three-grid procedure as a function and also prints the index for a
  formal order of one; `tools/replayGate.py` keeps the logged mean, rms and `Neff`.

- `tools/statConv.py` now follows the kernel's conventions exactly: mean sampling interval of
  the whole record, no truncation below eight batches, autocorrelation normalised by the
  unbiased variance, and a symmetric Hann window in the Welch estimate. It agrees with the
  compiled kernel to the printed precision.
- New `tools/verifyThirdParty.py`: a check that takes the autocorrelation from statsmodels, in
  addition to the `scipy.stats` quantiles and `scipy.signal.welch` periodogram that
  `tools/statConv.py` already uses.
- New `tools/replayGate.py`: replays the stopping rule from a solver log, with options for the
  tolerance, the gated signals, the settling time, the number of successive checks and a
  requirement that the minimum of the transient rule be interior.
- New `tools/genUduct.py`: the parametric mesh generator of the U-duct case.
- Documentation: the description of the reported quantities and of the stopping rule now
  matches the code, including the role of `assessAfter`, the tolerance scale
  `max(|mean|, refScale)`, the time units of `win` and the search-limit caveat; the limitations
  now describe the undercoverage on short records, the restart behavior, the zero signal of a
  missing field, the role of `writeControl` and of `refScale` in the stopping rule; code comments attribute the truncation of the autocorrelation sum to the
  first-zero-crossing rule; the README cites the Zenodo concept record.

## 1.0.0

First public release.

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
