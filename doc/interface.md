# Interface reference

Every dictionary entry the `statisticalConvergence` function object reads, and every quantity
it reports. Defaults are those compiled into the code.

## Top-level entries

| Entry | Type | Default | Meaning |
|---|---|---|---|
| `type` | word | — | `statisticalConvergence` |
| `libs` | list | — | `(cosamStatisticalConvergence)` |
| `sampleInterval` | label | `1` | Record each signal every `N` time steps. |
| `minSamples` | label | `100` | No statistics are formed below this many recorded samples. |
| `tolMean` | scalar | `0.01` | Target for the relative 95% half-width of the mean. Used only for the reported `samplingAdequate` flag. |
| `tolRms` | scalar | `0.05` | Target for the relative 95% interval of the rms. Used only for the reported flag. |
| `signals` | dict | — | Sub-dictionary; one entry per signal, keyed by the signal name. |

## Automatic stopping

| Entry | Type | Default | Meaning |
|---|---|---|---|
| `autoStop` | bool | `false` | Enable the plateau criterion and the write-and-stop. |
| `assessAfter` | scalar | `0` | Simulated time before which no assessment is made, so that a start-up transient cannot trigger a stop. |
| `compareWindow` | scalar | `10` | Length `w` of each averaging window, in time units. |
| `tolPlateau` | scalar | `0.05` | Two consecutive disjoint windows count as agreeing when their means differ by less than `tolPlateau * refScale`. |
| `refScale` | scalar | `1` | Scale the plateau tolerance is measured against, for example a bulk velocity. |
| `stopConsecutive` | label | `2` | Number of successive checks that must agree before the run is stopped. |

The gate becomes evaluable two comparison windows after the detected transient end, since two
disjoint windows of length `w` must both lie inside the stationary record.

Scale `tolPlateau * refScale` to the variability of the gated signal, not to a global
reference, if the plateau test is to be an active constraint rather than a permissive one. A
tolerance set far above the window-to-window scatter of the signal is satisfied at the first
evaluation and never fails, in which case the stopping time is fixed by the transient
detection and the two-window confirmation rather than by the plateau test.

## Per-signal entries

Each entry of `signals` is a dictionary whose keyword is the signal name.

| Entry | Type | Default | Meaning |
|---|---|---|---|
| `mode` | word | `point` | `point` for a single-cell probe, `volAverage` for a volume-weighted average over the whole domain. |
| `field` | word | — | Name of a registered `volScalarField` or `volVectorField`. |
| `component` | label | `0` | Vector component. In `volAverage` mode, `-1` selects the magnitude. |
| `location` | point | — | Probe position; required in `point` mode. The cell containing it is found once, at construction. |
| `autoStopGate` | bool | `true` | Whether this signal participates in the stop decision. Signals with `false` are still reported. |

In parallel, a probe value is recovered by a reduction across ranks, so exactly one rank owns
the cell and all ranks see the same value. A volume average is formed from the local
volume-weighted sums and reduced.

## Reported quantities

Written to the log at each write interval, per signal.

| Name | Meaning |
|---|---|
| `t` | Current time. |
| `nStat` | Samples retained after the transient is discarded. |
| `transientEnd` | Time at which the initial transient is judged to end. |
| `mean`, `rms` | Mean and fluctuation root mean square of the retained record. |
| `Neff` | Effective number of independent samples, `nStat / T0`. |
| `CImean` | Relative 95% half-width of the mean, from the Student-*t* quantile at `Neff - 1` degrees of freedom. |
| `CIrms` | Relative 95% interval of the rms, from the chi-squared distribution. Reported as an asymmetric pair, because the interval is not symmetric at small `Neff`. |
| `period` | Reciprocal of the lowest resolvable Welch frequency, that is the resolution limit of the spectral estimate on the present record, not a resolved flow period. |
| `win`, `plateauRel` | Length of the comparison window actually used, and the window-to-window difference as a fraction of `refScale`. |
| `converged` | Whether this signal currently satisfies the plateau criterion. |

## Method

1. **Transient.** The record is grouped into non-overlapping batches of `m = 5` samples with
   batch means `b_j`. The truncation index minimises the biased variance of the retained
   batches divided by their number, searched over the first half of the record only. The
   first `d* m` samples are discarded.
2. **Moments.** Mean and rms of the retained record, the rms with the unbiased normalisation.
3. **Integral time scale.** The autocorrelation is summed to the first non-positive lag,
   giving the statistical inefficiency `T0 = 1 + 2 sum rho(k)`; the physical integral time is
   `tau = T0 dt / 2`.
4. **Effective sample size.** `Neff = nStat / T0`.
5. **Intervals.** Student-*t* for the mean and chi-squared for the rms, both at `Neff - 1`
   degrees of freedom. The quantiles are evaluated by bisection on the incomplete beta and
   incomplete gamma functions respectively, so non-integer degrees of freedom are handled.
6. **Spectrum.** Welch periodogram with a Hann window, fifty per cent overlap and per-segment
   detrending. The segment length is `min(max(256, nStat/8), nStat/2)`, so the resolution
   improves as the record grows.

## Limitations

- The full history of every monitored signal is retained and the autocorrelation is recomputed
  at each evaluation, so memory and work grow with the record. This is negligible for the
  record lengths of a typical campaign but is not an asymptotic statement.
- The integral time scale uses a truncated autocorrelation sum, which carries a modest positive
  bias relative to an autoregressive fit.
- The adequacy tolerances are relative to the local mean, which is stringent and effectively
  unreachable for a near-zero-mean signal such as a reversed-flow probe. That is why the
  stopping gate uses a plateau on a low-variance signal rather than the point interval.
- The transient index returned by the marginal-standard-error rule is not invariant to the
  length of the record: on a longer record of the same signal the rule may select a later
  truncation. Report the transient index together with the record it was obtained from.
