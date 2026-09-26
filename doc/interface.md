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
| `tolMean` | scalar | `0.01` | Target for the relative 95% half-width of the mean. Used only for the `samplingAdequate` flag printed by `statConvCheck`; the function object reports the intervals themselves. |
| `tolRms` | scalar | `0.05` | Target for the nominal relative 95% half-width of the rms interval. Used only for the same flag. |
| `signals` | dict | — | Sub-dictionary; one entry per signal, keyed by the signal name. |

## Automatic stopping

| Entry | Type | Default | Meaning |
|---|---|---|---|
| `autoStop` | bool | `false` | Enable the plateau criterion and the write-and-stop. |
| `assessAfter` | scalar | `0` | Simulated time before which no stop can be issued. It does not constrain the comparison windows, which depend only on the detected transient. |
| `compareWindow` | scalar | `10` | Length `w` of each averaging window, in time units. |
| `tolPlateau` | scalar | `0.05` | Two consecutive windows count as agreeing when their means differ by less than `tolPlateau * max(|mean|, refScale)`. |
| `refScale` | scalar | `1` | Scale the plateau tolerance is measured against, for example a bulk velocity. |
| `stopConsecutive` | label | `2` | Number of successive checks that must agree before the run is stopped. |

A gated signal is ready when its retained record holds two windows; the two windows are the
last `2 w` of the record, converted to samples with the mean sampling interval of the whole
record, so the stop cannot come earlier than two windows after the detected transient end. A
counter increases at every evaluation at which `t >= assessAfter` and all gated signals agree,
is reset otherwise, and the run is stopped when it reaches `stopConsecutive`.

When the minimum of the transient rule lies on its half-record search limit, the rule has not
found the end of the transient. Version 1.x does not test for this before stopping. The log
contains everything needed to check it after the fact:
`python3 tools/replayGate.py log.pimpleFoam --interior` replays the rule with an interior
minimum required, and the other options of the script replay it with a different tolerance,
set of gated signals, settling time or number of successive checks. Use more than one gated
signal: a single signal whose transient is misjudged can otherwise trigger an early stop.

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

Written to the log at each evaluation, that is at the write interval of the function object, per signal.

| Name | Meaning |
|---|---|
| `t` | Current time. |
| `nStat` | Samples retained after the transient is discarded. |
| `transientEnd` | Time at which the initial transient is judged to end. |
| `mean`, `rms` | Mean and fluctuation root mean square of the retained record. |
| `Neff` | Effective number of independent samples, `nStat / T0`. |
| `CImean` | Relative 95% half-width of the mean, from the Student-*t* quantile at `Neff - 1` degrees of freedom. |
| `CIrms` | Relative 95% interval of the rms, from the chi-squared distribution. Reported as an asymmetric pair, because the interval is not symmetric at small `Neff`. |
| `period` | Reciprocal of the frequency of the largest Welch bin, excluding the zero-frequency and Nyquist bins. On records that are short compared with the slowest motion this is the first bin, that is the resolution limit, not a resolved flow period. |
| `win`, `plateauRel` | Span of the two comparison windows in time units, `2 n_w dt` (zero while the retained record holds fewer than two windows), and the window-to-window difference as a fraction of `max(|mean|, refScale)`. |
| `converged` | Whether this signal currently satisfies the plateau criterion. |

## Method

1. **Transient.** The record is grouped into non-overlapping batches of `m = 5` samples with
   batch means `b_j`. The truncation index minimises the biased variance of the retained
   batches divided by their number, searched over `d < n_b/2`; ties go to the smaller `d`,
   and no truncation is attempted below eight batches. The first `d* m` samples are
   discarded.
2. **Moments.** Mean and rms of the retained record, the rms with the unbiased normalisation.
3. **Integral time scale.** The autocorrelation, normalised by `nStat` times the unbiased
   variance, is summed from lag 1 to the last lag before it first becomes non-positive, with
   the lag capped at `min(nStat - 2, max(50, nStat/4))`. This gives the statistical
   inefficiency `T0 = 1 + 2 sum rho(k)`; the physical integral time is `tau = T0 dt / 2`, with
   `dt` the mean sampling interval of the whole record.
4. **Effective sample size.** `Neff = nStat / T0`.
5. **Intervals.** Student-*t* for the mean and chi-squared for the rms, both at `Neff - 1`
   degrees of freedom. The quantiles are evaluated by bisection on the incomplete beta and
   incomplete gamma functions respectively, so non-integer degrees of freedom are handled.
6. **Spectrum.** Welch periodogram with a symmetric Hann window, fifty per cent overlap and
   removal of each segment's mean. The segment length is `min(max(256, nStat/8), nStat/2)`,
   so the resolution improves as the record grows.

## Limitations

- The full history of every monitored signal is retained and the autocorrelation is recomputed
  at each evaluation, so memory and work grow with the record. This is negligible for the
  record lengths of a typical campaign but is not an asymptotic statement.
- On records that span only a few integral time scales the truncation of the autocorrelation
  sum underestimates `T0`, so `Neff` is overstated and the intervals undercover: with
  `tools/coverage.py`, records of 1 to 4.5 `T0` report a median `Neff` of 5.5 to 8.7, the 95%
  interval of the mean covers 62 to 83% of them (42 to 70% after MSER truncation), and the
  coverage reaches 93% at 50 `T0` and the nominal level only near 100 `T0`. On long records the
  chi-squared interval of the rms is conservative, because the variance estimate has about
  `2 Neff` effective degrees of freedom.
- The adequacy tolerances are relative to the local mean, which is stringent and effectively
  unreachable for a near-zero-mean signal such as a reversed-flow probe. The stopping rule is a
  plateau test and does not use the intervals.
- `stopConsecutive` counts evaluations, and the evaluation interval follows the function
  object's `writeControl`: with `timeStep` their physical spans depend on the time step, with
  `runTime` or `adjustableRunTime` they are fixed in simulated time.
- `component -1` (magnitude) is available for `volAverage` signals only; point signals need a
  component index of 0, 1 or 2.
- The histories are held in memory. They are not saved across restarts, and they start afresh
  when the function-object dictionary is re-read during a run.
- A missing field or a probe outside the mesh gives a constant zero signal, which passes the
  plateau test; only a location outside the mesh triggers a warning. Check the gated signals in
  the log at the start of a run.
- The transient index returned by the marginal-standard-error rule is not invariant to the
  length of the record: on a longer record of the same signal the rule may select a later or an
  earlier truncation, and its minimum may sit on the search limit. Report the transient index
  together with the record it was obtained from.
- The effective sample size itself depends on the record length: until the record spans many
  multiples of the slowest time scale of the signal, the integral time scale is under-estimated
  and `Neff` over-stated. Do not project the averaging time needed for a target precision from
  a short record.
