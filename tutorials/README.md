# Tutorials

| Case | What it shows | Runtime |
|---|---|---|
| `arSignal` | Verification of the statistical kernel on a synthetic AR(1) process whose integral time scale is known in closed form | < 1 s |
| `cavityProbe` | The function object attached to a small transient case, exercising both signal modes | a few seconds |

Run either with `./Allrun` and clean with `./Allclean`.

## arSignal

`makeSignal.py` writes a seeded AR(1) record, `x_i = 0.9 x_{i-1} + eps_i`, shifted to a mean
of 5 and contaminated by a linear ramp over the first 2000 samples so that the
transient-detection stage is exercised. `statConvCheck` then analyses it.

Expected output: transient index `1460`, mean `5.0479`, rms `2.3075`, `T0 = 21.1`,
`Neff = 1826`. The theoretical statistical inefficiency is `(1+phi)/(1-phi) = 19`; the estimate
of 21.1 lies within the scatter of the estimator at this record length, `19.5 +- 1.7` over 400
realizations (`tools/coverage.py --T0 19 --ratios 2028.4 --reps 400 --seed 2`).

## cavityProbe

The standard lid-driven cavity, with the function object monitoring a volume-weighted average
of `|U|` and a single-cell probe below the lid. It demonstrates the dictionary and the log
output on a case that runs in seconds. It is not a convergence demonstration: the cavity
settles to a steady state, so the reported fluctuation level falls towards zero.
