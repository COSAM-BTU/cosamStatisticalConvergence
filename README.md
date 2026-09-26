# cosamStatisticalConvergence

An OpenFOAM function object that assesses the **statistical convergence of time-averaged
quantities while the solver runs**, and can stop the run once the bulk statistics have
settled.

Unsteady simulations are commonly stopped when residuals stagnate, when a preset number of
time steps has elapsed, or when the field stops changing to the eye. None of these measures
the sampling uncertainty of the averages that are actually reported. This component measures
it in situ: for every monitored signal it locates the end of the initial transient, estimates
how many statistically independent samples the record really contains, and reports confidence
intervals for the mean and for the root mean square at every write.

## Features

- **Initial-transient detection** by the batched marginal-standard-error rule (MSER-5).
- **Integral time scale** from the autocorrelation, truncated at the first non-positive lag,
  and the **effective sample size** `Neff` that follows from it.
- **Small-sample confidence intervals**: Student-*t* for the mean and the (asymmetric)
  chi-squared interval for the root mean square. At the `Neff = O(5)` typical of turbulent
  records the Student-*t* interval on the mean is about half as wide again as the normal one
  (the normal quantile understates it by about a third), so the correct quantiles matter.
- **Welch periodogram** with a segment length that scales with the record.
- **Signals from the field**: either a single-cell probe or a volume-weighted spatial average
  of a component or of the magnitude. A spatial average has a far smaller variance, but the
  average of the velocity magnitude is nearly fixed by the flow rate, so gate on the quantities
  you report (see the notes on the stopping rule).
- **Optional auto-stop**: when two consecutive averaging windows of a gated signal agree to
  within a tolerance, over a number of successive checks, the function object writes and
  stops the run through the solver time control. This is a plateau test, not an error target;
  see the notes below before relying on it.
- **A standalone verifier**, `statConvCheck`, compiled from the *same* kernel header as the
  function object, so the code path that is verified is the code path the solver executes.

## Requirements

- OpenFOAM v2406 (ESI), double precision.
- Python 3 with NumPy and SciPy, for the reference implementation, the tools and the AR(1)
  tutorial; statsmodels in addition for `tools/verifyThirdParty.py`.

## Installation

```bash
git clone https://github.com/COSAM-BTU/cosamStatisticalConvergence.git
cd cosamStatisticalConvergence/src
./Allwmake
```

This builds `libcosamStatisticalConvergence.so` into `$FOAM_USER_LIBBIN` and `statConvCheck`
into `$FOAM_USER_APPBIN`.

## Usage

Add the function object to `system/controlDict`:

```
functions
{
    statisticalConvergence1
    {
        type            statisticalConvergence;
        libs            (cosamStatisticalConvergence);

        writeControl    runTime;
        writeInterval   1;

        sampleInterval  1;        // sample every N time steps
        minSamples      100;      // no assessment below this many samples
        tolMean         0.01;     // target 95% interval on the mean, relative to |mean|
        tolRms          0.05;     // target 95% interval on the rms, relative to rms

        // --- optional automatic stopping ---
        autoStop        true;
        assessAfter     15;       // do not assess before this time
        compareWindow   8;        // length of each averaging window, in time units
        tolPlateau      0.05;     // window-to-window agreement, times refScale
        refScale        1;        // scale the tolerance is measured against
        stopConsecutive 2;        // successive checks that must agree

        signals
        {
            volUmag
            {
                mode        volAverage;   // volume-weighted average over the domain
                field       U;
                component   -1;           // -1 = magnitude
            }

            probeUy
            {
                mode        point;        // single-cell probe (the default)
                field       U;
                component   1;
                location    (0.3 0.3 0.5);
                autoStopGate false;       // reported, but excluded from the stop decision
            }
        }
    }
}
```

Each signal is reported at every evaluation with the time, the retained sample count, the
transient end, the mean, the rms, `Neff`, the two intervals, the period of the largest Welch
bin and the quantities of the stopping test. `statConvCheck` additionally prints `T0`, the
integral time and the `samplingAdequate` flag. See `doc/interface.md` for every entry.

**Notes on the stopping rule.** The rule is a plateau test of the running mean; it does not
use the intervals.

- Gate on the quantities that will be reported, or on bulk quantities that respond to the flow
  features of interest, and use more than one gated signal. The volume average of `|U|` is
  nearly fixed by the flow rate, so its plateau says little about local statistics, and a
  single signal whose transient is misjudged can trigger an early stop.
- Scale `tolPlateau` to the variability of the gated signals. A tolerance far above their
  window-to-window scatter never fails, and the stop is then fixed by the transient detection.
- A stop issued while the minimum of the transient rule sits on its half-record search limit is
  not statistically supported. Version 1.x does not test for this; check it, and replay the
  rule with other settings, from the log:
  `python3 tools/replayGate.py log.pimpleFoam --interior`. An interior minimum is necessary but
  not sufficient: with signals sampled every time step the statistical inefficiency can be
  hundreds of samples, far more than the batch size of five.
- `stopConsecutive` counts evaluations, and the evaluation interval follows the function
  object's `writeControl`. With `writeControl timeStep` the physical spacing, and hence the
  stopping time, depends on the time step; use `runTime` or `adjustableRunTime` to evaluate in
  simulated time.
- The plateau tolerance is relative to `max(|mean|, refScale)`. For a signal whose mean is much
  smaller than `refScale` (for example a kinematic pressure with `refScale 1`) the test is
  effectively absolute and easily passed; choose `refScale` for each gated signal accordingly.
- The per-signal intervals continue to be reported after the gate has fired, so a run stopped on
  a bulk mean still documents where local statistics remain uncertain.

## Tutorials

| Case | What it shows | Runtime |
|---|---|---|
| `tutorials/arSignal` | Kernel verification on a synthetic AR(1) process with a known integral time scale | < 1 s |
| `tutorials/cavityProbe` | The function object attached to a small transient case, with both signal modes | a few seconds |

```bash
cd tutorials/arSignal && ./Allrun
```

`arSignal` generates a seeded AR(1) record (`phi = 0.9`, so the statistical inefficiency has
the closed form `T0 = (1+phi)/(1-phi) = 19`) with a ramp over the first 2000 samples, then
analyses it. The expected output is a transient index of `1460`, a mean of `5.0479`, an rms
of `2.3075`, `T0 = 21.1` and `Neff = 1826`. The difference from the theoretical `19` lies
within the scatter of the estimator at this record length: over 400 stationary realizations
`tools/coverage.py --T0 19 --ratios 2028.4 --reps 400 --seed 2` gives a ratio of estimated to
true `T0` of `1.027 +- 0.088`, that is `19.5 +- 1.7`.

```bash
cd tutorials/cavityProbe && ./Allrun
```

`cavityProbe` is the standard lid-driven cavity with the function object attached, monitoring
both a volume average of `|U|` and a single-cell probe. It exists to show the dictionary and
the log output on a case that runs in seconds, not to demonstrate convergence of a turbulent
flow.

## Repository layout

```
src/statisticalConvergence/   the function object; statConvMath.H is the shared kernel
src/statConvCheck/            standalone verifier built from the same kernel
tutorials/                    the two cases above
tools/statConv.py             separate Python implementation of the kernel
tools/verifyThirdParty.py     check against statsmodels (autocorrelation) and scipy
tools/coverage.py             coverage of the intervals on synthetic AR(1) records
tools/replayGate.py           replay and audit the stopping rule from a solver log
tools/genUduct.py             blockMeshDict generator for the U-duct case
tools/gci.py                  grid-convergence index (Celik et al. 2008)
tools/fftPeriod.py            offline amplitude spectrum of a probe record
doc/interface.md              every dictionary entry and every reported quantity
```

## Verification

`statConvCheck`, which executes the same code as the function object, is compared on identical
input with `tools/statConv.py`, a separate Python implementation by the same author that shares
no source with the C++ and takes its quantiles from `scipy.stats` and its periodogram from
`scipy.signal.welch`, and with `tools/verifyThirdParty.py`, which in addition takes the
autocorrelation from statsmodels. On the AR(1) record of `tutorials/arSignal` and on the
deposited U-duct probe record they agree on every reported quantity to about 3 parts in 10^9,
the precision to which `statConvCheck` prints. The MSER-5 truncation agrees with pyMSER
(Oliveira et al., J. Chem. Theory Comput. 20:8559, 2024) when its search is restricted to the
first half of the record.

Agreement between implementations does not show that the intervals cover at the nominal rate.
`tools/coverage.py` measures the coverage on stationary AR(1) records. With `T0 = 468` samples
and records of 1 to 4.5 `T0` the tool reports a median `Neff` of 5.5 to 8.7, and the 95%
interval of the mean covers 62 to 83% (42 to 70% after MSER truncation); the coverage reaches
93% at 50 `T0` and the nominal 95% only near 100 `T0`. Treat a small `Neff` as a likely
overestimate.

The U-duct case dictionaries, the probe records and the data behind the examples are deposited
separately on Zenodo (doi:10.5281/zenodo.22880129).

## Citation

If you use this software, please cite the Zenodo archive. A paper describing the software is in
preparation; this section will be updated when it is published.

```bibtex
@misc{aydinbakar26statconvcode,
    author       = {Ayd{\i}nbakar, Levent},
    title        = {{cosamStatisticalConvergence}: in-situ statistical-convergence
                    assessment for {OpenFOAM}},
    year         = {2026},
    howpublished = {Zenodo},
    note         = {Concept identifier, resolving to the latest version},
    doi          = {10.5281/zenodo.22880057},
}
```

## License

GNU General Public License v3.0; see `LICENSE`. This is the license of OpenFOAM itself, which
this library extends.

## Author

Levent Aydınbakar, Department of Mechanical Engineering, Bursa Technical University, Bursa,
Türkiye. ORCID [0000-0002-8820-9874](https://orcid.org/0000-0002-8820-9874).

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the
OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trademarks.
