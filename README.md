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
  of a component or of the magnitude. The spatial average has a far smaller variance and is
  what makes an automatic stop reproducible.
- **Optional auto-stop**: when two consecutive averaging windows of a gated signal agree to
  within a tolerance, over a number of successive checks, the function object writes and
  stops the run through the solver time control.
- **A standalone verifier**, `statConvCheck`, compiled from the *same* kernel header as the
  function object, so the code path that is verified is the code path the solver executes.

## Requirements

- OpenFOAM v2406 (ESI), double precision.
- Python 3 with NumPy and SciPy, for the reference implementation and the AR(1) tutorial.

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

Each signal is reported at every write with its transient end, retained sample count, mean,
rms, integral time scale, `Neff`, the two intervals, the first resolvable spectral frequency
and the stationarity and adequacy flags.

**A note on the stopping gate.** Put the gate on a low-variance signal, that is on a spatial
average rather than on a single point. A point signal in a separated flow keeps wandering,
its windowed mean never presents a clean plateau, and the resulting stop is erratic. The
per-signal intervals continue to be reported for every signal after the gate has fired, so a
run that is stopped on the bulk mean still documents where local statistics remain uncertain.

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
of `2.3075`, `T0 = 21.1` and `Neff = 1826`. The small positive bias of the truncated-sum
estimator relative to the theoretical `19` is a known property of the direct estimator.

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
tools/statConv.py             independent Python reference implementation
tools/gci.py                  grid-convergence index (Richardson/Celik)
tools/fftPeriod.py            offline amplitude spectrum of a probe record
doc/interface.md              every dictionary entry and every reported quantity
```

## Verification

`tools/statConv.py` is the independent reference against which the compiled kernel is
verified: it shares no source with the C++ and is compared on identical input. On the AR(1)
record of `tutorials/arSignal` and on the deposited U-duct probe record, the
two agree to nine significant figures on the mean and the rms, to about one part in two
thousand on the integral time scale and the effective sample size, and to two significant
figures on the Welch frequency.

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
