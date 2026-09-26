#!/usr/bin/env python3
"""
Coverage of the statConv intervals on synthetic AR(1) records.

Each realization is a stationary first-order autoregressive record
x_i = phi*x_{i-1} + a_i with unit Gaussian innovations a_i, started from the stationary
distribution. Its statistical inefficiency is T0 = (1 + phi)/(1 - phi), its mean is zero
and its standard deviation is (1 - phi^2)^(-1/2). The estimators of statConv.py are
applied to every record, with and without MSER-5 truncation, and the script reports for
each record length N/T0:
  - the median (and the mean and standard deviation) of the ratio of estimated to true T0,
  - the fraction of records whose 95% interval of the mean contains the true mean, with
    the estimated T0 (as the tool does) and with the true T0,
  - the same fraction for the chi-squared interval of the standard deviation,
  - the median effective sample size N/T0hat that the tool would report,
  - how often MSER truncates more than 20% of the stationary record with its minimum inside
    the search range, and how often the minimum lies on the search limit.
With --T0 1 the records are uncorrelated (phi = 0), which serves as a control.

The autocorrelation and the MSER statistic are evaluated with an FFT and with cumulative
sums of mean-removed batch means. Both reproduce statConv.acf and statConv.mser, which is
checked at start-up on AR(1) records, at a cost that permits thousands of realizations of
long records.

Usage:
  python3 coverage.py                                            # Table 2 (T0 = 468, seed 1)
  python3 coverage.py --T0 19 --ratios 2028.4 --reps 400 --seed 2  # AR(1) signal of Table 1
  python3 coverage.py --T0 1 --ratios 2106 --reps 2000 --seed 3    # uncorrelated control
Options: --T0, --ratios (list of N/T0), --reps, --seed, --csv <file>
"""
import argparse
import math
import sys

import numpy as np
from scipy.signal import lfilter
from scipy.stats import chi2, t as tdist

import statConv


def acf_fft(x, maxlag):
    """statConv.acf evaluated with an FFT (same normalization)."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    v = np.dot(x, x) / (n - 1)
    if v <= 0:
        return np.zeros(maxlag + 1)          # as statConv.acf for a constant record
    nfft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, nfft)
    r = np.fft.irfft(f * np.conj(f), nfft)[: maxlag + 1]
    return r / (v * n)


def t0_hat(x):
    """statConv.integral_time_scale (T0 only) with the FFT autocorrelation."""
    n = len(x)
    maxlag = min(n - 2, max(50, n // 4))
    r = acf_fft(x, maxlag)
    neg = np.nonzero(r[1:] <= 0)[0]
    k = neg[0] + 1 if len(neg) else len(r)
    return 1.0 + 2.0 * float(np.sum(r[1:k]))


def mser_fast(x, m=5):
    """statConv.mser with cumulative sums: minimizes var(b_d..)/(n_b - d), d < n_b//2."""
    nb = len(x) // m
    if nb < 8:
        return 0
    b = np.asarray(x[: nb * m], float).reshape(nb, m).mean(axis=1)
    b = b - b.mean()                       # conditioning of the one-pass variance
    s1 = np.cumsum(b[::-1])[::-1]          # sum of b[d:]
    s2 = np.cumsum((b * b)[::-1])[::-1]    # sum of b[d:]^2
    z = nb - np.arange(nb)
    var = np.maximum(s2 / z - (s1 / z) ** 2, 0.0)
    val = var[: nb // 2] / z[: nb // 2]
    return int(np.argmin(val)) * m         # argmin returns the first (smallest d) minimum


def self_check(rng):
    for n, phi in ((3000, 0.95), (1200, 0.99)):
        x = lfilter([1.0], [1.0, -phi], rng.standard_normal(n))
        lag = min(n - 2, max(50, n // 4))
        assert np.allclose(acf_fft(x, lag), statConv.acf(x, lag), rtol=0, atol=1e-12)
        assert abs(t0_hat(x) - statConv.integral_time_scale(x, 1.0)[0]) < 1e-9
        assert mser_fast(x) == statConv.mser(x)


def intervals(xs, T0):
    """Student-t interval of the mean and chi-squared interval of sigma, as statConv."""
    n = len(xs)
    mean = float(np.mean(xs))
    rms = float(np.std(xs, ddof=1))
    neff = n / T0
    df = max(neff - 1.0, 0.5)
    hw = tdist.ppf(0.975, df) * rms / math.sqrt(neff)
    lo = rms * math.sqrt(df / chi2.ppf(0.975, df))
    hi = rms * math.sqrt(df / chi2.ppf(0.025, df))
    return mean, hw, lo, hi


def run(T0, ratio, reps, rng):
    phi = (T0 - 1.0) / (T0 + 1.0)
    sigma = 1.0 / math.sqrt(1.0 - phi * phi)
    n = int(round(ratio * T0))
    res = {False: [], True: []}
    nb = n // 5
    for _ in range(reps):
        a = rng.standard_normal(n)
        a[0] *= sigma                       # stationary start
        x = lfilter([1.0], [1.0, -phi], a)
        for use_mser in (False, True):
            d = mser_fast(x) if use_mser else 0
            xs = x[d:]
            th = t0_hat(xs)
            mean, hw, lo, hi = intervals(xs, th)
            mean_t, hw_t, _, _ = intervals(xs, T0)
            res[use_mser].append((th / T0, abs(mean) <= hw, abs(mean_t) <= hw_t,
                                  lo <= sigma <= hi, len(xs) / th,
                                  d > 0.2 * n and d // 5 < nb // 2 - 1, d // 5 == nb // 2 - 1))
    out = []
    for use_mser in (False, True):
        a = np.array(res[use_mser], float)
        out.append(dict(T0=T0, ratio=ratio, N=n, mser=use_mser,
                        t0_median=float(np.median(a[:, 0])), t0_mean=float(np.mean(a[:, 0])),
                        t0_sd=float(np.std(a[:, 0], ddof=1)), cover_mean=float(np.mean(a[:, 1])),
                        cover_mean_trueT0=float(np.mean(a[:, 2])),
                        cover_sigma=float(np.mean(a[:, 3])),
                        neff_median=float(np.median(a[:, 4])),
                        trunc20_interior=float(np.mean(a[:, 5])), on_limit=float(np.mean(a[:, 6]))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--T0", type=float, default=468.0)
    ap.add_argument("--ratios", type=float, nargs="+", default=[1, 2, 3, 4.5, 10, 20, 50, 100])
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    self_check(rng)
    rows = []
    print(f"AR(1), T0 = {a.T0:g}, {a.reps} realizations, seed {a.seed}; Monte Carlo standard "
          f"error of a coverage p: sqrt(p(1-p)/{a.reps})")
    print("  N/T0  MSER  median N/T0hat  median T0hat/T0  mean+-sd T0hat/T0  cover(mean)"
          "  cover(mean, true T0)  cover(sigma)  trunc>20%  onLimit")
    for r in a.ratios:
        for row in run(a.T0, r, a.reps, rng):
            rows.append(row)
            print(f"{row['ratio']:6g}  {'yes' if row['mser'] else ' no'}   {row['neff_median']:13.2f}"
                  f"   {row['t0_median']:14.3f}   {row['t0_mean']:6.3f} +- {row['t0_sd']:5.3f}"
                  f"   {row['cover_mean']:10.3f}   {row['cover_mean_trueT0']:19.3f}"
                  f"   {row['cover_sigma']:11.3f}   {row['trunc20_interior']:9.3f}"
                  f"   {row['on_limit']:7.3f}")
    if a.csv:
        keys = list(rows[0].keys())
        with open(a.csv, "w") as f:
            f.write(",".join(keys) + "\n")
            for row in rows:
                f.write(",".join(str(row[k]) for k in keys) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
