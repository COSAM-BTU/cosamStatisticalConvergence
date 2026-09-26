#!/usr/bin/env python3
"""
Check of the statistical kernel against third-party routines.

The compiled kernel (src/statisticalConvergence/statConvMath.H, also used by the
standalone verifier statConvCheck) and the Python implementation tools/statConv.py were
written by the same author. tools/statConv.py already takes its quantiles from scipy.stats
and its periodogram from scipy.signal.welch; this script adds an autocorrelation computed
by statsmodels, so that every stage except the MSER-5 truncation is checked against a
third-party routine. The truncation index is borrowed from tools/statConv.py here; it
agrees with pyMSER (Oliveira et al., J. Chem. Theory Comput. 20:8559, 2024) when the pyMSER
search is restricted to the first half of the record, as in the kernel.

Conventions reproduced from the kernel, so that differences are genuine:
  * sampling interval dt = (t_last - t_first)/(n - 1);
  * autocorrelation summed from lag 1 to the first non-positive lag, capped at
    max(50, n/4) and n - 2; the kernel normalises the lag-k sum by n times the
    unbiased variance, statsmodels (adjusted=False) by the biased variance, so the
    two differ by the exact factor (n-1)/n, which is applied and reported;
  * Welch: symmetric Hann window, segment length L = min(max(256, n/8), n/2),
    50 % overlap, per-segment mean removal, largest bin excluding k = 0 and k = L/2.

Usage:  python3 verifyThirdParty.py <two-column t,value file> [--tolMean 0.01 --tolRms 0.05]
"""
import sys
import math

import numpy as np
from scipy import signal, stats
from statsmodels.tsa.stattools import acf

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from statConv import mser  # noqa: E402  (only the truncation index is borrowed)


def third_party(t, x):
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    n0 = len(x)
    dt = (t[-1] - t[0]) / (n0 - 1)
    d = mser(x)
    xs = x[d:]
    n = len(xs)
    mean = float(np.mean(xs))
    rms = float(np.std(xs, ddof=1))

    maxlag = min(n - 2, max(50, n // 4))
    rho_sm = acf(xs, nlags=maxlag, adjusted=False, fft=True)
    rho = rho_sm * (n - 1) / n            # kernel normalisation (unbiased variance)
    s = 0.0
    for k in range(1, len(rho)):
        if rho[k] <= 0:
            break
        s += rho[k]
    T0 = 1.0 + 2.0 * s
    s_sm = 0.0
    for k in range(1, len(rho_sm)):
        if rho_sm[k] <= 0:
            break
        s_sm += rho_sm[k]
    T0_sm = 1.0 + 2.0 * s_sm

    Neff = n / T0
    df = max(Neff - 1.0, 0.5)
    hw_mean = stats.t.ppf(0.975, df) * rms / math.sqrt(Neff)
    sig_lo = rms * math.sqrt(df / stats.chi2.ppf(0.975, df))
    sig_hi = rms * math.sqrt(df / stats.chi2.ppf(0.025, df))

    L = min(max(256, n // 8), n // 2)
    f, P = signal.welch(xs - mean, fs=1.0 / dt, window=signal.windows.hann(L, sym=True),
                        nperseg=L, noverlap=L - max(1, L // 2), detrend="constant",
                        scaling="density")
    kmax = L // 2                          # kernel ignores the Nyquist bin
    fdom = float(f[1 + int(np.argmax(P[1:kmax]))])

    return {
        "transientEndIdx": d, "nStat": n, "mean": mean, "rms": rms,
        "T0": T0, "T0_statsmodelsNorm": T0_sm, "integralTime": T0 * dt / 2.0,
        "Neff": Neff, "ciMeanRel": hw_mean / abs(mean),
        "ciRmsRelLo": (rms - sig_lo) / rms, "ciRmsRelHi": (sig_hi - rms) / rms,
        "fdom": fdom,
    }


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    data = np.loadtxt(args[0], comments="#")
    r = third_party(data[:, 0], data[:, 1])
    print(f"# third-party channel (statsmodels acf, scipy.stats, scipy.signal.welch): {args[0]}")
    for k, v in r.items():
        print(f"  {k:20s} {v}")
