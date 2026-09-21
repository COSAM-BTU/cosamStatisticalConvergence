#!/usr/bin/env python3
"""
statisticalConvergence -- REFERENCE (Python) implementation of the decision method
that the OpenFOAM C++ functionObject will implement. Given a scalar time signal it:
  1. detects the initial transient (MSER-m, White 1997 / MSER-5),
  2. computes the autocorrelation, integral time scale, effective sample size Neff,
  3. autocorrelation-corrected 95% CI for the MEAN and the RMS,
  4. dominant frequency via Welch PSD,
  5. reports stationary? / samplingAdequate? against user tolerances.

This is the independent reference implementation against which the compiled C++
kernel (src/statisticalConvergence/statConvMath.H) is verified: the two share no
source and are compared on identical input. Method references: Oliver et al.,
Phys. Fluids 26:035101 (2014); White, Simulation 69:323 (1997); Spratt (1998) and
White, Cobb & Spratt (2000) for the batched MSER-5 form; Geyer, Statist. Sci.
7:473 (1992) for the truncation rule; Welch, IEEE Trans. Audio Electroacoust.
15:70 (1967) for the periodogram.

Usage:
  self-test:   python3 statConv.py --selftest
  probe file:  python3 statConv.py <openfoam_probes_file> [--comp Uy] [--tolMean 0.01] [--tolRms 0.05]
"""
import sys, math
import numpy as np

try:
    from scipy.signal import welch
    from scipy.stats import t as _tdist, chi2 as _chi2
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


def mser(x, m=5):
    """MSER-m initial-transient truncation. Returns truncation index in ORIGINAL samples."""
    x = np.asarray(x, float)
    nb = len(x) // m
    if nb < 4:
        return 0
    b = x[:nb * m].reshape(nb, m).mean(axis=1)
    n = len(b)
    best_d, best_val = 0, math.inf
    for d in range(0, n // 2):          # don't truncate more than half
        y = b[d:]
        z = len(y)
        val = y.var(ddof=0) / z          # MSER(d) = sum(y-ybar)^2 / z^2 = var/z
        if val < best_val:
            best_val, best_d = val, d
    return best_d * m


def acf(x, maxlag):
    """Biased autocorrelation rho(0..maxlag)."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    v = np.dot(x, x) / n
    if v <= 0:
        return np.zeros(maxlag + 1)
    r = np.correlate(x, x, "full")[n - 1: n - 1 + maxlag + 1] / (v * n)
    return r


def integral_time_scale(x, dt):
    """T0 (statistical inefficiency, dimensionless) = 1 + 2*sum rho(k) up to first zero
    crossing; physical integral time tau = T0*dt/2. Returns (T0, tau)."""
    n = len(x)
    maxlag = min(n - 2, max(50, n // 4))
    r = acf(x, maxlag)
    s = 0.0
    for k in range(1, len(r)):
        if r[k] <= 0:                    # truncate at first non-positive lag (Sokal)
            break
        s += r[k]
    T0 = 1.0 + 2.0 * s
    return T0, T0 * dt / 2.0


def analyze(t, x, tol_mean=0.01, tol_rms=0.05):
    t = np.asarray(t, float); x = np.asarray(x, float)
    n0 = len(x)
    dt = np.median(np.diff(t)) if n0 > 1 else 1.0

    d = mser(x)                          # transient end (samples)
    xs = x[d:]; ts = t[d:]
    n = len(xs)
    if n < 20:
        return {"error": "too few stationary samples", "transientEnd_idx": d}

    mean = float(np.mean(xs))
    rms = float(np.std(xs, ddof=1))      # rms of fluctuation
    T0, tau = integral_time_scale(xs, dt)
    Neff = n / T0 if T0 > 0 else float(n)

    # ACF-corrected 95% CIs, small-sample quantiles (mirrors statConvMath.H exactly):
    # Student-t for the mean, chi-squared (asymmetric) for the rms.
    df = max(Neff - 1.0, 0.5)
    tq = _tdist.ppf(0.975, df)
    hw_mean = tq * rms / math.sqrt(max(Neff, 1e-9))
    sig_lo = rms * math.sqrt(df / _chi2.ppf(0.975, df))
    sig_hi = rms * math.sqrt(df / _chi2.ppf(0.025, df))
    hw_rms = 0.5 * (sig_hi - sig_lo)
    rel_mean = hw_mean / abs(mean) if abs(mean) > 1e-12 else float("inf")
    rel_rms = hw_rms / rms if rms > 1e-12 else float("inf")
    rel_rms_lo = (rms - sig_lo) / rms if rms > 1e-12 else float("inf")
    rel_rms_hi = (sig_hi - rms) / rms if rms > 1e-12 else float("inf")

    # dominant frequency (Welch)
    fdom = None
    if HAVE_SCIPY and n >= 64:
        fs = 1.0 / dt
        L = min(max(256, n // 8), n // 2)   # matches statConvMath.H
        f, P = welch(xs - mean, fs=fs, nperseg=L)
        if len(f) > 1:
            fdom = float(f[1 + int(np.argmax(P[1:]))])

    stationary = (d < n0 - 20)
    sampling_ok = (rel_mean < tol_mean) and (rel_rms < tol_rms)

    return {
        "n_total": n0, "dt": dt, "transientEnd_idx": d, "transientEnd_t": float(t[d]),
        "n_stationary": n, "mean": mean, "rms": rms,
        "T0_statIneff": T0, "integralTime": tau, "Neff": Neff,
        "CI95_mean_halfwidth": hw_mean, "CI95_mean_rel": rel_mean,
        "CI95_rms_halfwidth": hw_rms, "CI95_rms_rel": rel_rms,
        "CI95_rms_rel_lo": rel_rms_lo, "CI95_rms_rel_hi": rel_rms_hi,
        "dominantFreq": fdom,
        "stationary": stationary, "samplingAdequate": sampling_ok,
        "tol_mean": tol_mean, "tol_rms": tol_rms,
    }


def read_probe(path, comp="Uy"):
    """Read an OpenFOAM probes file. Returns (t, value) for probe 0, chosen component.
    Vector comps: Ux/Uy/Uz; scalar: value."""
    idx = {"Ux": 0, "Uy": 1, "Uz": 2, "mag": None, "value": 0}.get(comp, 1)
    t, v = [], []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            # replace ( ) with spaces to split vector components
            a = ln.replace("(", " ").replace(")", " ").split()
            if len(a) < 2:
                continue
            tt = float(a[0]); rest = [float(z) for z in a[1:]]
            t.append(tt)
            if comp == "mag" and len(rest) >= 3:
                v.append(math.sqrt(rest[0] ** 2 + rest[1] ** 2 + rest[2] ** 2))
            else:
                v.append(rest[idx] if idx < len(rest) else rest[0])
    return np.array(t), np.array(v)


def _selftest():
    # AR(1): x_i = phi x_{i-1} + eps.  Theory: T0 = (1+phi)/(1-phi), tau = T0*dt/2.
    rng = np.random.default_rng(0)
    phi, N, dt = 0.9, 40000, 0.01
    x = np.zeros(N)
    for i in range(1, N):
        x[i] = phi * x[i - 1] + rng.standard_normal()
    x += 5.0                              # mean = 5
    x[:2000] += np.linspace(10, 0, 2000)  # add an initial transient to test MSER
    t = np.arange(N) * dt
    r = analyze(t, x)
    T0_theory = (1 + phi) / (1 - phi)
    print("AR(1) self-test (phi=0.9):")
    print(f"  T0  estimated = {r['T0_statIneff']:.2f}   theory = {T0_theory:.2f}")
    print(f"  transient end idx = {r['transientEnd_idx']} (injected ~2000)")
    print(f"  mean = {r['mean']:.3f} (theory 5.0)   Neff = {r['Neff']:.0f}  of {r['n_stationary']}")
    print(f"  CI95 mean rel = {r['CI95_mean_rel']*100:.2f}%   rms rel = {r['CI95_rms_rel']*100:.2f}%")
    print(f"  stationary={r['stationary']} samplingAdequate={r['samplingAdequate']}  scipy={HAVE_SCIPY}")
    ok = abs(r["T0_statIneff"] - T0_theory) / T0_theory < 0.25 and abs(r["mean"] - 5) < 0.1
    print("  SELFTEST", "PASS" if ok else "CHECK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _selftest(); sys.exit(0)
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    comp = "Uy"; tolM, tolR = 0.01, 0.05
    for i, a in enumerate(sys.argv):
        if a == "--comp": comp = sys.argv[i + 1]
        if a == "--tolMean": tolM = float(sys.argv[i + 1])
        if a == "--tolRms": tolR = float(sys.argv[i + 1])
    t, v = read_probe(args[0], comp)
    r = analyze(t, v, tolM, tolR)
    print(f"# statConv  {args[0]}  comp={comp}")
    for k, val in r.items():
        print(f"  {k:22s} {val}")
