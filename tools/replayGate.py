#!/usr/bin/env python3
"""
Replay the stopping rule of the statisticalConvergence function object from a solver log.

Every evaluation of the function object writes one line per signal, for example this line
of the production run on the medium grid of the accompanying paper, written by the
pre-release build:

  statConv[vol_Umag] t=32.390808 nStat=3255 transientEnd=14.612272 mean=1.0320192
      rms=0.0022574605 Neff=5.5818482 CImean=0.18146774% CIrms=58.661351%
      period=1.2758911 win=15.998478 plateauRel=0.11792219% converged=true

The pre-release build used the normal quantile and printed CIrms as a single value; release
1.1 uses the Student-t quantile and prints CIrms as an asymmetric pair, here
CImean=0.24468% and CIrms=[-38.6%,+159.9%]. Both formats are read. The rule of Section "Stopping rule" of the
accompanying paper is: a gating signal is ready when its retained record holds two
comparison windows (win > 0 in the log) and converged when, in addition, plateauRel is
below the tolerance; a counter increases at every evaluation at which t >= assessAfter and
all gating signals are converged, and is reset otherwise; the run stops when the counter
reaches stopConsecutive. Because all of these quantities are logged, the rule and variants
of its tolerance, reference scale, settling time and counter, as well as stricter readiness
conditions, can be replayed exactly without repeating the simulation. Variants in the window
length, the batch size or the search range cannot, because only the window difference at the
configured window is logged, and none is logged while a signal is not ready.

Options
  --gates a,b        gating signals (default: vol_Umag,vol_p)
  --tol x            plateau tolerance in per cent (default 5)
  --assessAfter t    settling time before a stop may be issued (default 15)
  --consecutive n    consecutive converged evaluations required (default 2)
  --interior         additionally require that the MSER minimum of every gating signal is
                     interior, i.e. not on the half-record search limit. This needs the
                     number of samples at each evaluation, which is recovered from the
                     'Time = ' lines of the log (one sample per time step; use
                     --sampleInterval if the function object sampled every n-th step).
  --trace file.csv   write the per-evaluation trace

Usage:  python3 replayGate.py log.pimpleFoam [options]
"""
import argparse
import bisect
import re
import sys

NUM = r"([-+0-9.eE]+)"
ROW = re.compile(
    r"statConv\[(\w+)\] t=" + NUM + " nStat=" + NUM + " transientEnd=" + NUM + " mean=" + NUM
    + " rms=" + NUM + " Neff=" + NUM + " CImean=" + NUM + r"% CIrms=(?:\[-)?" + NUM
    + r"%(?:,\+" + NUM + r"%\])? period=" + NUM + " win=" + NUM + " plateauRel=" + NUM
    + "% converged=(true|false)")


def read_log(path):
    """Return (step_times, evaluations); each evaluation maps signal -> dict."""
    step_times, evals, cur = [], [], {}
    for line in open(path, errors="replace"):
        if line.startswith("Time = "):
            if cur:
                evals.append(cur)
                cur = {}
            step_times.append(float(line.split()[2]))
            continue
        m = ROW.search(line)
        if m:
            g = m.groups()
            cur[g[0]] = dict(t=float(g[1]), nStat=int(float(g[2])), ttr=float(g[3]),
                             mean=float(g[4]), rms=float(g[5]), Neff=float(g[6]),
                             ready=float(g[11]) > 0, plateau=float(g[12]))
    if cur:
        evals.append(cur)
    return step_times, evals


def n_samples(step_times, t, sample_interval):
    """Number of samples taken up to the evaluation at time t."""
    k = bisect.bisect_left(step_times, t)
    k = min(range(max(k - 1, 0), min(k + 2, len(step_times))),
            key=lambda j: abs(step_times[j] - t))
    return (k + 1 + sample_interval - 1) // sample_interval


def on_limit(n_total, n_stat, batch=5):
    """True when the MSER minimum lies on the half-record search limit."""
    nb = n_total // batch
    return n_total - n_stat == (nb // 2 - 1) * batch


def replay(step_times, evals, gates, tol, assess_after, consecutive, interior,
           sample_interval=1, trace=None):
    count = 0
    rows = []
    for ev in evals:
        if not all(g in ev for g in gates):
            count = 0
            continue
        t = ev[gates[0]]["t"]
        n_total = n_samples(step_times, t, sample_interval) if interior else None
        ok = t >= assess_after
        flags = []
        for g in gates:
            s = ev[g]
            limit = interior and on_limit(n_total, s["nStat"])
            conv = s["ready"] and s["plateau"] < tol and not limit
            flags.append((s["ttr"], s["ready"], s["plateau"], limit, conv))
            ok = ok and conv
        count = count + 1 if ok else 0
        rows.append((t, flags, count))
        if count >= consecutive:
            break
    if trace:
        with open(trace, "w") as f:
            hdr = ["t"] + [f"{g}_{k}" for g in gates
                           for k in ("ttr", "ready", "plateau", "onLimit", "converged")]
            f.write(",".join(hdr + ["counter"]) + "\n")
            for t, flags, c in rows:
                vals = [f"{t:.6g}"]
                for (ttr, rd, pl, lim, cv) in flags:
                    vals += [f"{ttr:.6g}", str(int(rd)), f"{pl:.6g}" if rd else "nan",
                             str(int(bool(lim))), str(int(cv))]
                f.write(",".join(vals + [str(c)]) + "\n")
    return rows[-1][0] if rows and rows[-1][2] >= consecutive else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("log")
    ap.add_argument("--gates", default="vol_Umag,vol_p")
    ap.add_argument("--tol", type=float, default=5.0)
    ap.add_argument("--assessAfter", type=float, default=15.0)
    ap.add_argument("--consecutive", type=int, default=2)
    ap.add_argument("--interior", action="store_true")
    ap.add_argument("--sampleInterval", type=int, default=1)
    ap.add_argument("--trace")
    a = ap.parse_args()
    gates = [g.strip() for g in a.gates.split(",") if g.strip()]
    step_times, evals = read_log(a.log)
    t_stop = replay(step_times, evals, gates, a.tol, a.assessAfter, a.consecutive,
                    a.interior, a.sampleInterval, a.trace)
    last = next((ev[gates[0]]["t"] for ev in reversed(evals) if gates[0] in ev), None)
    if t_stop is None:
        print(f"no stop within the logged record (last evaluation t = {last})")
        sys.exit(1)
    print(f"stop at t = {t_stop}")


if __name__ == "__main__":
    main()
