#!/usr/bin/env python3
"""Flow-development + FFT/period analysis of a URANS probe signal, following the
data-analysis of Aydinbakar et al. (2021): reference time T = L/U with
L = (3 + 0.65*pi + 6.5) D; FFT amplitude vs frequency in 1/T; the lowest-frequency
local maximum sets the characteristic period, from which an adequate averaging window
(integer multiples of that period) is chosen.

Usage: fftPeriod.py <signal.dat: t value> [t_transient] [outprefix]
Writes <prefix>_spectrum.dat (f[1/T], amp), <prefix>_devmean.dat (window_end[T], mean, rms),
and prints the characteristic frequency/period and the running-mean convergence.
"""
import sys
import numpy as np
from scipy import signal as sig
from scipy.interpolate import interp1d

L = 3.0 + 0.65*np.pi + 6.5     # 11.542 D  (inlet 3D + bend arc 0.65*pi + outlet 6.5D)
T = L                          # T = L/U_b, U_b = 1 -> T = 11.542 convective times

fn = sys.argv[1]
t_tr = float(sys.argv[2]) if len(sys.argv) > 2 else None
pref = sys.argv[3] if len(sys.argv) > 3 else "sig"

d = np.loadtxt(fn)
t, x = d[:, 0], d[:, 1]

# initial-transient removal (MSER-style value passed in, else 1 flow-through)
if t_tr is None:
    t_tr = min(T, 0.4*(t[-1]-t[0]))
m = t >= t_tr
t, x = t[m], x[m]
print(f"# T = {T:.4f} convective times;  record {t[0]:.3f}..{t[-1]:.3f} "
      f"= {(t[-1]-t[0])/T:.2f} T after transient t_tr={t_tr:.3f}")

# uniform resample (signal has adaptive dt)
dt_u = np.median(np.diff(t))
tu = np.arange(t[0], t[-1], dt_u)
xu = interp1d(t, x, kind="linear")(tu)
N = len(tu)

# --- FFT amplitude spectrum (2021-style), Hann window, single-sided ---
xw = (xu - xu.mean()) * sig.windows.hann(N)
amp = 2.0/np.sum(sig.windows.hann(N)) * np.abs(np.fft.rfft(xw))
f_conv = np.fft.rfftfreq(N, dt_u)          # cycles per convective time
f_T = f_conv * T                            # cycles per T
np.savetxt(pref+"_spectrum.dat", np.column_stack([f_T[1:], amp[1:]]),
           header="f[1/T] amplitude", comments="")

# lowest-frequency local maximum (skip f<0.2/T to avoid the DC skirt)
lo = f_T > 0.2
idx = np.where(lo)[0]
loc = [i for i in idx[1:-1] if amp[i] > amp[i-1] and amp[i] > amp[i+1]]
if loc:
    i0 = loc[0]
    fpk, Tp = f_T[i0], 1.0/f_T[i0]
    print(f"# lowest local-max: f = {fpk:.3f} /T  ->  period = {Tp:.3f} T "
          f"= {Tp*T:.2f} convective times")
    print(f"# suggested averaging window (>= 2x period): {int(np.ceil(2*Tp))} T")
else:
    print("# no clear local maximum above 0.2/T (broadband) -> use integral-time / CI criterion")

# --- flow development: running mean & rms vs averaging-window length ---
rows = []
for we in np.arange(t_tr+0.5*T, t[-1]+1e-9, 0.25*T):
    mm = tu <= we
    if mm.sum() < 10:
        continue
    seg = xu[mm]
    rows.append([(we-t_tr)/T, seg.mean(), seg.std()])
rows = np.array(rows)
np.savetxt(pref+"_devmean.dat", rows, header="window[T] mean rms", comments="")
if len(rows):
    print(f"# running mean over last windows: {rows[-3:,1]}")
    print(f"# final mean={rows[-1,1]:.4f}  rms={rows[-1,2]:.4f}  "
          f"window={rows[-1,0]:.2f} T")
