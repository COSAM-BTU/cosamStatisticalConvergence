#!/usr/bin/env python3
"""Generate the synthetic AR(1) test signal used to verify the statistical kernel.

    x_i = phi x_{i-1} + eps_i,   phi = 0.9

for which the statistical inefficiency has the closed form T0 = (1+phi)/(1-phi) = 19.
A linear ramp is added over the first 2000 samples so that the transient-detection
stage is exercised as well. The stream is seeded, so the file is reproducible.

Usage:  ./makeSignal.py [--out arSignal.dat] [--n 40000] [--phi 0.9] [--dt 0.01]
"""
import argparse

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--out", default="arSignal.dat")
p.add_argument("--n", type=int, default=40000)
p.add_argument("--phi", type=float, default=0.9)
p.add_argument("--dt", type=float, default=0.01)
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()

rng = np.random.default_rng(a.seed)
x = np.zeros(a.n)
for i in range(1, a.n):
    x[i] = a.phi * x[i - 1] + rng.standard_normal()
x += 5.0                                             # non-zero mean
x[:2000] += np.linspace(10, 0, 2000)                 # initial transient
t = np.arange(a.n) * a.dt

with open(a.out, "w") as f:
    f.write("# synthetic AR(1), phi=%g, theoretical T0=%g\n" % (a.phi, (1 + a.phi) / (1 - a.phi)))
    for ti, xi in zip(t, x):
        f.write("%.6f %.10f\n" % (ti, xi))
print("wrote %s (%d samples)" % (a.out, a.n))
