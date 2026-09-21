#!/usr/bin/env python3
"""
Roache Grid Convergence Index (GCI) for a 3-grid study (ASME V&V 20 / Celik et al. 2008).
Valid ONLY for iteratively-converged, steady solutions (here: k-epsilon, which converges;
NOT SST/SSG which stall on this unsteady-separation flow -> those need URANS time-averaging).

Usage: gci.py N1 f1  N2 f2  N3 f3   (1=finest, 3=coarsest; N = cell count, f = metric)
"""
import sys, math

def main():
    a = [float(x) for x in sys.argv[1:7]]
    N1, f1, N2, f2, N3, f3 = a
    # representative grid size h ~ N^(-1/3); refinement ratios
    r21 = (N1 / N2) ** (1.0 / 3.0)
    r32 = (N2 / N3) ** (1.0 / 3.0)
    e21 = f1 - f2
    e32 = f2 - f3
    s = math.copysign(1.0, e32 / e21) if e21 != 0 else 1.0
    # solve apparent order p iteratively (Celik 2008 eq., handles non-constant r)
    p = 2.0
    for _ in range(100):
        q = math.log((r21 ** p - s) / (r32 ** p - s)) if (r21 != r32) else 0.0
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-8:
            p = p_new; break
        p = p_new
    f_ext = (r21 ** p * f1 - f2) / (r21 ** p - 1.0)
    ea21 = abs((f1 - f2) / f1)                      # approx relative error
    eext21 = abs((f_ext - f1) / f_ext)              # extrapolated relative error
    gci21 = 1.25 * ea21 / (r21 ** p - 1.0)          # fine-grid GCI (Fs=1.25)
    gci32 = 1.25 * abs((f2 - f3) / f2) / (r32 ** p - 1.0)
    print("=== 3-grid GCI (Roache / Celik 2008) ===")
    print(f"  N (fine,med,coarse) = {N1:.0f}, {N2:.0f}, {N3:.0f}")
    print(f"  f (fine,med,coarse) = {f1:.6g}, {f2:.6g}, {f3:.6g}")
    print(f"  r21={r21:.4f}  r32={r32:.4f}")
    print(f"  apparent order p    = {p:.3f}")
    print(f"  Richardson extrap f = {f_ext:.6g}")
    print(f"  GCI_fine (21)       = {gci21*100:.3f} %")
    print(f"  GCI_medium (32)     = {gci32*100:.3f} %")
    print(f"  asymptotic ratio    = {gci32/(r21**p*gci21):.3f}  (should be ~1)")

if __name__ == "__main__":
    main()
