#!/usr/bin/env python3
"""
Grid Convergence Index (GCI) for a three-grid study, following the procedure of Celik et al.,
J. Fluids Eng. 130 (2008) 078001, with the safety factor 1.25 and the representative grid
size h ~ N^(-1/3). For time-averaged quantities of unsteady simulations the three values
should be averaged over a common window, and the sampling uncertainty of each should be
small compared with the grid-to-grid differences.

Usage: gci.py N1 f1  N2 f2  N3 f3   (1 = finest, 3 = coarsest; N = cell count, f = metric)
"""
import sys, math

def celik(N1, f1, N2, f2, N3, f3, p=None, Fs=1.25):
    """Celik et al. (2008) three-grid procedure. 1 = finest. Returns (p, f_ext, GCI_fine,
    GCI_medium) with the GCIs as fractions; if p is given, it is used instead of the observed
    order."""
    # representative grid size h ~ N^(-1/3); refinement ratios
    r21 = (N1 / N2) ** (1.0 / 3.0)
    r32 = (N2 / N3) ** (1.0 / 3.0)
    e21 = f1 - f2
    e32 = f2 - f3
    if p is None:
        s = math.copysign(1.0, e32 / e21) if e21 != 0 else 1.0
        # solve apparent order p iteratively (Celik 2008 eq., handles non-constant r)
        p = 2.0
        for _ in range(500):
            q = math.log((r21 ** p - s) / (r32 ** p - s)) if (r21 != r32) else 0.0
            p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
            if abs(p_new - p) < 1e-10:
                p = p_new
                break
            p = p_new
    f_ext = (r21 ** p * f1 - f2) / (r21 ** p - 1.0)
    gci21 = Fs * abs((f1 - f2) / f1) / (r21 ** p - 1.0)      # fine-grid GCI
    gci32 = Fs * abs((f2 - f3) / f2) / (r32 ** p - 1.0)      # medium-grid GCI
    return p, f_ext, gci21, gci32


def main():
    a = [float(x) for x in sys.argv[1:7]]
    N1, f1, N2, f2, N3, f3 = a
    r21 = (N1 / N2) ** (1.0 / 3.0)
    r32 = (N2 / N3) ** (1.0 / 3.0)
    p, f_ext, gci21, gci32 = celik(N1, f1, N2, f2, N3, f3)
    p1, _, gci21_p1, _ = celik(N1, f1, N2, f2, N3, f3, p=1.0)
    print("=== 3-grid GCI (Roache / Celik 2008) ===")
    print(f"  N (fine,med,coarse) = {N1:.0f}, {N2:.0f}, {N3:.0f}")
    print(f"  f (fine,med,coarse) = {f1:.6g}, {f2:.6g}, {f3:.6g}")
    print(f"  r21={r21:.4f}  r32={r32:.4f}")
    print(f"  apparent order p    = {p:.3f}")
    print(f"  Richardson extrap f = {f_ext:.6g}")
    print(f"  GCI_fine (21)       = {gci21*100:.3f} %")
    print(f"  GCI_medium (32)     = {gci32*100:.3f} %")
    print(f"  asymptotic ratio    = {gci32/(r21**p*gci21):.3f}  (should be ~1)")
    print(f"  GCI_fine with p = 1 = {gci21_p1*100:.3f} %")


if __name__ == "__main__":
    main()
