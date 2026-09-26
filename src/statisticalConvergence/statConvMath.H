/*---------------------------------------------------------------------------*\
  statConvMath.H -- shared, header-only statistical-convergence math used by BOTH
  the statisticalConvergence functionObject and the standalone statConvCheck
  verification app. Single source of truth -> C++ FO and the checker are identical,
  and both are verified against the Python implementation tools/statConv.py.
\*---------------------------------------------------------------------------*/
#ifndef statConvMath_H
#define statConvMath_H

#include "scalarField.H"
#include "SubList.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace statConvMath
{

struct Result
{
    label  nTotal, transientEnd, nStat;
    scalar mean, rms, T0, integralTime, Neff;
    scalar ciMeanRel, ciRmsRel, ciRmsRelLo, ciRmsRelHi, dominantFreq;
    bool   stationary, samplingAdequate;
};


// ---------------------------------------------------------------------------
// Small-sample quantiles. At the effective sample sizes reached in turbulent
// records (N_eff = O(5)) the normal quantile understates the interval on the
// mean by ~50% and the symmetric large-sample form misrepresents the interval
// on the rms, which is strongly asymmetric. Both are therefore taken from the
// Student-t and chi-squared distributions, with non-integer degrees of freedom.
// ---------------------------------------------------------------------------

// standard normal quantile (Acklam's rational approximation, |eps| < 1.15e-9)
inline scalar normalQuantile(const scalar p)
{
    static const scalar a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285104469687e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277459239e+00};
    static const scalar b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const scalar c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                 4.374664141464968e+00,  2.938163982698783e+00};
    static const scalar d[4] = { 7.784695709041462e-03,  3.224671290700398e-01,
                                 2.445134137142996e+00,  3.754408661907416e+00};
    const scalar pl = 0.02425;
    if (p <= 0 || p >= 1) return 0;
    if (p < pl)
    {
        const scalar q = Foam::sqrt(-2*Foam::log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])
             / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    if (p > 1 - pl)
    {
        const scalar q = Foam::sqrt(-2*Foam::log(1-p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])
              / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    const scalar q = p - 0.5, r = q*q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q
         / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
}

// incomplete beta continued fraction (Lentz), for the Student-t cdf
inline scalar betaCF(const scalar a, const scalar b, const scalar x)
{
    const scalar tiny = 1e-30;
    scalar qab = a + b, qap = a + 1, qam = a - 1;
    scalar c = 1, d = 1 - qab*x/qap;
    if (mag(d) < tiny) d = tiny;
    d = 1/d;
    scalar h = d;
    for (label m = 1; m <= 300; ++m)
    {
        const label m2 = 2*m;
        scalar aa = m*(b - m)*x/((qam + m2)*(a + m2));
        d = 1 + aa*d; if (mag(d) < tiny) d = tiny;
        c = 1 + aa/c; if (mag(c) < tiny) c = tiny;
        d = 1/d; h *= d*c;
        aa = -(a + m)*(qab + m)*x/((a + m2)*(qap + m2));
        d = 1 + aa*d; if (mag(d) < tiny) d = tiny;
        c = 1 + aa/c; if (mag(c) < tiny) c = tiny;
        d = 1/d;
        const scalar del = d*c;
        h *= del;
        if (mag(del - 1) < 1e-14) break;
    }
    return h;
}

// regularized incomplete beta I_x(a,b)
inline scalar betaI(const scalar a, const scalar b, const scalar x)
{
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const scalar bt = Foam::exp
    (
        Foam::lgamma(a + b) - Foam::lgamma(a) - Foam::lgamma(b)
      + a*Foam::log(x) + b*Foam::log(1 - x)
    );
    return (x < (a + 1)/(a + b + 2))
         ?     bt*betaCF(a, b, x)/a
         : 1 - bt*betaCF(b, a, 1 - x)/b;
}

// Student-t cdf and quantile (bisection); valid for non-integer nu > 0
inline scalar tCDF(const scalar t, const scalar nu)
{
    const scalar x = nu/(nu + t*t);
    const scalar ib = 0.5*betaI(0.5*nu, 0.5, x);
    return (t > 0 ? 1 - ib : ib);
}

inline scalar tQuantile(const scalar p, const scalar nu)
{
    if (nu <= 0 || p <= 0 || p >= 1) return 0;
    scalar lo = 0, hi = 10;
    while (tCDF(hi, nu) < p && hi < 1e6) hi *= 2;
    for (label i = 0; i < 200; ++i)
    {
        const scalar mid = 0.5*(lo + hi);
        (tCDF(mid, nu) < p ? lo : hi) = mid;
        if (hi - lo < 1e-10*max(scalar(1), hi)) break;
    }
    return 0.5*(lo + hi);
}

// regularized lower incomplete gamma P(a,x) (Numerical Recipes series / CF)
inline scalar gammaP(const scalar a, const scalar x)
{
    if (x <= 0 || a <= 0) return 0;
    const scalar gln = Foam::lgamma(a);
    if (x < a + 1)
    {
        scalar ap = a, del = 1.0/a, sum = del;
        for (label i = 0; i < 300; ++i)
        {
            ap += 1; del *= x/ap; sum += del;
            if (mag(del) < mag(sum)*1e-14) break;
        }
        return sum*Foam::exp(-x + a*Foam::log(x) - gln);
    }
    scalar bb = x + 1 - a, cc = 1e30, dd2 = 1.0/bb, h = dd2;
    for (label i = 1; i <= 300; ++i)
    {
        const scalar an = -scalar(i)*(scalar(i) - a);
        bb += 2;
        dd2 = an*dd2 + bb; if (mag(dd2) < 1e-30) dd2 = 1e-30;
        cc = bb + an/cc;   if (mag(cc)  < 1e-30) cc  = 1e-30;
        dd2 = 1.0/dd2;
        const scalar del = dd2*cc;
        h *= del;
        if (mag(del - 1.0) < 1e-14) break;
    }
    return 1.0 - Foam::exp(-x + a*Foam::log(x) - gln)*h;
}

// chi-squared quantile by bisection on P(nu/2, x/2) = p
inline scalar chi2Quantile(const scalar p, const scalar nu)
{
    if (nu <= 0 || p <= 0 || p >= 1) return 0;
    scalar lo = 1e-8, hi = max(scalar(10), 4*nu);
    while (gammaP(0.5*nu, 0.5*hi) < p && hi < 1e8) hi *= 2;
    for (label i = 0; i < 200; ++i)
    {
        const scalar mid = 0.5*(lo + hi);
        (gammaP(0.5*nu, 0.5*mid) < p ? lo : hi) = mid;
        if (hi - lo < 1e-10*max(scalar(1), hi)) break;
    }
    return 0.5*(lo + hi);
}

// MSER-m (m=5) initial-transient truncation index (in original samples)
inline label mserTruncation(const UList<scalar>& x)
{
    const label m = 5;
    const label nb = x.size()/m;
    if (nb < 8) return 0;
    scalarField b(nb, Zero);
    for (label j = 0; j < nb; ++j)
    {
        scalar s = 0;
        for (label k = 0; k < m; ++k) s += x[j*m + k];
        b[j] = s/m;
    }
    label bestD = 0; scalar bestVal = GREAT;
    for (label d = 0; d < nb/2; ++d)
    {
        const label z = nb - d;
        scalar mean = 0;
        for (label j = d; j < nb; ++j) mean += b[j];
        mean /= z;
        scalar var = 0;
        for (label j = d; j < nb; ++j) var += sqr(b[j] - mean);
        var /= z;
        const scalar val = var/z;
        if (val < bestVal) { bestVal = val; bestD = d; }
    }
    return bestD*m;
}

// Full analysis on a scalar series with mean sampling interval dt.
inline Result analyse
(
    const UList<scalar>& all,
    const scalar dt,
    const scalar tolMean,
    const scalar tolRms
)
{
    Result r;
    r.nTotal = all.size();
    r.transientEnd = mserTruncation(all);
    const label n = r.nTotal - r.transientEnd;
    r.nStat = n;
    r.mean = r.rms = r.T0 = r.integralTime = r.Neff = 0;
    r.ciMeanRel = r.ciRmsRel = r.ciRmsRelLo = r.ciRmsRelHi = r.dominantFreq = 0;
    r.stationary = false; r.samplingAdequate = false;
    if (n < 20) return r;

    const SubList<scalar> x(all, n, r.transientEnd);

    scalar mean = 0;
    for (const scalar v : x) mean += v;
    mean /= n;
    scalar var = 0;
    for (const scalar v : x) var += sqr(v - mean);
    var /= (n - 1);
    const scalar rms = Foam::sqrt(var);

    // integral time scale: sum the ACF up to the lag before its first non-positive
    // value (first-zero-crossing truncation, O'Neill et al. 2004; Mora & Obligado 2020)
    scalar T0 = 1.0;
    if (var > VSMALL)
    {
        const label maxLag = min(n-2, max(50, n/4));
        scalar s = 0;
        for (label k = 1; k <= maxLag; ++k)
        {
            scalar c = 0;
            for (label i = 0; i < n-k; ++i) c += (x[i]-mean)*(x[i+k]-mean);
            c /= (scalar(n)*var);
            if (c <= 0) break;
            s += c;
        }
        T0 = 1.0 + 2.0*s;
    }
    const scalar Neff = (T0 > 0 ? n/T0 : scalar(n));

    // small-sample intervals: Student-t for the mean, chi-squared for the rms
    const scalar nu = max(Neff - 1.0, scalar(0.5));
    const scalar tq = tQuantile(0.975, nu);
    const scalar hwMean = tq*rms/Foam::sqrt(max(Neff, SMALL));
    const scalar chiHi  = chi2Quantile(0.975, nu);   // -> lower bound on sigma
    const scalar chiLo  = chi2Quantile(0.025, nu);   // -> upper bound on sigma
    const scalar sigLo  = (chiHi > 0 ? rms*Foam::sqrt(nu/chiHi) : 0);
    const scalar sigHi  = (chiLo > 0 ? rms*Foam::sqrt(nu/chiLo) : GREAT);
    const scalar hwRms  = 0.5*(sigHi - sigLo);       // nominal half-width, reported for continuity

    // dominant frequency via Welch PSD (Hann window, 50% overlap, averaged
    // periodograms) -- matches scipy.signal.welch in the Python reference. Segment
    // detrended per Welch to suppress slow drift. Resolution = fs/L (first bin).
    const scalar twoPi = 2.0*constant::mathematical::pi;
    const scalar fs = 1.0/dt;
    scalar fdom = 0;
    // segment length scales with the record so that low-frequency content is
    // resolvable; a fixed cap would report the first bin for every long record
    const label L = min(max(label(256), n/8), n/2);
    if (L >= 8)
    {
        const label step = max(label(1), L/2);
        const label nfreq = L/2;
        scalarField w(L);
        for (label i = 0; i < L; ++i)
            w[i] = 0.5*(1.0 - Foam::cos(twoPi*i/(L-1)));
        scalarField psd(nfreq, Zero);
        label nSeg = 0;
        for (label s0 = 0; s0 + L <= n; s0 += step)
        {
            scalar segMean = 0;
            for (label i = 0; i < L; ++i) segMean += x[s0+i];
            segMean /= L;
            for (label k = 1; k < nfreq; ++k)
            {
                scalar re = 0, im = 0;
                for (label i = 0; i < L; ++i)
                {
                    const scalar xi = (x[s0+i]-segMean)*w[i];
                    const scalar ph = twoPi*k*i/L;
                    re += xi*Foam::cos(ph);
                    im -= xi*Foam::sin(ph);
                }
                psd[k] += re*re + im*im;
            }
            ++nSeg;
        }
        if (nSeg > 0)
        {
            scalar pmax = 0;
            for (label k = 1; k < nfreq; ++k)
                if (psd[k] > pmax) { pmax = psd[k]; fdom = fs*scalar(k)/L; }
        }
    }

    r.mean = mean; r.rms = rms; r.T0 = T0;
    r.integralTime = T0*dt/2.0; r.Neff = Neff;
    r.ciMeanRel = (mag(mean) > VSMALL ? hwMean/mag(mean) : GREAT);
    r.ciRmsRel  = (rms > VSMALL ? hwRms/rms : GREAT);
    r.ciRmsRelLo = (rms > VSMALL ? (rms - sigLo)/rms : GREAT);
    r.ciRmsRelHi = (rms > VSMALL ? (sigHi - rms)/rms : GREAT);
    r.dominantFreq = fdom;
    r.stationary = (r.transientEnd < r.nTotal - 20);
    r.samplingAdequate = (r.ciMeanRel < tolMean) && (r.ciRmsRel < tolRms);
    return r;
}

} // End namespace statConvMath
} // End namespace Foam

#endif

// ************************************************************************* //
