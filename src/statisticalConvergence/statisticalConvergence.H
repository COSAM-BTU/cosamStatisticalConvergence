/*---------------------------------------------------------------------------*\
  statisticalConvergence functionObject (COSAM)
  In-situ statistical-convergence / stationarity assessment of scalar signals
  (probe-cell values of a field component) for turbulent CFD.

  Per signal, at write time: MSER-m initial-transient detection, mean, rms,
  autocorrelation integral time scale (T0), effective sample size Neff,
  autocorrelation-corrected 95% CIs for mean AND rms, dominant frequency (DFT).

  Convergence decision (autoStop): plateau / range-insensitivity of the running
  time-mean. The two most-recent, consecutive, disjoint averaging windows of
  length compareWindow are compared; convergence is declared when their mean
  difference falls below tolPlateau (relative to refScale, a fixed velocity
  scale, so it is well defined even for near-zero-mean signals), sustained over
  stopConsecutive checks. This mirrors the averaging-range insensitivity
  criterion of Aydinbakar et al., Comput. Mech. (2021) -- "(20T,30T) vs
  (30T,40T) not significant". The 95% CIs are reported as honest uncertainty
  but are NOT the stop gate.

  Reference method: Oliver et al. Phys.Fluids 26:035101 (2014); White 1997 /
  Wang & Glynn ACM TOMACS (2016). Validated against Python reference statConv.py.
\*---------------------------------------------------------------------------*/
#ifndef statisticalConvergence_H
#define statisticalConvergence_H

#include "fvMeshFunctionObject.H"
#include "DynamicList.H"
#include "point.H"
#include "statConvMath.H"

namespace Foam
{
namespace functionObjects
{

class statisticalConvergence
:
    public fvMeshFunctionObject
{
    // Per-signal data
    struct signalData
    {
        word   name;
        word   mode;                    // "point" (cell probe) or "volAverage"
        bool   gate;                    // participate in the autoStop decision?
        point  location;
        word   field;
        label  component;               // vector comp; <0 -> magnitude
        label  cell;                    // owning cell on this rank (-1 if none)
        DynamicList<scalar> series;
    };

    // Per-signal analysis outcome (statConvMath result + plateau test)
    struct sigResult
    {
        statConvMath::Result r;
        scalar plateauRel;      // |mean(recent win) - mean(previous win)| / max(|mean|,refScale)
        scalar windowSpan;      // time spanned by the two comparison windows (2*compareWindow)
        bool   ready;           // enough stationary data for two comparison windows
        bool   converged;       // plateau below tolPlateau AND ready
    };

    // Private data
    DynamicList<signalData> signals_;
    DynamicList<scalar>     times_;      // sample times (shared)
    label   sampleInterval_;
    label   stepCounter_;
    scalar  tolMean_;
    scalar  tolRms_;
    label   minSamples_;

    // Auto-stop controls (plateau / range-insensitivity)
    bool    autoStop_;
    scalar  compareWindow_;    // length (time units) of EACH of the two consecutive windows
    scalar  tolPlateau_;       // two-window mean-agreement tol (fraction of refScale/|mean|)
    scalar  refScale_;         // reference velocity scale U_b (robust for zero-mean)
    scalar  assessAfter_;      // only assess convergence / allow autoStop after this time
    label   stopConsecutive_;  // sustained consecutive converged checks
    label   consecCount_;

    // Private members
    scalar    sampleValue(const signalData&) const;
    sigResult analyse(const signalData&) const;   // statConvMath + plateau test

public:

    TypeName("statisticalConvergence");

    statisticalConvergence
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~statisticalConvergence() = default;

    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
