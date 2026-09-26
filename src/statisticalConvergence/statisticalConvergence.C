/*---------------------------------------------------------------------------*\
  statisticalConvergence functionObject (COSAM) -- implementation.
  Uses the kernel of statConvMath.H, verified against the Python implementation
  tools/statConv.py, and adds an autoStop plateau / range-insensitivity rule
  (Aydinbakar et al. 2021).
\*---------------------------------------------------------------------------*/
#include "statisticalConvergence.H"
#include "statConvMath.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(statisticalConvergence, 0);
    addToRunTimeSelectionTable(functionObject, statisticalConvergence, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::statisticalConvergence::statisticalConvergence
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    signals_(),
    times_(),
    sampleInterval_(1),
    stepCounter_(0),
    tolMean_(0.01),
    tolRms_(0.05),
    minSamples_(100),
    autoStop_(false),
    compareWindow_(10.0),
    tolPlateau_(0.05),
    refScale_(1.0),
    assessAfter_(0.0),
    stopConsecutive_(2),
    consecCount_(0)
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::functionObjects::statisticalConvergence::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    sampleInterval_ = dict.getOrDefault<label>("sampleInterval", 1);
    tolMean_        = dict.getOrDefault<scalar>("tolMean", 0.01);
    tolRms_         = dict.getOrDefault<scalar>("tolRms", 0.05);
    minSamples_     = dict.getOrDefault<label>("minSamples", 100);

    // Auto-stop (plateau / range-insensitivity) controls
    autoStop_        = dict.getOrDefault<bool>("autoStop", false);
    compareWindow_   = dict.getOrDefault<scalar>("compareWindow", 10.0);
    tolPlateau_      = dict.getOrDefault<scalar>("tolPlateau", 0.05);
    refScale_        = dict.getOrDefault<scalar>("refScale", 1.0);
    assessAfter_     = dict.getOrDefault<scalar>("assessAfter", 0);
    stopConsecutive_ = dict.getOrDefault<label>("stopConsecutive", 2);
    consecCount_     = 0;

    signals_.clear();
    times_.clear();         // histories restart together when the dictionary is re-read
    const dictionary& sigs = dict.subDict("signals");
    for (const entry& e : sigs)
    {
        if (!e.isDict()) continue;
        const dictionary& sd = e.dict();
        signalData s;
        s.name      = e.keyword();
        s.mode      = sd.getOrDefault<word>("mode", "point");
        s.field     = sd.get<word>("field");
        s.component = sd.getOrDefault<label>("component", 0);
        s.gate      = sd.getOrDefault<bool>("autoStopGate", true);

        if (s.mode == "volAverage")
        {
            // spatial (volume-weighted) average over the whole domain: a
            // low-variance monitor signal -> clean plateau / reliable autoStop
            // even in high-turbulence-intensity flows (point signals are noisy).
            s.location = point::zero;
            s.cell     = -1;
        }
        else
        {
            s.location = sd.get<point>("location");
            s.cell     = mesh_.findCell(s.location);
            const label owners =
                returnReduce(label(s.cell >= 0 ? 1 : 0), sumOp<label>());
            if (owners == 0)
            {
                WarningInFunction
                    << "signal " << s.name << ": location " << s.location
                    << " not found in any cell." << endl;
            }
        }
        signals_.append(s);
    }

    Info<< type() << " " << name() << ": monitoring " << signals_.size()
        << " signal(s), sampleInterval=" << sampleInterval_
        << ", autoStop=" << (autoStop_ ? "true" : "false")
        << ", compareWindow=" << compareWindow_
        << ", tolPlateau=" << tolPlateau_ << " (refScale=" << refScale_
        << "), assessAfter=" << assessAfter_
        << ", stopConsecutive=" << stopConsecutive_ << endl;

    return true;
}


Foam::scalar Foam::functionObjects::statisticalConvergence::sampleValue
(
    const signalData& s
) const
{
    // ---- spatial (volume-weighted) average over the whole domain ----
    if (s.mode == "volAverage")
    {
        const scalarField& Vc = mesh_.V();
        scalar sumFV = 0, sumV = 0;
        if (mesh_.foundObject<volScalarField>(s.field))
        {
            const volScalarField& f = mesh_.lookupObject<volScalarField>(s.field);
            forAll(f, c) { sumFV += f[c]*Vc[c]; sumV += Vc[c]; }
        }
        else if (mesh_.foundObject<volVectorField>(s.field))
        {
            const volVectorField& f = mesh_.lookupObject<volVectorField>(s.field);
            if (s.component < 0)                 // magnitude
            {
                forAll(f, c) { sumFV += mag(f[c])*Vc[c]; sumV += Vc[c]; }
            }
            else
            {
                forAll(f, c) { sumFV += f[c][s.component]*Vc[c]; sumV += Vc[c]; }
            }
        }
        sumFV = returnReduce(sumFV, sumOp<scalar>());
        sumV  = returnReduce(sumV,  sumOp<scalar>());
        return (sumV > VSMALL ? sumFV/sumV : scalar(0));
    }

    // ---- point (single-cell probe) ----
    scalar local = 0;
    if (s.cell >= 0)
    {
        if (mesh_.foundObject<volScalarField>(s.field))
        {
            local = mesh_.lookupObject<volScalarField>(s.field)[s.cell];
        }
        else if (mesh_.foundObject<volVectorField>(s.field))
        {
            local =
                mesh_.lookupObject<volVectorField>(s.field)[s.cell][s.component];
        }
    }
    // exactly one rank owns the cell -> sum recovers the value on all ranks
    return returnReduce(local, sumOp<scalar>());
}


bool Foam::functionObjects::statisticalConvergence::execute()
{
    if ((stepCounter_++ % sampleInterval_) != 0)
    {
        return true;
    }
    times_.append(mesh_.time().value());
    forAll(signals_, i)
    {
        signals_[i].series.append(sampleValue(signals_[i]));
    }
    return true;
}


Foam::functionObjects::statisticalConvergence::sigResult
Foam::functionObjects::statisticalConvergence::analyse
(
    const signalData& sig
) const
{
    sigResult out;
    out.plateauRel = GREAT;
    out.windowSpan = 0;
    out.ready      = false;
    out.converged  = false;

    // default-initialise the embedded statConvMath::Result
    statConvMath::Result& r = out.r;
    r.nTotal = sig.series.size();
    r.transientEnd = 0; r.nStat = 0;
    r.mean = r.rms = r.T0 = r.integralTime = r.Neff = 0;
    r.ciMeanRel = r.ciRmsRel = r.dominantFreq = 0;
    r.stationary = false; r.samplingAdequate = false;

    if (sig.series.size() < minSamples_) return out;

    scalar dt = 1;
    if (times_.size() > 1)
    {
        dt = (times_[times_.size()-1] - times_[0])/(times_.size()-1);
    }

    // shared math (identical to standalone statConvCheck + validated vs statConv.py)
    r = statConvMath::analyse(sig.series, dt, tolMean_, tolRms_);
    if (r.nStat < 20) return out;

    // ---- plateau / range-insensitivity: compare the mean of the two most-recent
    //      consecutive disjoint windows, each of length compareWindow (time).
    //      Mirrors Aydinbakar 2021 "(20T,30T) vs (30T,40T) not significant". ----
    const label n   = r.nStat;
    const label off = r.transientEnd;
    label Lwin = label(compareWindow_/max(dt, SMALL) + 0.5);
    if (Lwin < 1) Lwin = 1;

    if (n >= 2*Lwin)
    {
        const label base = off + n;            // one past the last stationary sample
        scalar mRecent = 0, mPrev = 0;
        for (label i = base - Lwin;   i < base;        ++i) mRecent += sig.series[i];
        for (label i = base - 2*Lwin; i < base - Lwin; ++i) mPrev   += sig.series[i];
        mRecent /= Lwin;
        mPrev   /= Lwin;

        const scalar denom = max(mag(r.mean), refScale_);
        out.plateauRel = mag(mRecent - mPrev)/max(denom, SMALL);
        out.windowSpan = 2*Lwin*dt;
        out.ready      = true;
        out.converged  = (out.plateauRel < tolPlateau_);
    }

    const scalar tEnd =
        (r.transientEnd < times_.size() ? times_[r.transientEnd] : scalar(0));

    Info<< "  statConv[" << sig.name << "] t=" << mesh_.time().value()
        << " nStat=" << r.nStat << " transientEnd=" << tEnd
        << " mean=" << r.mean << " rms=" << r.rms
        << " Neff=" << r.Neff
        << " CImean=" << r.ciMeanRel*100 << "%"
        << " CIrms=[-" << r.ciRmsRelLo*100 << "%,+" << r.ciRmsRelHi*100 << "%]"
        << " period=" << (r.dominantFreq > 0 ? 1.0/r.dominantFreq : scalar(0))
        << " win=" << out.windowSpan
        << " plateauRel=" << out.plateauRel*100 << "%"
        << " converged=" << (out.converged ? "true" : "false") << endl;

    return out;
}


bool Foam::functionObjects::statisticalConvergence::write()
{
    Info<< type() << " " << name() << " write:" << endl;

    label nSig = 0, nConv = 0;
    forAll(signals_, i)
    {
        const sigResult s = analyse(signals_[i]);   // report every signal
        if (!signals_[i].gate) continue;            // but gate only on flagged
        ++nSig;
        if (s.converged) ++nConv;
    }

    // convergence is only assessed (and autoStop allowed) after assessAfter_,
    // so a warm-start settling phase cannot trigger a premature stop.
    const bool assessing = (mesh_.time().value() >= assessAfter_);
    const bool allConverged = assessing && (nSig > 0 && nConv == nSig);
    if (allConverged) ++consecCount_; else consecCount_ = 0;

    if (autoStop_)
    {
        if (!assessing)
        {
            Info<< "  statConv: settling phase (t < assessAfter=" << assessAfter_
                << "), convergence assessment not yet active." << endl;
        }
        else if (allConverged && consecCount_ >= stopConsecutive_)
        {
            Info<< "  statConv: PLATEAU CONVERGED -- all " << nSig
                << " signal(s), " << consecCount_ << " consecutive checks"
                << " (two consecutive " << compareWindow_ << "-window means agree"
                << " < " << tolPlateau_*100 << "% of max(|mean|, refScale))."
                << " autoStop -> writing fields and ending run at t="
                << mesh_.time().value() << endl;

            const_cast<Foam::Time&>(mesh_.time()).stopAt(Foam::Time::saWriteNow);
        }
        else
        {
            Info<< "  statConv: not converged yet (converged signals "
                << nConv << "/" << nSig << ", consecutive "
                << consecCount_ << "/" << stopConsecutive_ << ")." << endl;
        }
    }

    return true;
}

// ************************************************************************* //
