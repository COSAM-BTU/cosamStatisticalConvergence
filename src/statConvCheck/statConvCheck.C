/*---------------------------------------------------------------------------*\
  statConvCheck -- standalone verification app for the statisticalConvergence
  math. Reads a 2-column text file (time value), runs the SHARED statConvMath::
  analyse (identical to the functionObject), and prints the result. Compared
  bit-for-bit against the Python reference scripts/statConv.py on the same file.

  usage: statConvCheck <file> [tolMean] [tolRms]
\*---------------------------------------------------------------------------*/
#include "statConvMath.H"
#include "DynamicList.H"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace Foam;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: statConvCheck <file> [tolMean] [tolRms]\n";
        return 1;
    }
    const scalar tolMean = (argc > 2 ? std::atof(argv[2]) : 0.01);
    const scalar tolRms  = (argc > 3 ? std::atof(argv[3]) : 0.05);

    DynamicList<scalar> t, v;
    std::ifstream f(argv[1]);
    std::string line;
    while (std::getline(f, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double a, b;
        if (ss >> a >> b) { t.append(a); v.append(b); }
    }

    const scalar dt =
        (t.size() > 1 ? (t[t.size()-1] - t[0])/(t.size()-1) : scalar(1));

    const statConvMath::Result r = statConvMath::analyse(v, dt, tolMean, tolRms);

    std::cout.precision(9);
    std::cout
        << "nTotal "          << r.nTotal
        << " transientEndIdx "<< r.transientEnd
        << " nStat "          << r.nStat
        << " mean "           << r.mean
        << " rms "            << r.rms
        << " T0 "             << r.T0
        << " integralTime "   << r.integralTime
        << " Neff "           << r.Neff
        << " ciMeanRel "      << r.ciMeanRel
        << " ciRmsRel "       << r.ciRmsRel
        << " ciRmsRelLo "     << r.ciRmsRelLo
        << " ciRmsRelHi "     << r.ciRmsRelHi
        << " fdom "           << r.dominantFreq
        << " samplingAdequate " << (r.samplingAdequate ? 1 : 0)
        << "\n";
    return 0;
}
