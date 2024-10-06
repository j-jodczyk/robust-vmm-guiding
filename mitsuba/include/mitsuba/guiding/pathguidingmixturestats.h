/*
    This file is part of the implementation of the SIGGRAPH 2020 paper
    "Robust Fitting of Parallax-Aware Mixtures for Path Guiding".
    The implementation extends Mitsuba, a physically based rendering system.

    Copyright (c) 2020 Lukas Ruppert, Sebastian Herholz.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PATHGUIDINGMIXTURESTATS_H
#define PATHGUIDINGMIXTURESTATS_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/serialization.h>
#include <mitsuba/core/spectrum.h>

#include <mitsuba/guiding/incrementalcovariance2d.h>
#include <mitsuba/guiding/incrementaldistance.h>
#include <mitsuba/guiding/incrementalpearsonchisquared.h>

#include <mitsuba/guiding/guiding.h>
#include <mitsuba/guiding/SampleStats.h>

GUIDING_NAMESPACE_BEGIN

template<typename TDistribution>
struct PathGuidingMixtureStats
{
    typedef TDistribution DistributionType;

    Guiding::SampleStats lastFitStats;
    IncrementalDistance<TDistribution>          incrementalDistance;
    IncrementalPearsonChiSquared<TDistribution> incrementalPearsonChiSquared;
    IncrementalCovariance2D<TDistribution>      incrementalCovariance2D;
    Spectrum flux {0.0f};
    size_t samplesSinceLastMerge {0};

    PathGuidingMixtureStats() = default;

    explicit PathGuidingMixtureStats(Stream* stream)
        : lastFitStats{stream}, incrementalDistance{stream}, incrementalPearsonChiSquared{stream}, incrementalCovariance2D{stream}, flux{stream}, samplesSinceLastMerge{stream->readSize()}
    {

    }

    void serialize(Stream* stream) const
    {
        lastFitStats.serialize(stream);
        incrementalDistance.serialize(stream);
        incrementalPearsonChiSquared.serialize(stream);
        incrementalCovariance2D.serialize(stream);
        flux.serialize(stream);
        stream->writeSize(samplesSinceLastMerge);
    }

    ///decay the weights of the contained information
    void operator *=(const float factor)
    {
        incrementalDistance          *= factor;
        incrementalPearsonChiSquared *= factor;
        incrementalCovariance2D      *= factor;
    }
    ///clear all information
    void clear()
    {
        *this = std::move(PathGuidingMixtureStats<TDistribution>{});
    }

    std::string toString() const
    {
        std::ostringstream oss;

        oss << "PathGuidingMixtureStats [\n"
            << "  lastFitStats = " << lastFitStats.toString() << '\n'
            << "  incrementalDistance = " << incrementalDistance.toString() << '\n'
            << "  incrementalPearsonChiSquared = " << incrementalPearsonChiSquared.toString() << '\n'
            << "  incrementalCovariance2D = " << incrementalCovariance2D.toString() << '\n'
            << "  flux = " << flux.toString() << '\n'
            << "  samplesSinceLastMerge = " << samplesSinceLastMerge << '\n'
            << ']';

        return oss.str();
    }
};

GUIDING_NAMESPACE_END

#endif // PATHGUIDINGMIXTURESTATS_H
