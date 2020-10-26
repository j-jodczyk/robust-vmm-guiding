/*
    This file is part of the implementation of the SIGGRAPH 2020 paper
    "Robust Fitting of Parallax-Aware Mixtures for Path Guiding",
    as well as the updated implementation of the ACM TOG 2019 paper
    "Volume Path Guiding Based on Zero-Variance Random Walk Theory".
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

#pragma once

#include "guiding.h"

#include "SampleStats.h"

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/core/stream.h>

GUIDING_NAMESPACE_BEGIN

template<typename TDistribution, typename TStatistics>
struct alignas(alignof(TDistribution)) GuidingRegion
{
    typedef TDistribution DistributionType;
    typedef TStatistics   StatisticsType;
    
    GuidingRegion() = default;

    explicit GuidingRegion(Stream *stream)
    {
        valid = stream->readBool();
        spatialSplitFlag = stream->readBool();
        size_t dataStreamSize = stream->readSize();
        char dataStreamString[dataStreamSize];
        stream->read(dataStreamString, dataStreamSize);
        std::istringstream dataSerializationStream{std::move(std::string{dataStreamString, dataStreamSize})};
        distribution.deserialize(dataSerializationStream);
        statistics = StatisticsType(stream);

        statsSinceLastSplit = SampleStats(stream);
        statsUntilLastSplit = SampleStats(stream);
        dStart = stream->readUInt();
        nSamples = stream->readUInt();
        bound = AABB(stream);
    }

    void serialize(Stream* stream) const
    {
        stream->writeBool(valid);
        stream->writeBool(spatialSplitFlag);
        std::ostringstream dataSerializationStream;
        distribution.serialize(dataSerializationStream);
        const std::string dataStreamString = dataSerializationStream.str();
        const size_t dataStreamSize = dataStreamString.size();
        stream->writeSize(dataStreamSize);
        stream->write(dataStreamString.data(), dataStreamSize);
        statistics.serialize(stream);

        statsSinceLastSplit.serialize(stream);
        statsUntilLastSplit.serialize(stream);
        stream->writeUInt(dStart);
        stream->writeUInt(nSamples);
        bound.serialize(stream);
    }

    std::string toString() const
    {
        std::ostringstream oss;
        oss.setf(std::ios_base::boolalpha);

        oss << "GuidingRegion[" << endl
            << "  distribution = " << distribution.toString() << ",\n"
            << "  statistics = " << statistics.toString() << ",\n"
            << "  valid = " << valid << ",\n"
            << "  spatialSplitFlag = " << spatialSplitFlag << ",\n"
            << "  dStart = " << dStart << ",\n"
            << "  nSamples = " << nSamples << ",\n"
            << "  statsSinceLastSplit = " << statsSinceLastSplit.toString() << ",\n"
            << "  statsUntilLastSplit = " << statsUntilLastSplit.toString() << ",\n"
            << "  bound = " << bound.toString() << '\n'
            << "]";
        return oss.str();
    }

    bool valid {false};
    /// set on spatial split, as an indicator for the fitting procedure
    bool spatialSplitFlag {false};

    TDistribution distribution;
    TStatistics statistics;

    uint32_t dStart {0};
    uint32_t nSamples {0};

    //updated as samples are added, reset on split
    SampleStats statsSinceLastSplit;
    //updated on splits
    SampleStats statsUntilLastSplit;

    //NOTE: not really needed, but nice to have and not that large in relation to the whole region
    AABB bound;
};

GUIDING_NAMESPACE_END
