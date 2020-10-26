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

#include <mitsuba/core/aabb.h>
#include <mitsuba/core/stream.h>

GUIDING_NAMESPACE_BEGIN

struct SampleStats
{
    Point mean {Point(0.0f)};
    /// variance times (numSamples-1)
    Vector unnormalizedVariance {Vector(0.0f)};
    uint32_t numSamples {0};
    // for flux estimation when using the forward training
    uint32_t numZeroValuedSamples {0};
    AABB sampleBound;

    SampleStats() = default;

    explicit SampleStats(Stream* stream)
        : mean(stream), unnormalizedVariance(stream), numSamples(stream->readUInt()), numZeroValuedSamples(stream->readUInt()), sampleBound(stream)
    {

    }

    void serialize(Stream* stream) const {
        mean.serialize(stream);
        unnormalizedVariance.serialize(stream);
        stream->writeUInt(numSamples);
        stream->writeUInt(numZeroValuedSamples);
        sampleBound.serialize(stream);
    }

    void clear() {
        *this = SampleStats{};
    }

    Vector computeVariance() const {
        if (numSamples > 1)
            return unnormalizedVariance/static_cast<Float>(numSamples);
        else
            return Vector{0.0f};
    }

    SampleStats operator+(const Point& sample) const {
        SampleStats resultStats{*this};
        resultStats += sample;

        return resultStats;
    }

    void operator+=(const Point& sample) {
        const Vector diffVector = sample-mean;
        mean += diffVector*(1.0f/static_cast<Float>(++numSamples));
        unnormalizedVariance.x += diffVector.x*(sample.x-mean.x);
        unnormalizedVariance.y += diffVector.y*(sample.y-mean.y);
        unnormalizedVariance.z += diffVector.z*(sample.z-mean.z);
        sampleBound.expandBy(sample);
    }

    SampleStats operator+(const SampleStats& other) const {
        SampleStats resultStats{*this};
        resultStats += other;

        return resultStats;
    }

    void operator+=(const SampleStats& other) {
        if (other.numSamples == 0)
            return;
        else if (numSamples == 0)
        {
            *this = other;
            return;
        }

        const Float weightSelf = numSamples;
        const Float weightOther = other.numSamples;
        const Float sumWeight = numSamples+other.numSamples;

        const Point meanSelfSqr {mean.x*mean.x, mean.y*mean.y, mean.z*mean.z};
        const Point meanOtherSqr {other.mean.x*other.mean.x, other.mean.y*other.mean.y, other.mean.z*other.mean.z};

        mean = (mean*weightSelf+other.mean*weightOther)/sumWeight;

        const Point combinedMeanSqr {mean.x*mean.x, mean.y*mean.y, mean.z*mean.z};

        unnormalizedVariance = weightSelf*meanSelfSqr+unnormalizedVariance+weightOther*meanOtherSqr+other.unnormalizedVariance-sumWeight*combinedMeanSqr;
        numSamples += other.numSamples;
        numZeroValuedSamples += other.numZeroValuedSamples;
        sampleBound.expandBy(other.sampleBound);
    }

    SampleStats operator*(const float factor) const {
        SampleStats resultStats{*this};
        resultStats *= factor;

        return resultStats;
    }

    void operator*=(const float factor) {
        if (EXPECT_NOT_TAKEN(!(factor >= 0.0f)))
        {
            SLog(EError, "invalid scaling factor %f.", factor);
            return;
        }

        numSamples *= factor;
        numZeroValuedSamples *= factor;
        unnormalizedVariance *= factor;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SampleStats[" << endl
            << "  mean = " << mean.toString() << "," << endl
            << "  variance = " << computeVariance().toString() << "," << endl
            << "  numSamples = " << numSamples << "," << endl
            << "  numZeroValuedSamples = " << numZeroValuedSamples << "," << endl
            << "  sampleBound = " << sampleBound.toString() << endl
            << "]";
        return oss.str();
    }
};

GUIDING_NAMESPACE_END
