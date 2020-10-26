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

#ifndef PATHGUIDINGINTERSECTIONDATA_H
#define PATHGUIDINGINTERSECTIONDATA_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/shape.h>

GUIDING_NAMESPACE_BEGIN

struct ContributionAndThroughput
{
    Spectrum contribution {0.0f};
    Spectrum throughput   {1.0f};
    ContributionAndThroughput() {}
    explicit ContributionAndThroughput(const Spectrum& contribution, const Spectrum& throughput)
        : contribution{contribution}, throughput{throughput} {}

    Spectrum getClamped(const float maxThroughput)
    {
        if (throughput.max() > maxThroughput)
            return contribution*throughput*(maxThroughput/throughput.max());
        return contribution*throughput;
    }
};

struct PathGuidingIntersectionData
{
    Point pos;
    float distance;
    float cosThetaO;

    Vector wiWorld {0.0f};
    float pdfWiWorld;
    float cosThetaI;
    unsigned int sampledType {0};

    //direct light from next its, bsdf*miWeight
    ContributionAndThroughput bsdfDirectLight;
    float bsdfMiWeight     {1.0f};
    //direct light from next event estimation, bsdf*miWeight
    ContributionAndThroughput neeDirectLight;
    //self emission
    ContributionAndThroughput emission;

    //bsdfWeight/rrSurvivalProb, applies to following bounces only
    Spectrum throughputFactors {1.0f};

    float eta {1.0f};
    float roughness {1.0f};

    AABB regionBounds;
    AABB splattingVolume;

    PathGuidingIntersectionData(const Intersection& its)
        : pos{its.p}, distance{its.t}, cosThetaO{its.wi.z} {}

    PathGuidingIntersectionData& operator*=(Spectrum factor)
    {
        bsdfDirectLight.throughput *= factor;
        neeDirectLight.throughput  *= factor;
        emission.throughput        *= factor;

        return *this;
    }
};

GUIDING_NAMESPACE_END

#endif
