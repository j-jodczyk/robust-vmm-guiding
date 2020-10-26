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

#ifndef GUIDINGFIELDFACTORY_H
#define GUIDINGFIELDFACTORY_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/statistics.h>

#include <pmm/pmm.h>

#include <mitsuba/guiding/guiding.h>
#include <mitsuba/guiding/GuidingRegion.h>
#include <mitsuba/guiding/GuidingFieldFactoryFwd.h>
#include <mitsuba/guiding/GuidedBSDF.h>

#include <string>
#include <sstream>
#ifdef GUIDING_DETAILED_STATISTICS
#include <chrono>
#endif

GUIDING_NAMESPACE_BEGIN

#ifdef GUIDING_DETAILED_STATISTICS
//NOTE: too costly to keep enabled by default
//static StatsCounter unassignedSamples("GuidingFieldFactory", "Unassigned samples", EAverage);
static StatsCounter avgFitDuration("GuidingFieldFactory", "Average fitting duration [us]", EAverage);
#endif

template<typename TDistribution, typename TStatistics>
class GuidingFieldFactory
{
private:
    using TDistributionFactory = GuidingDistributionFactory<TDistribution>;
    using TStatisticsFactory   = GuidingStatsFactory<TDistribution, TStatistics>;

public:
    using TGuidingRegion = GuidingRegion<TDistribution, TStatistics>;
    typedef TDistributionFactory DistributionFactoryType;
    typedef TStatisticsFactory StatisticsFactoryType;

    GuidingFieldFactory(const Properties& props=Properties{})
        : m_distributionFactory{TDistributionFactory(props)}, m_statisticsFactory{props},
          m_decayOnSpatialSplit{props.getFloat("decayOnSpatialSplit", 0.25f)} {}

    explicit GuidingFieldFactory(Stream* stream)
        : m_distributionFactory{stream}, m_statisticsFactory{stream},
          m_decayOnSpatialSplit{stream->readSingle()} {}

    void serialize(Stream* stream) const
    {
        m_distributionFactory.serialize(stream);
        m_statisticsFactory.serialize(stream);
        stream->writeSingle(m_decayOnSpatialSplit);
    }

    template<typename TSampleRange>
    void fit(TGuidingRegion& region, TSampleRange& samples) const
    {
        region.spatialSplitFlag = false;

        region.statistics.clear();
        m_distributionFactory.initialize(region.distribution);

        if (samples.size() < region.distribution.getK()*2)
        {
            region.valid = false;
            return;
        }

        SampleStats combinedStats = region.statsUntilLastSplit+region.statsSinceLastSplit;
        combinedStats.sampleBound = region.statsSinceLastSplit.sampleBound;

        m_statisticsFactory.preFit(region.statistics, region.distribution, samples, combinedStats);
#ifdef GUIDING_DETAILED_STATISTICS
            const std::chrono::high_resolution_clock::time_point fitStart = std::chrono::high_resolution_clock::now();
            avgFitDuration.incrementBase();
#endif
        m_distributionFactory.fit(region.distribution, samples);
#ifdef GUIDING_DETAILED_STATISTICS
        const std::chrono::high_resolution_clock::time_point fitEnd = std::chrono::high_resolution_clock::now();
        avgFitDuration += static_cast<size_t>(std::chrono::duration_cast<std::chrono::microseconds>(fitEnd-fitStart).count());
#endif
        region.valid = m_distributionFactory.checkValidity(region.distribution);

        if (EXPECT_TAKEN(region.valid))
        {
            m_statisticsFactory.postFit(region.statistics, region.distribution, samples, combinedStats, m_distributionFactory);
            region.valid = m_distributionFactory.checkValidity(region.distribution);
        }
    }

    template<typename TSampleRange>
    void updateFit(TGuidingRegion& region, TSampleRange& samples) const
    {
        region.distribution.m_numEMIterations = 0;

        if (samples.size() < region.distribution.getK()*2)
            return;

        if (EXPECT_NOT_TAKEN(!region.valid || !region.distribution.m_numSamples))
        {
            SLog(EWarn, "tried to update an invalid mixture, falling back to regular fit.");

            fit(region, samples);
            return;
        }

        if (region.spatialSplitFlag)
        {
            m_distributionFactory.onSpatialSplit(region.distribution);
            m_statisticsFactory.onSpatialSplit(region.statistics);

            region.spatialSplitFlag = false;
        }

#ifdef GUIDING_DETAILED_STATISTICS
        /*
        unassignedSamples.incrementBase(samples.size());
        unassignedSamples += std::count_if(samples.begin(), samples.end(), [&region](const typename TSampleRange::value_type& sample) -> bool {return !(region.distribution.pdf(sample.direction) > PMM_EPSILON);});
        */
#endif

        SampleStats combinedStats = region.statsUntilLastSplit+region.statsSinceLastSplit;
        combinedStats.sampleBound = region.statsSinceLastSplit.sampleBound;

        m_statisticsFactory.preFit(region.statistics, region.distribution, samples, combinedStats);

#ifdef GUIDING_DETAILED_STATISTICS
        const std::chrono::high_resolution_clock::time_point fitStart = std::chrono::high_resolution_clock::now();
        avgFitDuration.incrementBase();
#endif
        m_distributionFactory.updateFit(region.distribution, samples);
#ifdef GUIDING_DETAILED_STATISTICS
        const std::chrono::high_resolution_clock::time_point fitEnd = std::chrono::high_resolution_clock::now();
        avgFitDuration += static_cast<size_t>(std::chrono::duration_cast<std::chrono::microseconds>(fitEnd-fitStart).count());
#endif
        region.valid = m_distributionFactory.checkValidity(region.distribution);

        if (EXPECT_TAKEN(region.valid))
        {
            m_statisticsFactory.postFit(region.statistics, region.distribution, samples, combinedStats, m_distributionFactory);
            region.valid = m_distributionFactory.checkValidity(region.distribution);
        }
    }

    void prepareGuidedBSDFForSampling(GuidedBSDF<TDistribution, TStatistics>& gBSDF, const Intersection& its, const Vector wiWorld) const
    {
        Guiding_Assert(gBSDF.m_region->valid);

        gBSDF.m_guidingData.m_liDistribution = gBSDF.m_region->distribution;

        m_statisticsFactory.preSample(gBSDF.m_region->statistics, gBSDF.m_guidingData.m_liDistribution, its.p);
        m_distributionFactory.prepareGuidedBSDFForSampling(gBSDF, its, wiWorld);
    }

    std::string toString() const
    {
        std::ostringstream oss;

        oss.setf(std::ios_base::boolalpha);

        oss << "GuidingFieldFactory[\n"
            << "  distributionFactory: " << m_distributionFactory.toString() << '\n'
            << "  statisticsFactory: " << m_statisticsFactory.toString() << '\n'
            << "  decayOnSpatialSplit: " << m_decayOnSpatialSplit << '\n'
            << ']';

        return oss.str();
    }

public:
    TDistributionFactory m_distributionFactory;
    TStatisticsFactory m_statisticsFactory;

    float m_decayOnSpatialSplit;
};

GUIDING_NAMESPACE_END

#endif // GUIDINGFIELDFACTORY_H
