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
#include "GuidingTree.h"
#include "GuidingRegion.h"
#include "Range.h"

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/thread.h>
#include <mitsuba/core/util.h>

#include <atomic>

#include "GuidedBSDF.h"

// factories for the default types.
// make sure to include your own before instantiating the GuidingField
// if you use different distributions or statistics.
#include "GuidingFieldFactory.h"
#include "GuidingDistributionFactoryVMM.h"

GUIDING_NAMESPACE_BEGIN

template<typename TDistribution, typename TStatistics>
class GuidingField {

public:

    typedef TDistribution DistributionType;
    typedef TStatistics   StatisticsType;

    typedef GuidingRegion<TDistribution, TStatistics> GuidingRegionType;

    typedef GuidingTree<GuidingRegionType> GuidingTreeType;

    typedef GuidingFieldFactory<TDistribution, TStatistics>     GuidingFieldFactoryType;
    typedef GuidedBSDF<TDistribution, TStatistics>              GuidedBSDFType;

    GuidingField() = default;

    explicit GuidingField(const Properties& props)
        : m_guidingTree{typename GuidingTreeType::Settings{props}}, m_factory{props} {}

    explicit GuidingField(Stream *stream)
        : m_guidingTree{stream}, m_factory{stream},
          m_iteration{stream->readUInt()}, m_totalSPP{stream->readUInt()} {}

    void serialize(Stream* stream) const{
        m_guidingTree.serialize(stream);
        m_factory.serialize(stream);

        stream->writeUInt(m_iteration);
        stream->writeUInt(m_totalSPP);
    }

    void configureGuidedBSDF(GuidedBSDFType& gBSDF, const Intersection& its, const Vector wiWorld, const BSDF* surfaceBSDF) const {
        Guiding_Assert(its.isValid());
        Guiding_Assert(isValid(wiWorld));
        gBSDF.m_bsdf = surfaceBSDF;

        if (gBSDF.bsdfTypeAllowsGuiding() && m_guidingTree.getBounds().contains(its.p))
        {
            gBSDF.m_region = &m_guidingTree.getRegion(its.p);
            if (EXPECT_TAKEN(gBSDF.m_region->valid))
                m_factory.prepareGuidedBSDFForSampling(gBSDF, its, wiWorld);
            else
                gBSDF.m_guidingData.clear();
        }
        else
        {
            gBSDF.m_region = nullptr;
            gBSDF.m_guidingData.clear();
        }

        gBSDF.configure();
    }

    template<typename TSampleContainer>
    void buildField(const AABB &bounds, TSampleContainer& samples){
        std::vector<Point> zeroValuedSamples;
        initField(bounds);
        updateField(samples, zeroValuedSamples);
    }

    void initField(const AABB& bounds){
        m_guidingTree.initTree(bounds);
        m_iteration = 0;
        m_totalSPP  = 0;
    }

    template<typename TSampleContainer, typename TPointContainer>
    void updateField(TSampleContainer& samples, TPointContainer& zeroValuedSamples, const RenderJob* jobForProgress = 0){
        const bool updateFit = m_iteration;

#if defined(MTS_OPENMP)
        ref<Scheduler> scheduler = Scheduler::getInstance();
        size_t nCores = scheduler->getCoreCount();

        Thread::initializeOpenMP(nCores);
#endif
        Range<TSampleContainer> sampleRange{samples};
        Range<TPointContainer> zeroValuedSampleRange{zeroValuedSamples};

        std::atomic<uint32_t> updatedRegions{0};
        ProgressReporter treeUpdateProgress("Updating Guiding Field", m_guidingTree.getNumRegions()+samples.size()/m_guidingTree.m_settings.maxSamples, jobForProgress);
        auto regionUpdate = [this, &samples, updateFit, &treeUpdateProgress, &updatedRegions](GuidingRegionType& region) -> void
        {
            Range<TSampleContainer> sampleRange{samples.begin(), region.dStart, region.nSamples};
            if (updateFit)
                m_factory.updateFit(region, sampleRange);
            else
                m_factory.fit(region, sampleRange);

            treeUpdateProgress.update(++updatedRegions);
        };

        ref<Timer> treeTimer = new Timer();
        SLog(EInfo, "Begin Tree update");
        m_guidingTree.updateTree(sampleRange, zeroValuedSampleRange, nCores, regionUpdate);
        treeUpdateProgress.finish();
        SLog(EInfo, "Tree update time: %s", timeString(treeTimer->getSeconds(), true).c_str());

        SLog(EInfo, "Tree: %s\n", m_guidingTree.toString().c_str());
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "GuidingField [" << endl
            << "  BSPTree = " << m_guidingTree.toString() << ",\n"
            << "  factory = " << m_factory.toString() << ",\n"
            << "  iteration = " << m_iteration << ",\n"
            << "  totalSPP = " << m_totalSPP << ",\n"
            << "]";
        return oss.str();
    }

    void addTrainingIteration(uint32_t spp) {
        m_totalSPP += spp;
        ++m_iteration;
    }

public:
    GuidingTreeType m_guidingTree;
    GuidingFieldFactoryType m_factory;

    uint32_t m_iteration {0};
    uint32_t m_totalSPP  {0};
};

GUIDING_NAMESPACE_END
