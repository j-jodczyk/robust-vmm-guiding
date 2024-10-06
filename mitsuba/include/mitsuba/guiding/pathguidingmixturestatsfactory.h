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

#ifndef PATHGUIDINGMIXTURESTATSFACTORY_H
#define PATHGUIDINGMIXTURESTATSFACTORY_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/matrix.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/sampler.h>

#include <mitsuba/guiding/incrementalcovariance2d.h>
#include <mitsuba/guiding/incrementaldistance.h>
#include <mitsuba/guiding/incrementalpearsonchisquared.h>
#include <mitsuba/guiding/pathguidingmixturestats.h>

#include <pmm/pmm.h>
#include <pmm/ParametricMixtureModel.h>
#include <pmm/VMFKernel.h>

#include <mitsuba/guiding/guiding.h>
#include <mitsuba/guiding/SampleStats.h>
#include <mitsuba/guiding/GuidingRegion.h>
#include <mitsuba/guiding/GuidingFieldFactoryFwd.h>

#ifdef GUIDING_DETAILED_STATISTICS
#include <chrono>
#endif

GUIDING_NAMESPACE_BEGIN

#ifdef GUIDING_DETAILED_STATISTICS
static StatsCounter preventedMerges("GuidingStatsFactory<PathGuidingMixtureStats>", "Merges prevented by 3-point safety", EPercentage);
static StatsCounter avgNumMerges("GuidingStatsFactory<PathGuidingMixtureStats>", "Average num merges", EAverage);
static StatsCounter avgNumSplits("GuidingStatsFactory<PathGuidingMixtureStats>", "Average num splits", EAverage);
static StatsCounter avgSplitDuration("GuidingStatsFactory<PathGuidingMixtureStats>", "Average split duration (combined) [us]", EAverage);
static StatsCounter avgMergeDuration("GuidingStatsFactory<PathGuidingMixtureStats>", "Average merge duration (combined) [ns]", EAverage);
#endif

template<typename TScalar, uint32_t NKernels>
class GuidingStatsFactory<lightpmm::ParametricMixtureModel<lightpmm::VMFKernel<TScalar>, NKernels>, PathGuidingMixtureStats<lightpmm::ParametricMixtureModel<lightpmm::VMFKernel<TScalar>, NKernels>>>
{
private:
    using TKernel = lightpmm::VMFKernel<TScalar>;
    using TDistribution = lightpmm::ParametricMixtureModel<TKernel, NKernels>;
    using TDistributionFactory = GuidingDistributionFactory<TDistribution>;
    using TStatistics = PathGuidingMixtureStats<TDistribution>;
public:
    GuidingStatsFactory(const Properties& props=Properties{})
    {
        m_splitMinDivergence            = props.getFloat(  "splitMinDivergence",            m_splitMinDivergence);
        m_mergeMaxDivergence            = props.getFloat(  "mergeMaxDivergence",            m_mergeMaxDivergence);
        m_splitAndMerge                 = props.getBoolean("splitAndMerge",                 m_splitAndMerge);
        m_minSamplesForMerging          = props.getSize(   "minSamplesForMerging",          m_minSamplesForMerging);
        m_minSamplesForSplitting        = props.getSize(   "minSamplesForSplitting",        m_minSamplesForSplitting);
        m_minSamplesForPostSplitFitting = props.getSize(   "minSamplesForPostSplitFitting", m_minSamplesForPostSplitFitting);
        m_decayOnSpatialSplit           = props.getFloat(  "decayOnSpatialSplit",           m_decayOnSpatialSplit);

        m_parallaxCompensation          = props.getBoolean("parallaxCompensation",          m_parallaxCompensation);
        m_safetyMerge                   = props.getBoolean("safetyMerge",                   m_safetyMerge);
    }

    explicit GuidingStatsFactory(Stream* stream)
          : m_splitMinDivergence{stream->readSingle()}, m_mergeMaxDivergence{stream->readSingle()}, m_splitAndMerge{stream->readBool()},
            m_minSamplesForMerging{stream->readSize()}, m_minSamplesForSplitting{stream->readSize()}, m_minSamplesForPostSplitFitting{stream->readSize()},
            m_decayOnSpatialSplit{stream->readSingle()}, m_parallaxCompensation{stream->readBool()}, m_safetyMerge{stream->readBool()}
    {

    }
    void serialize(Stream* stream) const
    {
        stream->writeSingle(m_splitMinDivergence);
        stream->writeSingle(m_mergeMaxDivergence);
        stream->writeBool(m_splitAndMerge);
        stream->writeSize(m_minSamplesForMerging);
        stream->writeSize(m_minSamplesForSplitting);
        stream->writeSize(m_minSamplesForPostSplitFitting);
        stream->writeSingle(m_decayOnSpatialSplit);
        stream->writeBool(m_parallaxCompensation);
        stream->writeBool(m_safetyMerge);
    }

    ///to be called on spatial splits (only once per tree update, even for repeated splits)
    void onSpatialSplit(TStatistics& stats) const
    {
        stats *= m_decayOnSpatialSplit;
        stats.samplesSinceLastMerge = 0;
    }

    ///to be called before the PMM is fit
    template<typename TSampleRange>
    void preFit(TStatistics& stats, TDistribution& pmm, TSampleRange& samples, const Guiding::SampleStats statistics) const
    {
        if (m_parallaxCompensation)
        {
            //parallax shift samples to new mean
            {
                const Point parallaxMean = statistics.mean;

                for (auto& sample : samples)
                    repositionSample(sample, parallaxMean);
            }
            //parallax shift mixture to new mean
            if (stats.lastFitStats.numSamples > 0)
            {
                const Point oldParallaxMean = stats.lastFitStats.mean;
                const Point newParallaxMean = statistics.mean;
                stats.incrementalDistance.reposition(pmm, oldParallaxMean-newParallaxMean);
            }
        }

        stats.lastFitStats = statistics;
    }
    ///to be called after the PMM has been fit
    template<typename TSampleRange>
    void postFit(TStatistics& stats, TDistribution& pmm, TSampleRange& samples, const Guiding::SampleStats statistics, const TDistributionFactory& pmmFactory) const
    {
        if (samples.size() > 0){
            if (statistics.numSamples+statistics.numZeroValuedSamples > 0)
                stats.flux = Spectrum{pmm.m_sampleWeight*static_cast<float>(statistics.numSamples)/(pmm.m_numSamples*static_cast<float>(statistics.numSamples+statistics.numZeroValuedSamples))};
            else
                stats.flux = Spectrum{0.0f};
        }

        if (m_splitAndMerge)
        {
#ifdef GUIDING_DETAILED_STATISTICS
            avgNumMerges.incrementBase();
            avgNumSplits.incrementBase();
            avgSplitDuration.incrementBase();
#endif
            stats.samplesSinceLastMerge += samples.size();

            if (stats.samplesSinceLastMerge >= m_minSamplesForMerging)
            {
#ifdef GUIDING_DETAILED_STATISTICS
                avgMergeDuration.incrementBase();
                const std::chrono::high_resolution_clock::time_point mergeStart = std::chrono::high_resolution_clock::now();
#endif

                pmm.removeWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);
                //while (merge(stats, pmm));
                mergeAll(stats, pmm);
                pmm.applyWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);

                stats.samplesSinceLastMerge = 0;

#ifdef GUIDING_DETAILED_STATISTICS
                const std::chrono::high_resolution_clock::time_point mergeEnd = std::chrono::high_resolution_clock::now();
                avgMergeDuration += static_cast<size_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(mergeEnd-mergeStart).count());
#endif
            }

#ifdef GUIDING_DETAILED_STATISTICS
            const std::chrono::high_resolution_clock::time_point splitStart = std::chrono::high_resolution_clock::now();
#endif

            //update incremental statistics after merging
            stats.incrementalPearsonChiSquared.updateDivergence(pmm, samples);
            stats.incrementalCovariance2D.updateStatistics(pmm, samples);

            const bool firstFit = samples.size() == pmm.m_totalNumSamples;
            const bool firstFitOrEnoughSamplesForSplitting = firstFit || samples.size() > m_minSamplesForPostSplitFitting;

            pmm.removeWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);
            //while (split(stats, pmm, samples, pmmFactory, firstFitOrEnoughSamplesForSplitting));
            splitAll(stats, pmm, samples, pmmFactory, true, firstFitOrEnoughSamplesForSplitting);
            pmm.applyWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);

#ifdef GUIDING_DETAILED_STATISTICS
            const std::chrono::high_resolution_clock::time_point splitEnd = std::chrono::high_resolution_clock::now();
            avgSplitDuration += static_cast<size_t>(std::chrono::duration_cast<std::chrono::microseconds>(splitEnd-splitStart).count());
#endif
        }

        if (m_parallaxCompensation)
            stats.incrementalDistance.updateDistances(pmm, samples);
    }

    ///to be called before the PMM is being sampled
    void preSample(const TStatistics& stats, TDistribution& pmm, const Point samplePos) const
    {
        //parallax shift mixture to sample position
        if (m_parallaxCompensation)
            stats.incrementalDistance.repositionDistribution(pmm, stats.lastFitStats.mean-samplePos);
    }

    std::string toString() const
    {
        std::ostringstream oss;

        oss << "GuidingStatsFactory<PathGuidingMixtureStats> [\n"
            << "  splitMinDivergence = " << m_splitMinDivergence << '\n'
            << "  mergeMaxDivergence = " << m_mergeMaxDivergence << '\n'
            << "  splitAndMerge = " << m_splitAndMerge << '\n'
            << "  minSamplesForMerging = " << m_minSamplesForMerging << '\n'
            << "  minSamplesForSplitting = " << m_minSamplesForSplitting << '\n'
            << "  minSamplesForPostSplitFitting = " << m_minSamplesForPostSplitFitting << '\n'
            << "  decayOnSpatialSplit = " << m_decayOnSpatialSplit << '\n'
            << "  parallaxCompensation = " << m_parallaxCompensation << '\n'
            << "  safetyMerge = " << m_safetyMerge << '\n'
            << ']';

        return oss.str();
    }

public:
    float m_splitMinDivergence {0.5f};
    float m_mergeMaxDivergence {0.025f};
    bool m_splitAndMerge  {true};
    size_t m_minSamplesForMerging {8192};
    size_t m_minSamplesForSplitting {4096};
    size_t m_minSamplesForPostSplitFitting {4096};
    float m_decayOnSpatialSplit {0.25f};

    bool m_parallaxCompensation {true};
    bool m_safetyMerge {true};

    FINLINE static void mergeComponents(TStatistics& stats, TDistribution &pmm, const uint32_t componentA, const uint32_t componentB)
    {
#ifdef GUIDING_VALIDATE_INTERMEDIATE_MIXTURE_STATES
        const uint32_t numComponents = pmm.getK();

        if (!pmm.valid(pmmFactory.m_vmmFactory.computeMinComponentWeight(numComponents)))
            SLog(EWarn, "invalid mixture state before merging components %u and %u", componentA, componentB);
#endif

        stats.incrementalDistance.merge(pmm, componentA, componentB);
        stats.incrementalPearsonChiSquared.merge(pmm, componentA, componentB);
        stats.incrementalCovariance2D.merge(pmm, componentA, componentB);

        pmm.mergeComponents(componentA, componentB);

#ifdef GUIDING_VALIDATE_INTERMEDIATE_MIXTURE_STATES
        if (!pmm.valid(pmmFactory.m_vmmFactory.computeMinComponentWeight(numComponents-1)))
            SLog(EWarn, "error occured after merging components %u and %u", componentA, componentB);
#endif
    }

    bool mergeOne(TStatistics& stats, TDistribution& pmm) const
    {
        const uint32_t numComponents = pmm.getK();

        if (EXPECT_NOT_TAKEN(numComponents <= 1))
            return false;

        const uint32_t numSimilarityValues = numComponents*(numComponents-1)/2;

        const std::array<float, TDistribution::MaxK::value*(TDistribution::MaxK::value-1)/2> mergeMetric = computePearsonChiSquaredMergeMetric(pmm);
        const auto minElementIterator = std::min_element(mergeMetric.begin(), mergeMetric.begin()+numSimilarityValues);
        if (*minElementIterator > m_mergeMaxDivergence)
            return false;

        const uint32_t maxElementPos = static_cast<uint32_t>(std::distance(mergeMetric.begin(), minElementIterator));

        uint32_t componentA = 0;
        uint32_t componentBInPairsPerComponentA = maxElementPos;

        for (uint32_t stride=numComponents-1; componentBInPairsPerComponentA >= stride; componentBInPairsPerComponentA-=stride--)
            componentA++;

        const uint32_t componentB = componentA+componentBInPairsPerComponentA+1;

        mergeComponents(stats, pmm, componentA, componentB);

#ifdef GUIDING_DETAILED_STATISTICS
        ++avgNumMerges;
#endif

        return true;
    }

    uint32_t mergeAll(TStatistics& stats, TDistribution& pmm, const bool iterative=true) const
    {
        uint32_t totalNumMerges = 0;

        do
        {
            const uint32_t numComponents = pmm.getK();
            uint32_t numMerges = 0;

            if (EXPECT_NOT_TAKEN(numComponents <= 1))
                return false;

            const uint32_t numSimilarityValues = numComponents*(numComponents-1)/2;

            const std::array<float, TDistribution::MaxK::value*(TDistribution::MaxK::value-1)/2> mergeMetric = computePearsonChiSquaredMergeMetric(pmm);

            std::array<std::pair<float, std::pair<uint32_t, uint32_t>>, TDistribution::MaxK::value*(TDistribution::MaxK::value-1)/2> mergeCandidates;
            for (uint32_t componentA=0, offsetA=0; componentA < numComponents; ++componentA, offsetA+=numComponents-componentA)
            {
                for (uint32_t componentB=componentA+1, offsetB=0; componentB < numComponents; ++componentB, ++offsetB)
                {
                    mergeCandidates[offsetA+offsetB] = std::make_pair(mergeMetric[offsetA+offsetB], std::make_pair(componentA, componentB));
                }
            }

            if (m_safetyMerge)
            {
#ifdef GUIDING_DETAILED_STATISTICS
                preventedMerges.incrementBase(std::count_if(mergeCandidates.begin(), mergeCandidates.begin()+numSimilarityValues, [this](std::pair<float, std::pair<uint32_t, uint32_t>> candidate) -> bool {
                    return candidate.first <= m_mergeMaxDivergence;
                }));
#endif
                const Vector extents = stats.lastFitStats.sampleBound.getExtents();

                // for performance reasons 3 pre-sampled random vectors are are stored in a const array.
                const Vector randomOffsets[] = {Vector( 0.487672f,  0.394023f,  0.437635f),
                                                Vector(-0.262328f, -0.355977f, -0.312365f),
                                                Vector( 0.237672f,  0.144023f,  0.187635f)};

                /*
                // correct code for initializing and using a random sampler
                Properties samplerProps = Properties("ldsampler");
                samplerProps.setSize("sampleCount", 3);
                ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(samplerProps));
                sampler->configure();
                sampler->generate(Point2i(0,0));

                Vector randomOffsets[3];
                randomOffsets[0] =  Vector(0.5f-sampler->next1D(),0.5f-sampler->next1D(),0.5f-sampler->next1D());
                sampler->advance();
                randomOffsets[1] =  Vector(0.5f-sampler->next1D(),0.5f-sampler->next1D(),0.5f-sampler->next1D());
                sampler->advance();
                randomOffsets[2] =  Vector(0.5f-sampler->next1D(),0.5f-sampler->next1D(),0.5f-sampler->next1D());
                sampler->advance();
                */
                for (int i=0; i<3; ++i)
                {
                    TDistribution pmmCopy{pmm};
                    stats.incrementalDistance.repositionDistribution(pmmCopy, Vector(extents.x*randomOffsets[i].x,extents.y*randomOffsets[i].y,extents.z*randomOffsets[i].z));
                    const std::array<float, TDistribution::MaxK::value*(TDistribution::MaxK::value-1)/2> mergeMetric = computePearsonChiSquaredMergeMetric(pmmCopy);

                    for (uint32_t k=0; k < numSimilarityValues; ++k)
                    {
                        if (mergeMetric[k] > mergeCandidates[k].first)
                        {
#ifdef GUIDING_DETAILED_STATISTICS
                            if (mergeMetric[k] > m_mergeMaxDivergence && mergeCandidates[k].first <= m_mergeMaxDivergence)
                                ++preventedMerges;
#endif
                            mergeCandidates[k].first = mergeMetric[k];
                        }
                    }
                }
            }

            const auto validCandidatesEnd = std::partition(mergeCandidates.begin(), mergeCandidates.begin()+numSimilarityValues, [this](std::pair<float, std::pair<uint32_t, uint32_t>> candidate) -> bool {
                return candidate.first <= m_mergeMaxDivergence;
            });

            std::sort(mergeCandidates.begin(), validCandidatesEnd, [](std::pair<float, std::pair<uint32_t, uint32_t>> a, std::pair<float, std::pair<uint32_t, uint32_t>> b) -> bool {
                return a.first < b.first;
            });

            std::array<uint32_t, (TDistribution::MaxK::value+31)/32> bitmask;
            std::fill(bitmask.begin(), bitmask.end(), 0U);

            std::array<std::pair<uint32_t, uint32_t>, TDistribution::MaxK::value> merges;

            for (auto it=mergeCandidates.begin(); it != validCandidatesEnd; ++it)
            {
                const uint32_t componentA = std::min(it->second.first, it->second.second);
                const uint32_t componentB = std::max(it->second.first, it->second.second);

                if ((bitmask[componentA/32]&(1<<(componentA%32))) != 0 || (bitmask[componentB/32]&(1<<(componentB%32))) != 0)
                    continue;

                merges[numMerges++] = std::make_pair(componentA, componentB);
                bitmask[componentA/32] |= (1<<(componentA%32));
                bitmask[componentB/32] |= (1<<(componentB%32));
            }

            if (numMerges == 0)
                break;

            //this is needed to avoid updating indices after individual merges
            std::sort(merges.begin(), merges.begin()+numMerges, [](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t> b) -> bool {
                return a.second > b.second;
            });

            for (uint32_t i=0; i<numMerges; ++i)
            {
                SAssert(merges[i].second > merges[i].first);
                mergeComponents(stats, pmm, merges[i].first, merges[i].second);
            }

            totalNumMerges += numMerges;
        }
        while (iterative);

#ifdef GUIDING_DETAILED_STATISTICS
        avgNumMerges += totalNumMerges;
#endif

        return totalNumMerges;
    }

    static void splitComponentUsingPCA(TStatistics& stats, TDistribution& pmm, const uint32_t component, const float maxKappa)
    {
        if (EXPECT_NOT_TAKEN(pmm.getK() == TDistribution::MaxK::value))
        {
            SLog(EWarn, "Skipping split of component, since the mixture has reached its maximum size.");
            return;
        }

#ifdef GUIDING_VALIDATE_INTERMEDIATE_MIXTURE_STATES
        if (!pmm.valid(pmmFactory.m_vmmFactory.computeMinComponentWeight(pmm.getK())))
            SLog(EWarn, "invalid mixture state before splitting component %u", component);
#endif

        const uint32_t sourceIndex = component;
        const uint32_t targetIndex = pmm.getK();

        TKernel& sourceComponent = pmm.getComponent(sourceIndex/TScalar::Width::value);
        TKernel& targetComponent = pmm.getComponent(targetIndex/TScalar::Width::value);

        //Log(EDebug, "splitting component %zu", index);

        const Frame frame{sourceComponent.getMu(component%TScalar::Width::value)};

        const Matrix2x2 covariance = stats.incrementalCovariance2D.computeCovarianceMatrix(component);

        Matrix2x2 eigVec;
        Float eigVal[2];
        covariance.symEig(eigVec, eigVal);

        const int maxEigValIndex = eigVal[1] > eigVal[0];
        const float oneHalfSigmaOffset = std::min(1.0f, 0.5f*sqrt(eigVal[maxEigValIndex]));

        Vector2 principalComponentDir = eigVec.col(maxEigValIndex);
        principalComponentDir *= oneHalfSigmaOffset;

        const float z = std::sqrt(1.0f-oneHalfSigmaOffset*oneHalfSigmaOffset);

        const Vector splitMuA {frame.toWorld(Vector{ principalComponentDir.x,  principalComponentDir.y, z})};
        const Vector splitMuB {frame.toWorld(Vector{-principalComponentDir.x, -principalComponentDir.y, z})};

        pmm.setK(targetIndex+1);

        const float sourceAvgCosine = sourceComponent.getR()[sourceIndex%TScalar::Width::value];
        const float sourceWeight = sourceComponent.getWeight(sourceIndex%TScalar::Width::value);

        const float splitWeight = sourceWeight*0.5f;
        //using this concentration, an immediate merge would recreate the original component
        const float maxAvgCosine = lightpmm::kappaToMeanCosine(maxKappa);
        const float splitAvgCosine = (z > 0.0f) ? std::min(maxAvgCosine, sourceAvgCosine/z) : maxAvgCosine;
        const float splitKappa = lightpmm::meanCosineToKappa(splitAvgCosine);

        sourceComponent.setKappaAndR(sourceIndex%TScalar::Width::value, splitKappa, splitAvgCosine);
        sourceComponent.setMu(sourceIndex%TScalar::Width::value,        splitMuA);
        sourceComponent.setWeight(sourceIndex%TScalar::Width::value,    splitWeight);

        targetComponent.setKappaAndR(targetIndex%TScalar::Width::value, splitKappa, splitAvgCosine);
        targetComponent.setMu(targetIndex%TScalar::Width::value,        splitMuB);
        targetComponent.setWeight(targetIndex%TScalar::Width::value,    splitWeight);

        stats.incrementalDistance.split(pmm, component);
        stats.incrementalPearsonChiSquared.split(pmm, component);
        stats.incrementalCovariance2D.split(pmm, component);

#ifdef GUIDING_VALIDATE_INTERMEDIATE_MIXTURE_STATES
        if (!pmm.valid(pmmFactory.m_vmmFactory.computeMinComponentWeight(pmm.getK())))
            SLog(EWarn, "error occured after splitting component %u", component);
#endif
    }

    template<typename TSampleRange>
    bool splitOne(TStatistics& stats, TDistribution& pmm, const TSampleRange& samples, const TDistributionFactory& pmmFactory, const bool fit=true) const
    {
        if (EXPECT_NOT_TAKEN(pmm.getK() >= TDistribution::MaxK::value))
            return false;

        const bool firstFit = samples.size() == pmm.m_totalNumSamples;
        const uint32_t numActiveKernels = (pmm.getK()+TScalar::Width::value-1)/TScalar::Width::value;

        const std::array<TScalar, TDistribution::NumKernels::value> divergence = stats.incrementalPearsonChiSquared.computeDivergence(pmm);
        std::array<float, TDistribution::MaxK::value> splitCriterion;
        for (uint32_t k=0; k<numActiveKernels; ++k)
        {
            const typename TScalar::BooleanType seenEnoughSamples = stats.incrementalPearsonChiSquared.numSamples[k] > m_minSamplesForSplitting;
            const TScalar weightedDivergence = lightpmm::ifthen(firstFit || seenEnoughSamples, divergence[k]*pmm.getComponent(k).m_weights, 0.0f);
            for (uint32_t i=0; i<TScalar::Width::value; ++i)
                splitCriterion[k*TScalar::Width::value+i] = fabs(weightedDivergence[i]);
        }

        const auto maxElementIterator = std::max_element(splitCriterion.begin(), splitCriterion.begin()+pmm.getK());

        const float highestSplitCriterion = *maxElementIterator;
        if (!(highestSplitCriterion >= m_splitMinDivergence))
            return false;

        const uint32_t indexOfHighestSplitCriterion = static_cast<uint32_t>(std::distance(splitCriterion.begin(), maxElementIterator));

        splitComponentUsingPCA(stats, pmm, indexOfHighestSplitCriterion, pmmFactory.m_vmmFactory.m_properties.maxKappa);

        //repeated splitting requires fitting after each split unless all splits are determined beforehand.
        //otherwise, the same component is split over and over again and fitting all the split components becomes increasingly difficult

        //using the masked fit to fit only the splitted components
        if (fit)
        {
            pmm.applyWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);
            pmmFactory.m_vmmFactory.maskedFit(samples.begin(), samples.end(), pmm, computeMaskFromComponentIndexVector(std::vector<uint32_t>{{indexOfHighestSplitCriterion, pmm.getK()-1}}));
            pmm.removeWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);
        }

#ifdef GUIDING_DETAILED_STATISTICS
        ++avgNumSplits;
#endif

        return true;
    }

    template<typename TSampleRange>
    uint32_t splitAll(TStatistics& stats, TDistribution& pmm, const TSampleRange& samples, const TDistributionFactory& pmmFactory, const bool iterative=true, const bool fit=true) const
    {
        if (EXPECT_NOT_TAKEN(pmm.getK() >= TDistribution::MaxK::value))
            return false;

        const bool firstFit = samples.size() == pmm.m_totalNumSamples;
        uint32_t totalNumSplits = 0;
        uint32_t numSplits = 0;

        do
        {
            const uint32_t numActiveKernels = (pmm.getK()+TScalar::Width::value-1)/TScalar::Width::value;

            const std::array<TScalar, TDistribution::NumKernels::value> divergence = stats.incrementalPearsonChiSquared.computeDivergence(pmm);
            std::array<std::pair<float, uint32_t>, TDistribution::MaxK::value> splitCandidates;
            for (uint32_t k=0; k<numActiveKernels; ++k)
            {
                const typename TScalar::BooleanType seenEnoughSamples = stats.incrementalPearsonChiSquared.numSamples[k] > m_minSamplesForSplitting;
                const TScalar weightedDivergence = lightpmm::ifthen(firstFit || seenEnoughSamples, divergence[k]*pmm.getComponent(k).m_weights, 0.0f);
                for (uint32_t i=0; i<TScalar::Width::value; ++i)
                    splitCandidates[k*TScalar::Width::value+i] = std::make_pair(fabs(weightedDivergence[i]), k*TScalar::Width::value+i);
            }

            const auto lastCandidateIterator = std::partition(splitCandidates.begin(), splitCandidates.begin()+pmm.getK(),
                                                              [this](const std::pair<float, uint32_t> divergenceAndIndex) -> bool { return divergenceAndIndex.first >= m_splitMinDivergence; });

            const uint32_t numSplitCandidates = std::distance(splitCandidates.begin(), lastCandidateIterator);
            const uint32_t maxNumSplits = TDistribution::MaxK::value-pmm.getK();
            numSplits = std::min(numSplitCandidates, maxNumSplits);

            if (numSplits == 0)
                break;

            if (numSplitCandidates > numSplits)
                std::partial_sort(splitCandidates.begin(), splitCandidates.begin()+numSplits, lastCandidateIterator,
                                  [](const std::pair<float, uint32_t> a, const std::pair<float, uint32_t> b) -> bool { return a.first > b.first; });

            std::array<typename TScalar::BooleanType, TDistribution::NumKernels::value> modifedComponentMask;
            std::fill(modifedComponentMask.begin(), modifedComponentMask.end(), typename TScalar::BooleanType{false});

            for (uint32_t i=0; i<numSplits; ++i)
            {
                modifedComponentMask[splitCandidates[i].second/TScalar::Width::value].insert(splitCandidates[i].second%TScalar::Width::value, true);
                modifedComponentMask[pmm.getK()/TScalar::Width::value].insert(pmm.getK()%TScalar::Width::value, true);

                splitComponentUsingPCA(stats, pmm, splitCandidates[i].second, pmmFactory.m_vmmFactory.m_properties.maxKappa);
            }

            //repeated splitting requires fitting after each split unless all splits are determined beforehand.
            //otherwise, the same component is split over and over again and fitting all the split components becomes increasingly difficult

            //using the masked fit to fit only the splitted components
            if (fit)
            {
                pmm.applyWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);
                pmmFactory.m_vmmFactory.maskedFit(samples.begin(), samples.end(), pmm, modifedComponentMask);
                stats.incrementalPearsonChiSquared.updateDivergenceMasked(pmm, samples, modifedComponentMask);
                stats.incrementalCovariance2D.updateStatisticsMasked(pmm, samples, modifedComponentMask);
                pmm.removeWeightPrior(pmmFactory.m_vmmFactory.m_properties.vPrior);
            }

            totalNumSplits += numSplits;
        }
        while (iterative);

#ifdef GUIDING_DETAILED_STATISTICS
            avgNumSplits += totalNumSplits;
#endif

        return totalNumSplits;
    }

    static std::array<float, TDistribution::MaxK::value*(TDistribution::MaxK::value-1)/2> computePearsonChiSquaredMergeMetric(const TDistribution& distribution)
    {
        const uint32_t numComponents = distribution.getK();
        const uint32_t numActiveKernels = (numComponents+TScalar::Width::value-1)/TScalar::Width::value;

        std::array<TKernel, TDistribution::NumKernels::value> selfProduct;
        for (uint32_t k=0; k<numActiveKernels; ++k)
        {
            selfProduct[k] = distribution.getComponent(k);
            selfProduct[k].product(selfProduct[k]);
        }

        std::array<TScalar, TDistribution::NumKernels::value*TDistribution::MaxK::value> componentToKernelSimilarityValues;

        for (uint32_t i=0; i<numComponents-1; ++i)
        {
            const size_t kernelIndexI = i/TScalar::Width::value;
            const size_t inKernelIndexI = i%TScalar::Width::value;

            const TKernel componentI = distribution.getComponent(kernelIndexI).extract(inKernelIndexI);
            const TKernel componentISqr = selfProduct[kernelIndexI].extract(inKernelIndexI);

            for (uint32_t j=(i+1)/TScalar::Width::value; j<numActiveKernels; ++j)
            {
                TKernel productIJ{componentI};
                productIJ.product(distribution.getComponent(j));

                const TKernel kernelJ{distribution.getComponent(j)};
                TKernel merged{componentI};

                {
                    TKernel kJForMerge{kernelJ};
                    for (uint32_t l=0; l<TScalar::Width::value; ++l)
                        merged.mergeComponent(l, l, kJForMerge);
                }

                const TScalar quotientISqrMerged = TKernel{componentISqr}.division(merged);
                const TScalar quotientIJMerged = TKernel{productIJ}.division(merged);
                const TScalar quotientJSqrMerged = TKernel{selfProduct[j]}.division(merged);

                const TScalar divergence = quotientISqrMerged+2.0f*quotientIJMerged+quotientJSqrMerged-merged.m_weights;

                componentToKernelSimilarityValues[i*TDistribution::NumKernels::value+j] = divergence;
            }
        }

        std::array<float, TDistribution::MaxK::value*(TDistribution::MaxK::value-1)/2> similarityValues;

        uint32_t index = 0;
        for (uint32_t i=0; i<numComponents-1; ++i)
            for (size_t j=i+1; j<numComponents; ++j, ++index)
                similarityValues[index] = componentToKernelSimilarityValues[i*TDistribution::NumKernels::value+j/TScalar::Width::value][j%TScalar::Width::value];

        return similarityValues;
    }

    template<typename TSampleData>
    FINLINE static void repositionSample(TSampleData& sample, const Point targetPosition)
    {
        if (std::isinf(sample.distance))
        {
            sample.position = targetPosition;
            return;
        }
        else if (EXPECT_NOT_TAKEN(!(sample.distance > 0.0f)))
        {
            SLog(EWarn, "invalid sample distance %f", sample.distance);
            return;
        }
        const Point lightSource = sample.position+sample.direction*sample.distance;
        const Vector targetToLightSource = lightSource-targetPosition;
        const float newDistance = targetToLightSource.length();

        sample.position  = targetPosition;
        sample.distance  = newDistance;
        sample.direction = targetToLightSource/newDistance;
    }

    FINLINE static std::array<typename TScalar::BooleanType, TDistribution::NumKernels::value> computeMaskFromComponentIndexVector(const std::vector<uint32_t> &activeComponentIds)
    {
        using TBool = typename TScalar::BooleanType;
        std::array<TBool, TDistribution::NumKernels::value> mask;
        std::fill(mask.begin(), mask.end(), TBool{false});

        for (uint32_t component : activeComponentIds)
            mask[component/TScalar::Width::value].insert(component%TScalar::Width::value, true);

        return mask;
    }
};

GUIDING_NAMESPACE_END

#endif // PATHGUIDINGMIXTURESTATSFACTORY_H
