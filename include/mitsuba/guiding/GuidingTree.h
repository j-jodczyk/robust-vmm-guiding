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

#include "BSPTree.h"

#include <mitsuba/core/properties.h>
#include <mitsuba/core/aabb.h>

#ifdef GUIDING_USE_PARALLEL
#include <omp.h>
#include "AtomicallyGrowingVector.h"
#include <immintrin.h>
#endif
#include <vector>
#include <numeric>
#include <functional>

#include "GuidingRegion.h"
#include "SampleStats.h"

GUIDING_NAMESPACE_BEGIN

template<typename _TGuidingRegion>
class GuidingTree {
private:
    //redefinition to aid autocompletion tools
    using TGuidingRegion = GuidingRegion<typename _TGuidingRegion::DistributionType, typename _TGuidingRegion::StatisticsType>;
    static_assert (std::is_same<TGuidingRegion, _TGuidingRegion>::value, "GuidingRegion is required as a wrapper of the guiding data.");
public:
    typedef TGuidingRegion GuidingRegionType;

    enum SplitType{
        kNodeBoundCenterXYZ = 0,            //useful for photons
        kNodeBoundCenterAndMaxExtend,       //useful for photons
        kSampleBoundCenterAndMaxExtend,     //useful for photons
        kSampleMeanXYZ,                     //not useful for photons
        kSampleMeanAndVariance              //not useful for photons
    };

    struct Settings{

        uint32_t maxSamples{2000U};
        uint32_t minSamples{100U};
        uint32_t maxDepth{32U};
        SplitType splitType{kNodeBoundCenterXYZ};
        float decayOnSpatialSplit {0.25f};

        Settings(const Properties &props = Properties()) {
            maxSamples = props.getSize("tree.maxSamples", maxSamples);
            minSamples = props.getSize("tree.minSamples", minSamples);
            maxDepth   = props.getSize("tree.maxDepth",   maxDepth);
            splitType  = static_cast<SplitType>(props.getInteger("tree.splitType", static_cast<int>(splitType)));
            decayOnSpatialSplit = props.getFloat("decayOnSpatialSplit", 0.25f);
        }

        Settings(Stream *stream)
            : maxSamples{stream->readUInt()}, minSamples{stream->readUInt()},
              maxDepth{stream->readUInt()}, splitType{static_cast<SplitType>(stream->readInt())},
              decayOnSpatialSplit{stream->readSingle()}
        {
        }

        void serialize(Stream* stream) const{
            stream->writeUInt(maxSamples);
            stream->writeUInt(minSamples);
            stream->writeUInt(maxDepth);
            stream->writeInt(static_cast<int>(splitType));
            stream->writeSingle(decayOnSpatialSplit);
        }

        std::string toString() const{
            std::ostringstream oss;
            oss << "BSPTreeSettings[" << '\n'
                << "    maxSamples = " << maxSamples << '\n'
                << "    minSamples = " << minSamples << '\n'
                << "    maxDepth   = " << maxDepth << '\n'
                << "    splitType  = " << static_cast<int>(splitType) << '\n'
                << "]"  << '\n';
            return oss.str();
        }
    };

    GuidingTree() = default;
    explicit GuidingTree(Settings settings)
        : m_settings{settings}
    {
    }

    explicit GuidingTree(Stream *stream)
        : m_tree{stream}, m_bounds{stream}, m_totalNumSamples{stream->readSize()}, m_maxDepth{stream->readUInt()}, m_settings{stream}
    {
        size_t numRegions = stream->readSize();
        m_regions.reserve(numRegions);
        for (size_t i=0; i<numRegions; ++i)
            m_regions.emplace_back(stream);
    }

    void serialize(Stream* stream) const{
        m_tree.serialize(stream);
        m_bounds.serialize(stream);

        stream->writeSize(m_totalNumSamples);
        stream->writeUInt(m_maxDepth);

        m_settings.serialize(stream);

        size_t numRegions = m_regions.size();
        stream->writeSize(numRegions);
        for (size_t i=0; i<numRegions; ++i)
            m_regions[i].serialize(stream);
    }

    void initTree(AABB bounds){
        m_totalNumSamples = 0;
        m_bounds = bounds;
        m_maxDepth = 0;
        m_regions.clear();
        m_regions.resize(1);
        m_regions[0].bound = bounds;

        m_tree.clear();
    }

    template<typename TSampleRange, typename TPointRange>
    void updateTree(TSampleRange& samples, TPointRange& zeroValuedSamples, size_t numThreads = 1, std::function<void(TGuidingRegion& region)> regionUpdate = nullptr)
    {
        size_t currentNumSamples = samples.size();
        m_totalNumSamples += currentNumSamples;

#ifdef GUIDING_USE_PARALLEL
        //estimate number of leaf nodes
        int numEstLeafs = m_regions.size()+(currentNumSamples*2)/m_settings.maxSamples+32;
        m_tree.reserve(4*numEstLeafs);
        m_regions.reserve(2*numEstLeafs);
#endif

        if (m_tree.empty())
        {
            for (const auto& sample : samples)
                m_regions[m_tree[0].childIdx].statsSinceLastSplit += sample.position;
        }

#ifdef GUIDING_USE_PARALLEL
#pragma omp parallel num_threads(numThreads)
#pragma omp single nowait
#endif
        updateTreeNode(0, 0, m_bounds, samples.begin(), samples, zeroValuedSamples, regionUpdate);
    }

    template<typename TSampleRange, typename TPointRange>
    void updateTreeNode(uint32_t nodeIdx, uint32_t depth, AABB bound,
                        typename TSampleRange::iterator globalBeginSamples,
                        TSampleRange& samples, TPointRange& zeroValuedSamples,
                        std::function<void(TGuidingRegion& region)>& regionUpdate)
    {
        const BSPTree::Node& node = m_tree[nodeIdx];

        if (node.isLeafNode())
        {
            //set the data in the leaf node or continue with splitting

            const uint32_t dataIdx = node.childIdx;
            TGuidingRegion& region = m_regions[dataIdx];

            if (depth < m_settings.maxDepth && region.statsSinceLastSplit.numSamples >= m_settings.maxSamples)
            {
                //decide split dimension and position
                std::pair<uint8_t, float> splitDimAndPoint = estimateSplitDimensionAndSplitPoint(depth, bound, region.statsSinceLastSplit);

                //initialize child nodes
                try {
                    const uint32_t idxLeft = m_tree.splitNode(nodeIdx, splitDimAndPoint.second, splitDimAndPoint.first);
                    const uint32_t idxRight = idxLeft+1;

                    //decay old sample statistics and add half of the new statistics to each child node
                    region.statsUntilLastSplit *= m_settings.decayOnSpatialSplit;
                    region.statsUntilLastSplit += region.statsSinceLastSplit*0.5f;
                    region.statsSinceLastSplit.clear();
                    region.spatialSplitFlag = true;

#ifdef GUIDING_USE_PARALLEL
                    const uint32_t dataIdxRight = std::distance(m_regions.begin(), m_regions.back_insert(region));
#else
                    const uint32_t dataIdxRight = m_regions.size();
                    m_regions.emplace_back(region);
#endif

                    m_tree.setDataIndex(idxRight, dataIdxRight);

                    if (depth+1 >= m_settings.maxDepth){
                        SLog(EWarn, "!!! BSPTree: tree depth reached max depth: %d at node %d containing region %s !!!", m_maxDepth, nodeIdx, region.toString().c_str());
                        m_maxDepth = m_settings.maxDepth;
                    }
                    else
                        m_maxDepth = depth+1;
                }
                catch (std::length_error& e)
                {
                    SLog(EWarn, "!!! BSPTree: tree ran out of preallocated nodes: %s !!!", e.what());
                }
            }

            if (node.isLeafNode())
            {
                //if node is not being split, store the final data and call the region update function

                region.bound = bound;
                region.dStart = std::distance(globalBeginSamples, samples.begin());
                region.nSamples = samples.size();

                if (regionUpdate)
                    regionUpdate(region);
            }
        }

        //processsed afterwards to simplify splitting
        if (node.isInnerNode())
        {
            //sort sample data as long as we are at an inner node

            const uint8_t  splitDimension = node.splitDim;
            const float    pivot          = node.splitPos;
            const uint32_t leftChildIdx   = node.childIdx;
            const uint32_t rightChildIdx  = leftChildIdx+1;

            const BSPTree::Node& leftNode  = m_tree[leftChildIdx];
            const BSPTree::Node& rightNode = m_tree[rightChildIdx];

            TGuidingRegion* leftRegion  = leftNode.isLeafNode()  ? &m_regions[m_tree[leftChildIdx].childIdx]  : nullptr;
            TGuidingRegion* rightRegion = rightNode.isLeafNode() ? &m_regions[m_tree[rightChildIdx].childIdx] : nullptr;

            AABB leftBound {bound}, rightBound {bound};
            leftBound.max[splitDimension] = pivot;
            rightBound.min[splitDimension] = pivot;

            typename TPointRange::iterator centerZeroValuedSamples = std::partition(zeroValuedSamples.begin(), zeroValuedSamples.end(), [splitDimension, pivot](const Point& p) -> bool
            {
                return p[splitDimension] < pivot;
            });

            typename TSampleRange::iterator centerSamples;

            if (leftNode.isLeafNode() || rightNode.isLeafNode())
            {
                SampleStats leftStats, rightStats;

                centerSamples = std::partition(samples.begin(), samples.end(), [splitDimension, pivot, &leftStats, &rightStats](const typename TSampleRange::value_type& sample) -> bool
                {
                    const bool left = sample.position[splitDimension] < pivot;
                    const bool splatted = (sample.flags&ESplatted);
                    if (!splatted)
                    {
                        if (left)
                            leftStats  += sample.position;
                        else
                            rightStats += sample.position;
                    }

                    return left;
                });
                //statistics are only computed on non-splatted samples, scale to the total number of samples
                //set some reasonable value for the mean if all samples are splatted out
                if (leftRegion)
                {
                    if (leftStats.numSamples > 0)
                    {
                        leftStats *= static_cast<Float>(std::distance(samples.begin(), centerSamples))/static_cast<Float>(leftStats.numSamples);
                    }
                    else
                    {
                        leftStats.numSamples = std::distance(samples.begin(), centerSamples);
                        leftStats.mean = leftRegion->statsSinceLastSplit.mean;
                        leftStats.unnormalizedVariance = leftRegion->statsSinceLastSplit.unnormalizedVariance*m_settings.decayOnSpatialSplit;
                        leftStats.unnormalizedVariance[splitDimension] *= 0.5f;
                        leftStats.sampleBound.expandBy(leftStats.mean);
                        leftStats.mean[splitDimension] = std::min(pivot, leftStats.mean[splitDimension]);
                    }
                    //add some small safety to the bounds
                    leftStats.sampleBound.min -= (leftStats.sampleBound.max-leftStats.sampleBound.min)*Epsilon+Vector{Epsilon};
                    leftStats.sampleBound.max += (leftStats.sampleBound.max-leftStats.sampleBound.min)*Epsilon+Vector{Epsilon};
                    leftRegion->statsSinceLastSplit += leftStats;
                    leftRegion->statsSinceLastSplit.numZeroValuedSamples += std::distance(zeroValuedSamples.begin(), centerZeroValuedSamples);
                }
                if (rightRegion)
                {
                    if (rightStats.numSamples > 0)
                    {
                        rightStats *= static_cast<Float>(std::distance(centerSamples, samples.end()))/static_cast<Float>(rightStats.numSamples);
                    }
                    else
                    {
                        rightStats.numSamples = std::distance(centerSamples, samples.end());
                        rightStats.mean = rightRegion->statsSinceLastSplit.mean;
                        rightStats.unnormalizedVariance = rightRegion->statsSinceLastSplit.unnormalizedVariance*m_settings.decayOnSpatialSplit;
                        rightStats.unnormalizedVariance[splitDimension] *= 0.5f;
                        rightStats.sampleBound.expandBy(rightStats.mean);
                        rightStats.mean[splitDimension] = std::max(pivot, rightStats.mean[splitDimension]);
                    }
                    //add some small safety to the bounds
                    rightStats.sampleBound.min -= (rightStats.sampleBound.max-rightStats.sampleBound.min)*Epsilon+Vector{Epsilon};
                    rightStats.sampleBound.max += (rightStats.sampleBound.max-rightStats.sampleBound.min)*Epsilon+Vector{Epsilon};
                    rightRegion->statsSinceLastSplit += rightStats;
                    rightRegion->statsSinceLastSplit.numZeroValuedSamples += std::distance(centerZeroValuedSamples,   zeroValuedSamples.end());
                }
            }
            else
            {
                centerSamples = std::partition(samples.begin(), samples.end(), [splitDimension, pivot](const typename TSampleRange::value_type& sample) -> bool
                {
                    return sample.position[splitDimension] < pivot;
                });
            }

            TSampleRange leftSamples{samples.begin(), centerSamples}, rightSamples{centerSamples, samples.end()};
            TPointRange  leftZeroValuedSamples{zeroValuedSamples.begin(), centerZeroValuedSamples}, rightZeroValuedSamples{centerZeroValuedSamples, zeroValuedSamples.end()};

#ifdef GUIDING_USE_PARALLEL
#pragma omp task mergeable
#endif
            updateTreeNode(rightChildIdx, depth+1, rightBound, globalBeginSamples, rightSamples, rightZeroValuedSamples, regionUpdate);
            updateTreeNode(leftChildIdx,  depth+1, leftBound,  globalBeginSamples, leftSamples,  leftZeroValuedSamples,  regionUpdate);
        }
    }

    GUIDING_INLINE std::pair<uint8_t, float> estimateSplitDimensionAndSplitPoint(const uint32_t depth, const AABB &bound, const SampleStats &splitStats) const {
        uint8_t splitDim = 0;
        float splitPoint = 0.0f;

        auto maxDimension = [](const Vector& v) -> uint8_t
        {
            return v[v[1] > v[0]] > v[2] ? v[1] > v[0] : 2;
        };

        switch(m_settings.splitType){
            case kNodeBoundCenterXYZ:
                splitDim = depth % 3;
                splitPoint = bound.getCenter()[splitDim];
                break;
            case kNodeBoundCenterAndMaxExtend:
                splitDim = maxDimension(bound.getExtents());
                splitPoint = bound.getCenter()[splitDim];
                break;
            case kSampleBoundCenterAndMaxExtend:
                splitDim = maxDimension(splitStats.sampleBound.getExtents());
                splitPoint = splitStats.sampleBound.getCenter()[splitDim];
                break;
            case kSampleMeanXYZ:
                splitDim = depth % 3;
                splitPoint = splitStats.mean[splitDim];
                break;
            case kSampleMeanAndVariance:
                //the variance does not need to be normalized to determine the dimension with the largest extent
                splitDim = maxDimension(splitStats.unnormalizedVariance);
                splitPoint = splitStats.mean[splitDim];
                break;
        }
        
        return std::make_pair(splitDim, splitPoint);
    }

    bool findRegionIndexAndComputeBounds(const Point position, uint32_t &regionIdx, AABB &bounds) const {
        if (EXPECT_NOT_TAKEN(!m_bounds.contains(position))){
            return false;
        }

        uint32_t rootIdx = 0;
        bounds = m_bounds;
        regionIdx = getRegionIndexAndComputeBounds(position,  rootIdx, bounds);

        Guiding_Assert(regionIdx < m_regions.size());
        Guiding_Assert(m_regions[regionIdx].bound.contains(position));

        return true;
    }

    uint32_t getRegionIndexAndComputeBounds(const Point position, uint32_t nodeIdx, AABB& bounds) const {
        std::pair<uint32_t, AABB> nodeIdxAndBounds = m_tree.getNodeIndexAndRefineBounds(position, nodeIdx, bounds);
        bounds = nodeIdxAndBounds.second;
        Guiding_Assert(bounds.contains(position));
        return m_tree.getDataIndex(nodeIdxAndBounds.first);
    }

    uint32_t getRegionIndex(const Point position) const {
        Guiding_Assert(m_bounds.contains(position));
        return m_tree.findDataIndex(position);
    }

    const BSPTree& getTree() const
    {
        return m_tree;
    }

    size_t getNumRegions() const
    {
        return m_regions.size();
    }

    TGuidingRegion& getRegion(const size_t index)
    {
        return m_regions[index];
    }

    const TGuidingRegion& getRegion(const size_t index) const
    {
        return m_regions[index];
    }

    TGuidingRegion& getRegion(const Point position)
    {
        return m_regions[getRegionIndex(position)];
    }

    const TGuidingRegion& getRegion(const Point position) const
    {
        return m_regions[getRegionIndex(position)];
    }

    AABB getBounds() const {
        return m_bounds;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "GuidingTree[" << '\n'
            << "  totalNumSamples = " << m_totalNumSamples << ",\n"
            << "  bound = " << m_bounds.toString() << ",\n"
            << "  tree = " << m_tree.toString() << ",\n"
            << "  nData = " << m_regions.size() << ",\n"
            << "  maxDepth = " << m_maxDepth << ",\n"
            << "  settings = " << m_settings.toString() << ",\n"
            << "]";
        return oss.str();
    }

    struct Statistics{
        uint32_t numNodes {0};
        uint32_t reservedNodes {0};

        uint32_t numRegions {0};
        uint32_t reservedRegions {0};

        uint32_t minNumSamples {0};
        uint32_t maxNumSamples {0};
        float    avgNumSamples {0};
        float    stdNumSamples {0};

        uint32_t numEmptyRegions {0};

        uint32_t minNumTrainingIterations {0};
        uint32_t maxNumTrainingIterations {0};
        float    avgNumTrainingIterations {0.0f};
        float    stdNumTrainingIterations {0.0f};

        uint32_t minNumComponents {0};
        uint32_t maxNumComponents {0};
        float    avgNumComponents {0.0f};
        float    stdNumComponents {0.0f};

        uint32_t numRegionsWithMaxComponents {0};

        uint32_t numTrainingSamples {0};

        template<typename TSampleData>
        std::string toString() const {
            std::ostringstream oss;
            oss << "GuidingTree::Statistics[" << '\n'
                << "-- Nodes --" << '\n'
                << "   numNodes                         = " << numNodes << '\n'
                << "   sizePerNode                      = " << sizeof(BSPTree::Node) << " B\n"
                << "   usedNodeStorage                  = " << memString(numNodes*sizeof(BSPTree::Node)) << '\n'
                << "   reservedNodeStorage              = " << memString(reservedNodes*sizeof(BSPTree::Node)) << '\n'
                << "-- Regions --" << '\n'
                << "   numRegions                       = " << numRegions << '\n'
                << "   sizePerRegion                    = " << sizeof(TGuidingRegion) << " B\n"
                << "   usedRegionStorage                = " << memString(numRegions*sizeof(TGuidingRegion)) << '\n'
                << "   reservedRegionStorage            = " << memString(reservedRegions*sizeof(TGuidingRegion)) << '\n'
                << "   minNumSamples                    = " << minNumSamples << '\n'
                << "   maxNumSamples                    = " << maxNumSamples << '\n'
                << "   avgNumSamples                    = " << avgNumSamples << '\n'
                << "   stdNumSamples                    = " << stdNumSamples << '\n'
                << "   numEmptyRegions                  = " << numEmptyRegions << '\n'
                << "-- EM-Fit --" << '\n'
                << "   minNumTrainingIterations         = " << minNumTrainingIterations << '\n'
                << "   maxNumTrainingIterations         = " << maxNumTrainingIterations << '\n'
                << "   avgNumTrainingIterations         = " << avgNumTrainingIterations << '\n'
                << "   stdNumTrainingIterations         = " << stdNumTrainingIterations << '\n'
                << "-- Mixtures --" << '\n'
                << "   minNumComponents                 = " << minNumComponents << '\n'
                << "   maxNumComponents                 = " << maxNumComponents << '\n'
                << "   avgNumComponents                 = " << avgNumComponents << '\n'
                << "   stdNumComponents                 = " << stdNumComponents << '\n'
                << "   numRegionsWithMaxComponents      = " << numRegionsWithMaxComponents << '\n'
                << "-- TrainingSamples --" << '\n'
                << "   numTrainingSamples               = " << numTrainingSamples << '\n'
                << "   sizePerTrainingSample            = " << sizeof(TSampleData) << " B\n"
                << "   usedTrainingSampleStorage        = " << memString(numTrainingSamples*sizeof(TSampleData)) << '\n'
                << "]";
            return oss.str();
        }

        Statistics(const GuidingTree& tree){
            numNodes = tree.m_tree.size();
            reservedNodes = tree.m_tree.capacity();

            numRegions = tree.m_regions.size();
            reservedRegions = tree.m_regions.capacity();

            minNumSamples = std::numeric_limits<uint32_t>::max();
            maxNumSamples = std::numeric_limits<uint32_t>::min();
            avgNumSamples = 0.0f;
            stdNumSamples = 0.0f;

            minNumTrainingIterations = std::numeric_limits<uint32_t>::max();
            maxNumTrainingIterations = std::numeric_limits<uint32_t>::min();
            avgNumTrainingIterations = 0.0f;
            stdNumTrainingIterations = 0.0f;

            minNumComponents = std::numeric_limits<uint32_t>::max();
            maxNumComponents = std::numeric_limits<uint32_t>::min();
            avgNumComponents = 0.0f;
            stdNumComponents = 0.0f;

            size_t numValidRegions = 0;

#ifdef GUIDING_USE_PARALLEL
#pragma omp parallel for reduction(min: minNumSamples, minNumTrainingIterations, minNumComponents) \
                         reduction(max: maxNumSamples, maxNumTrainingIterations, maxNumComponents) \
                         reduction(+:   avgNumSamples, avgNumTrainingIterations, avgNumComponents, \
                                        stdNumSamples, stdNumTrainingIterations, stdNumComponents, \
                                        numValidRegions, numEmptyRegions, numRegionsWithMaxComponents)
#endif
            for (uint32_t i = 0; i < tree.m_regions.size(); i++){
                const uint32_t nSamples = tree.m_regions[i].nSamples;

                minNumSamples = std::min(minNumSamples, nSamples);
                maxNumSamples = std::max(maxNumSamples, nSamples);
                avgNumSamples += static_cast<float>(nSamples);
                stdNumSamples += static_cast<float>(nSamples)*static_cast<float>(nSamples);

                const bool valid = tree.m_regions[i].valid;
                const bool empty = (nSamples == 0);

                numValidRegions += valid;
                numEmptyRegions += empty;

                const uint32_t numTrainingIterations = valid ? tree.m_regions[i].distribution.m_numEMIterations : 0;

                minNumTrainingIterations = std::min(minNumTrainingIterations, numTrainingIterations);
                maxNumTrainingIterations = std::max(maxNumTrainingIterations, numTrainingIterations);
                avgNumTrainingIterations += static_cast<float>(numTrainingIterations);
                stdNumTrainingIterations += static_cast<float>(numTrainingIterations)*static_cast<float>(numTrainingIterations);

                const uint32_t numComponents = valid ? tree.m_regions[i].distribution.m_K : 0;

                minNumComponents = std::min(minNumComponents, numComponents);
                maxNumComponents = std::max(maxNumComponents, numComponents);
                avgNumComponents += static_cast<float>(numComponents);
                stdNumComponents += static_cast<float>(numComponents)*static_cast<float>(numComponents);

                const bool hasMaxComponents = (numComponents == TGuidingRegion::DistributionType::MaxK::value);
                numRegionsWithMaxComponents += hasMaxComponents;
            }

            numTrainingSamples = avgNumSamples;

            const float invNumRegions = 1.0f/static_cast<float>(tree.m_regions.size());
            avgNumSamples *= invNumRegions;
            stdNumSamples = math::safe_sqrt(stdNumSamples*invNumRegions-avgNumSamples*avgNumSamples);

            if (numValidRegions > 0)
            {
                const float invNumValidRegions = 1.0f/static_cast<float>(numValidRegions);

                avgNumTrainingIterations *= invNumValidRegions;
                stdNumTrainingIterations = math::safe_sqrt(stdNumTrainingIterations*invNumValidRegions-avgNumTrainingIterations*avgNumTrainingIterations);

                avgNumComponents *= invNumValidRegions;
                stdNumComponents = math::safe_sqrt(stdNumComponents*invNumValidRegions-avgNumComponents*avgNumComponents);
            }
        }
    };

private:
    BSPTree m_tree;
#ifdef GUIDING_USE_PARALLEL
    AtomicallyGrowingVector<TGuidingRegion> m_regions;
#else
    std::vector<TGuidingRegion> m_regions;
#endif

    AABB m_bounds;

    size_t m_totalNumSamples {0};
    uint32_t m_maxDepth {0};

public:
    Settings m_settings;
};


GUIDING_NAMESPACE_END
