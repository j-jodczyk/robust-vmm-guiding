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

#ifndef BSPTREE_H
#define BSPTREE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
#include "guiding.h"
#ifdef GUIDING_USE_PARALLEL
#include "AtomicallyGrowingVector.h"
#endif

GUIDING_NAMESPACE_BEGIN

/**
 * @brief The BSPTree class implements a simple 3-dimensional kd-tree.
 * It uses only 8 bytes per node and supports atomic refinement.
 * Each spatial region is mapped to a 30-bit integer,
 * which can be used externally to index corresponding data.
 */
class BSPTree {
public:
    struct Node {
        Node() = default;

        explicit Node(Stream *stream){
            splitPos = stream->readFloat();
            const uint32_t splitDimAndChildIdx = stream->readUInt();
            splitDim = splitDimAndChildIdx >> 30;
            childIdx = splitDimAndChildIdx & 0x3FFFFFFFU;
        }

        void serialize(Stream* stream) const {
            stream->writeFloat(splitPos);
            const uint32_t splitDimAndChildIdx = (static_cast<uint32_t>(splitDim)<<30) | childIdx;
            stream->writeUInt(splitDimAndChildIdx);
        }

        FINLINE void setInnerNode(float splitPos, uint8_t splitDim, uint32_t childIdx) {
            Guiding_Assert((childIdx & ~0x3FFFFFFFU) == 0);
            this->splitPos = splitPos;
            this->splitDim = splitDim;
            this->childIdx = childIdx;
        }

        FINLINE void setLeafNode(uint32_t regionIdx) {
            Guiding_Assert((regionIdx & ~0x3FFFFFFFU) == 0);
            splitDim = 3;
            childIdx = regionIdx;
        }

        FINLINE bool isLeafNode() const {
            return splitDim == 3;
        }

        FINLINE bool isInnerNode() const {
            return !isLeafNode();
        }

        float    splitPos;
        uint8_t  splitDim :  2;
        uint32_t childIdx : 30;
    };

    BSPTree()
    {
        clear();
    }

    explicit BSPTree(Stream *stream)
    {
        size_t nNodes = stream->readSize();
        m_nodes.reserve(nNodes);
        for (size_t i=0; i<nNodes; ++i)
            m_nodes.emplace_back(stream);
    }

    void serialize(Stream* stream) const{
        size_t nNodes = m_nodes.size();
        stream->writeSize(nNodes);
        for (size_t i=0; i<nNodes; ++i)
            m_nodes[i].serialize(stream);
    }

    void clear()
    {
        m_nodes.clear();
        m_nodes.resize(1);
        m_nodes[0].setLeafNode(0);
    }

    FINLINE void reserve(size_t capacity)
    {
        m_nodes.reserve(capacity);
    }

    FINLINE bool empty() const
    {
        return m_nodes[0].isLeafNode();
    }

    FINLINE size_t size() const
    {
        return m_nodes.size();
    }

    FINLINE size_t capacity() const
    {
        return m_nodes.capacity();
    }

    std::pair<uint32_t, AABB> getNodeIndexAndRefineBounds(const Point position, uint32_t startingNodeIdx, AABB initialBounds) const {
        std::pair<uint32_t, AABB> nodeAndBounds;
        uint32_t& nodeIdx = nodeAndBounds.first;
        nodeIdx = startingNodeIdx;
        AABB& bounds = nodeAndBounds.second;
        bounds = initialBounds;

        while (!m_nodes[nodeIdx].isLeafNode())
        {
            uint32_t leftIdx = m_nodes[nodeIdx].childIdx;
            uint8_t splitDim = m_nodes[nodeIdx].splitDim;
            float pivot = m_nodes[nodeIdx].splitPos;

            if (position[splitDim] < pivot)
            {
                bounds.max[splitDim] = pivot;
                nodeIdx              = leftIdx;
            }
            else
            {
                bounds.min[splitDim] = pivot;
                nodeIdx              = leftIdx+1;
            }
        }

        return nodeAndBounds;
    }

    uint32_t getNodeIndex(const Point position) const {
        uint32_t nodeIdx = 0;

        while (!m_nodes[nodeIdx].isLeafNode())
            nodeIdx = m_nodes[nodeIdx].childIdx+(position[m_nodes[nodeIdx].splitDim] >= m_nodes[nodeIdx].splitPos);

        return nodeIdx;
    }

    FINLINE const Node& operator[](uint32_t nodeIdx) const
    {
        return m_nodes[nodeIdx];
    }

    FINLINE uint32_t getDataIndex(uint32_t nodeIdx) const
    {
        Guiding_Assert(m_nodes[nodeIdx].isLeafNode());
        return m_nodes[nodeIdx].childIdx;
    }

    FINLINE uint32_t findDataIndex(const Point position) const
    {
        return getDataIndex(getNodeIndex(position));
    }

    FINLINE void setDataIndex(uint32_t nodeIdx, uint32_t dataIdx)
    {
        Guiding_Assert(m_nodes[nodeIdx].isLeafNode());
        m_nodes[nodeIdx].childIdx = dataIdx;
    }

    /**
     * @brief splitNode
     * split leaf node at index \ref nodeIdx, creating two new child nodes.
     * previous data index is propagated to both child nodes.
     * @param nodeIdx index of the node to be split
     * @param splitPos spatial split position
     * @param splitDim spatial split dimension
     * @return index of the left child node. right child node is placed at the following index.
     */
    uint32_t splitNode(uint32_t nodeIdx, float splitPos, uint8_t splitDim)
    {
        Guiding_Assert(m_nodes[nodeIdx].isLeafNode());

        const uint32_t dataIndex = getDataIndex(nodeIdx);
#ifdef GUIDING_USE_PARALLEL
        const uint32_t childNodesIndex = std::distance(m_nodes.begin(), m_nodes.back_insert(2, Node()));
#else
        const uint32_t childNodesIndex = m_nodes.size();
        m_nodes.insert(m_nodes.end(), 2, Node());
#endif
        m_nodes[nodeIdx].setInnerNode(splitPos, splitDim, childNodesIndex);

        m_nodes[childNodesIndex].setLeafNode(dataIndex);
        m_nodes[childNodesIndex+1].setLeafNode(dataIndex);

        return childNodesIndex;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "BSPTree[" << '\n'
            << "  numNodes = " << m_nodes.size() << ",\n"
            << "]";
        return oss.str();
    }

private:
#ifdef GUIDING_USE_PARALLEL
    AtomicallyGrowingVector<Node> m_nodes;
#else
    std::vector<Node> m_nodes;
#endif
};

GUIDING_NAMESPACE_END

#endif
