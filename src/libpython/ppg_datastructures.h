/*
    This file is part of the implementation of the SIGGRAPH 2020 paper
    "Robust Fitting of Parallax-Aware Mixtures for Path Guiding".
    The implementation extends Mitsuba, a physically based rendering system.

    This file contains fragments of the "Practical Path Guiding" code from
    https://github.com/Tom94/practical-path-guiding.git
    It is used solely for visualization of their guiding distributions for comparison.

    Copyright (c) 2020 Lukas Ruppert, Sebastian Herholz.
    Copyright (c) 2017 by ETH Zurich, Thomas Mueller.

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

#ifndef PPG_DATASTRUCTURES_H
#define PPG_DATASTRUCTURES_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/matrix.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/core/logger.h>

#include <string>
#include <iostream>
#include <fstream>
#include <type_traits>
#include <array>
#include <vector>
#include <ostream>

MTS_NAMESPACE_BEGIN

namespace PPG
{
    class BlobReader {
    public:
        BlobReader(const std::string& filename) : f(filename, std::ios::in | std::ios::binary) {}

        template <typename Type>
        typename std::enable_if<std::is_standard_layout<Type>::value, BlobReader&>::type
            operator >> (Type& Element) {
            Read(&Element, 1);
            return *this;
        }

        // CAUTION: This function may break down on big-endian architectures.
        //          The ordering of bytes has to be reverted then.
        template <typename T>
        void Read(T* Dest, size_t Size) {
            f.read(reinterpret_cast<char*>(Dest), Size * sizeof(T));
        }

        bool isValid() const {
            return (bool)(f);
        }

    private:
        std::ifstream f;
    };

    static Vector2f dirToCanonical(const Vector& d) {
        if (!std::isfinite(d.x) || !std::isfinite(d.y) || !std::isfinite(d.z)) {
            return {0, 0};
        }

        const Float cosTheta = std::min(std::max(d.z, -1.0f), 1.0f);
        Float phi = std::atan2(d.y, d.x);
        while (phi < 0)
            phi += 2.0 * M_PI;

        return {(cosTheta + 1) / 2, phi / (2 * M_PI)};
    }

    struct QuadTreeNode {
        std::array<float, 4> data;
        std::array<uint16_t, 4> children;

        inline bool isLeaf(int index) const {
            return children[index] == 0;
        }

        int computeDepth(const std::vector<QuadTreeNode>& nodes) const {
            int maxDepth = 0;
            for (int i = 0; i < 4; ++i) {
                if (!isLeaf(i)) {
                    maxDepth = std::max(maxDepth, nodes[children[i]].computeDepth(nodes) + 1);
                }
            }

            return maxDepth;
        }

        float computeMax(const std::vector<QuadTreeNode>& nodes) const {
            float maximum = 0;
            for (int i = 0; i < 4; ++i) {
                if (!isLeaf(i)) {
                    maximum = std::max(maximum, nodes[children[i]].computeMax(nodes));
                } else {
                    maximum = std::max(maximum, data[i]);
                }
            }

            return 4 * maximum;
        }

        int childIndex(Vector2f& p) const {
            int res = 0;
            for (int i = 0; i < Vector2f::dim; ++i) {
                if (p[i] < 0.5f) {
                    p[i] *= 2;
                } else {
                    p[i] = (p[i] - 0.5f) * 2;
                    res |= 1 << i;
                }
            }

            return res;
        }

        Float pdf(Vector2f& p, const std::vector<QuadTreeNode>& nodes) const {
            SAssert(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1);
            const int c = childIndex(p);
            if (!(data[c] > 0)) {
                return 0;
            }

            const Float factor = 4 * data[c] / (data[0] + data[1] + data[2] + data[3]);
            if (isLeaf(c)) {
                return factor;
            } else {
                return factor * nodes[children[c]].pdf(p, nodes);
            }
        }
    };

    class DTree {
    public:
        const Point3f pos() const { return mAABB.getCenter(); }

        const AABB& aabb() const { return mAABB; }

        float mean() const { return mMean; }

        float step() const {
            int dim = 1 << mDepth;
            return 1.0f / dim;
        }

        size_t nSamples() const {
            return mNumSamples;
        }

        bool read(BlobReader& blob) {
            uint64_t numNodes;
            uint64_t numSamples;
            blob >> mAABB.min.x >> mAABB.min.y >> mAABB.min.z >> mAABB.max.x >> mAABB.max.y >> mAABB.max.z >> mMean >> numSamples >> numNodes;
            if (!blob.isValid()) {
                return false;
            }

            mNumSamples = (size_t)numSamples;

            if (!std::isfinite(mMean)) {
                cerr << "INVALID MEAN: " << mMean << endl;
            }

            mAABB.max += mAABB.min;

            mNodes.resize(numNodes);
            for (size_t i = 0; i < mNodes.size(); ++i) {
                auto& n = mNodes[i];
                for (size_t j = 0; j < 4; ++j) {
                    blob >> n.data[j];
                    if (!std::isfinite(n.data[j])) {
                        cerr << "INVALID NODE: " << n.data[j] << endl;
                    }
                    blob >> n.children[j];
                }
            }

            mDepth = computeDepth();
            mMax = computeMax();

            return true;
        }

        Float pdfCanonical(Vector2f p) const {
            if (!(mean() > 0)) {
                return 1 / (4 * M_PI);
            }

            return mNodes[0].pdf(p, mNodes) / (4 * M_PI);
        }

        Float pdf(Vector dir) const {
            return pdfCanonical(dirToCanonical(dir));
        }

    private:
        int computeDepth() const {
            return mNodes[0].computeDepth(mNodes) + 1;
        }

        float computeMax() const {
            if (mNumSamples == 0) {
                return 0;
            }
            const float factor = 1 / (float)(4 * M_PI * mNumSamples);
            return factor * mNodes[0].computeMax(mNodes);
        }

        int depth() const { return mDepth; }

        AABB mAABB;

        std::vector<QuadTreeNode> mNodes;
        float mMean;
        size_t mNumSamples;

        int mDepth;
        float mMax;
    };

    struct STree {
        std::vector<DTree> dTrees;
        AABB aabb;

        const DTree& findDTree(const Point& pos) {
            //brute-force search due to missing spatial tree structure
            for (const DTree& dTree : dTrees)
            {
                if (dTree.aabb().contains(pos))
                {
                    return dTree;
                }
            }
            SLog(EError, "Point %s is not contained in any DTree's bounds.", pos.toString().c_str());
            return dTrees.at(-1);
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "STree ["
                << "\n  " << aabb.toString()
                << "\n  num DTrees: " << dTrees.size()
                << "\n]";
            return oss.str();
        }
    };
}

MTS_NAMESPACE_END

#endif
