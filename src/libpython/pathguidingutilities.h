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

#ifndef PATHGUIDINGUTILITIES_H
#define PATHGUIDINGUTILITIES_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/emitter.h>

#include <mitsuba/guiding/guiding.h>
#include <mitsuba/guiding/guiding_types.h>

#include "ppg_datastructures.h"

#include <vector>
#include <array>

GUIDING_NAMESPACE_BEGIN

/**
 * utility class for the Python API
 * this class is not needed for rendering in any way.
 */
class PathGuidingUtilities
{
public:
    static Vector2i renderSize;

    static inline Point2i dirVectorToSphericalPlanePoint(const Vector dir)
    {
        const float theta = acos(dir.z);
        const float phi   = atan2(dir.y, dir.x);

        return Point2i{std::min(std::max(static_cast<int>((M_PI-phi)*static_cast<float>(renderSize.x)*INV_TWOPI), 0), renderSize.x-1),
                       std::min(std::max(static_cast<int>(    theta *static_cast<float>(renderSize.y)*INV_PI),    0), renderSize.y-1)};
    }

    static inline Vector sphericalPlanePointToDirVector(const Point2i p)
    {
        const float theta =      (static_cast<float>(p.y+0.5f)/static_cast<float>(renderSize.y))*M_PI;
        const float phi   = M_PI-(static_cast<float>(p.x+0.5f)/static_cast<float>(renderSize.x))*2.0f*M_PI;

        return Vector{static_cast<float>(sin(theta)*cos(phi)), static_cast<float>(sin(theta)*sin(phi)), static_cast<float>(cos(theta))};
    }

    //envmap utilities
    static ref<Emitter> loadEnvmapForSampling(const std::string& filename);
    static std::vector<PathGuidingSampleData> sampleEnvMap(const Emitter* envMap, unsigned int numSamples, bool useImportanceSampling, const Point p=Point(0.0f), const Point2i seed=Point2i(0.0f));
    static ref<Bitmap> renderEnvmapPDF(const Emitter* envMap);

    //scene intersection based on pixel
    static Intersection getFirstSmoothSurfaceInteraction(const Scene* scene, Point2 pixel);

    //ground truth data generation
    static ref<Bitmap> renderSphericalView(const Scene* scene, const Point p, const std::string &integratorName, int numSamples);
    static ref<Bitmap> renderSphericalViewDistance(const Scene* scene, const Point p, int numSamples);
    static ref<Bitmap> renderBSDF(const Intersection& its);
    static ref<Bitmap> normalizeSphericalView(const ref<Bitmap> bitmap);

    //exporting samples
    static void exportSamplesAsCSV(const std::string& filename, const std::vector<PathGuidingSampleData>& samples);
    static void exportSamplesAsOBJ(const std::string& filename, const std::vector<PathGuidingSampleData>& samples);
    static ref<Bitmap> renderSampleBitmap(const std::vector<PathGuidingSampleData>& samples);
    static ref<Bitmap> renderSamplePDFBitmap(const std::vector<PathGuidingSampleData>& samples);
    static ref<Bitmap> renderSampleDistanceBitmap(const std::vector<PathGuidingSampleData>& samples);

    //importing samples
    static std::vector<PathGuidingSampleData> loadSamples(const std::string& filename);
    static std::vector<PathGuidingSampleData> loadSampleRange(const std::string& filename, size_t startIndex, size_t numSamples);

    //visualizing soft assignment
    static std::vector<PathGuidingSampleData> reweightSamplesWithComponentSoftAssignment(const VMFMixture& vmm, size_t component, const std::vector<PathGuidingSampleData>& samples);
    static std::vector<PathGuidingSampleData> reweightSamplesWithRelativePDF(const VMFMixture& vmm, const std::vector<PathGuidingSampleData>& samples);

    //exporting mixtures
    static void exportVMMAsCSV(const std::string& filename, const VMFMixture& vmm);
    static ref<Bitmap> renderVMMPDF(const VMFMixture& vmm);
    static ref<Bitmap> renderVMMDistance(const VMFMixture& vmm, const std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>& distances);
    static ref<Bitmap> renderVMFPDF(const VMFMixture& vmm, size_t component);

    //colored mixture plot
    static ref<Bitmap> renderVMMPDFColored(const VMFMixture& vmm, const std::array<Spectrum, VMFMixture::MaxK::value>& colors);

    //mixture likelihood and divergence
    static float computeLogLikelihood(const VMFMixture &distribution, const std::vector<PathGuidingSampleData>& samples);
    static float computePearsonChiSquaredDivergence(const VMFMixture &distribution, const std::vector<PathGuidingSampleData>& samples);

    //exporting the guiding field
    static void exportGuidingTreeASCII(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename);
    static void exportGuidingTreeOBJ(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename);
    static void exportGuidingTreeSampleBoundsOBJ(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename);
    static void exportGuidingTreeCentersOBJ(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename);

    //PPG
    static PPG::STree loadPPGSDTree(const std::string& filename);
    static ref<Bitmap> renderDTreePDF(const PPG::DTree& dTree);
    static void exportSTreeOBJ(const PPG::STree& sTree, const std::string &filename);
};

GUIDING_NAMESPACE_END

#endif // PATHGUIDINGUTILITIES_H
