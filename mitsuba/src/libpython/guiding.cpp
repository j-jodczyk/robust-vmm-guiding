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

#include "base.h"

#include <mitsuba/render/scene.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/renderjob.h>

#include <boost/algorithm/string.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
//#include <boost/python/numeric.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/implicit.hpp>

#include "pathguidingutilities.h"
#include <mitsuba/guiding/incrementalcovariance2d.h>
#include <mitsuba/guiding/incrementaldistance.h>
#include <mitsuba/guiding/incrementalpearsonchisquared.h>

using namespace mitsuba;
using namespace mitsuba::Guiding;

static std::array<typename VMFMixture::KernelType::ScalarType::BooleanType, VMFMixture::NumKernels::value> computeMaskFromComponentIndexVector(bp::list componentMask)
{
    std::vector<uint32_t> componentMaskVector;
    componentMaskVector.reserve(bp::len(componentMask));

    for (int i=0; i<bp::len(componentMask); ++i)
    {
        bp::extract<size_t> component(componentMask[i]);
        if (component.check())
            componentMaskVector.push_back(component());
        else
            SLog(EWarn, "invalid value encountered during component mask conversion");
    }

    return Guiding::GuidingStatsFactory<VMFMixture, PathGuidingMixtureStats<VMFMixture>>::computeMaskFromComponentIndexVector(componentMaskVector);
}

static bp::list pathGuidingMixtureStatsComputePearsonChiSquaredMergeMetric(VMFMixture& vmm)
{
    std::array<float, VMFMixture::MaxK::value*(VMFMixture::MaxK::value-1)/2> metric = GuidingFieldType::GuidingFieldFactoryType::StatisticsFactoryType::computePearsonChiSquaredMergeMetric(vmm);

    const uint32_t numComponents = vmm.getK();

    bp::list pythonList;
    for (size_t i=0; i<numComponents*(numComponents-1)/2; ++i)
        pythonList.append(metric[i]);

    return pythonList;
}

static bp::list incrementalPearsonChiSquaredComputeDivergence(const IncrementalPearsonChiSquared<VMFMixture>& incrementalPearsonChiSquared, VMFMixture& vmm)
{
    std::array<typename VMFMixture::ScalarType, VMFMixture::NumKernels::value> divergence = incrementalPearsonChiSquared.computeDivergence(vmm);

    const uint32_t numComponents = vmm.getK();

    bp::list pythonList;
    for (size_t i=0; i<numComponents; ++i)
        pythonList.append(divergence[i/VMFMixture::ScalarType::Width::value][i%VMFMixture::ScalarType::Width::value]);

    return pythonList;
}

static mitsuba::ref<Bitmap> PathGuidingUtilitiesRenderVMMPDFColored(const VMFMixture& vmm, bp::list pythonColors)
{
    std::array<Spectrum, VMFMixture::MaxK::value> colors;
    unsigned int count = 0;
    std::generate(colors.begin(), colors.end(), [pythonColors, &count](){
        Spectrum s;
        bp::extract<bp::list> colorList(pythonColors[count++]);
        if (colorList.check())
        {
            for (int i=0;i<Spectrum::dim && i<bp::len(colorList); ++i)
            {
                bp::extract<Float> value(colorList()[i]);
                if (value.check())
                    s[i] = value();
            }
        }
        return s;
    });

    return PathGuidingUtilities::renderVMMPDFColored(vmm, colors);
}

FINLINE static bool vmmRangeCheck(const VMFMixture& vmm, const size_t component)
{
    if (component < vmm.getK())
        return true;

    SLog(EWarn, "cannot access component %zu of mixture containing only %zu components", component, vmm.getK());
    return false;
}

static float vmfMixtureGetKappa(const VMFMixture& vmm, size_t component)
{
    if (vmmRangeCheck(vmm, component))
        return vmm.getComponent(component/VMFMixture::KernelType::ScalarType::Width::value).getKappa(component%VMFMixture::KernelType::ScalarType::Width::value);

    return -1.0f;
}

///range checked access to mixture parameters that also abstracts away the vectorization
static void vmfMixtureSetKappa(VMFMixture& vmm, size_t component, float kappa)
{
    if (vmmRangeCheck(vmm, component))
        vmm.getComponent(component/VMFMixture::KernelType::ScalarType::Width::value).setKappa(component%VMFMixture::KernelType::ScalarType::Width::value, kappa);
}

static Vector vmfMixtureGetMu(const VMFMixture& vmm, size_t component)
{
    if (vmmRangeCheck(vmm, component))
        return vmm.getComponent(component/VMFMixture::KernelType::ScalarType::Width::value).getMu(component%VMFMixture::KernelType::ScalarType::Width::value);

    return Vector{-1.0f};
}

static void vmfMixtureSetMu(VMFMixture& vmm, size_t component, Vector mu)
{
    if (vmmRangeCheck(vmm, component))
        vmm.getComponent(component/VMFMixture::KernelType::ScalarType::Width::value).setMu(component%VMFMixture::KernelType::ScalarType::Width::value, mu);
}

static float vmfMixtureGetWeight(const VMFMixture& vmm, size_t component)
{
    if (vmmRangeCheck(vmm, component))
        return vmm.getComponent(component/VMFMixture::KernelType::ScalarType::Width::value).getWeight(component%VMFMixture::KernelType::ScalarType::Width::value);

    return -1.0f;
}

static void vmfMixtureSetWeight(VMFMixture& vmm, size_t component, float weight)
{
    if (vmmRangeCheck(vmm, component))
        vmm.getComponent(component/VMFMixture::KernelType::ScalarType::Width::value).setWeight(component%VMFMixture::KernelType::ScalarType::Width::value, weight);
}

static float vmfMixturePDFComponent(VMFMixture& vmm, size_t component, Vector direction)
{
    if (vmmRangeCheck(vmm, component))
        return vmm.getComponent(component/VMFMixture::KernelType::ScalarType::Width::value).pdf(lightpmm::Vec<VMFMixture::KernelType::ScalarType>(direction))[component%VMFMixture::KernelType::ScalarType::Width::value];

    return -1.0f;
}

static float ScalarArray_getItem(const std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>& a, size_t index)
{
    if (index < VMFMixture::MaxK::value)
        return a[index/VMFMixture::KernelType::ScalarType::Width::value][index%VMFMixture::KernelType::ScalarType::Width::value];

    SLog(EWarn, "cannot access component %zu of array containing only %zu components", index, VMFMixture::NumKernels::value);
    return -1.0f;
}

static void ScalarArray_setItem(std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>& a, size_t index, float value)
{
    if (index < VMFMixture::MaxK::value)
        a[index/VMFMixture::KernelType::ScalarType::Width::value].insert(index%VMFMixture::KernelType::ScalarType::Width::value, value);
    else
        SLog(EWarn, "cannot access component %zu of array containing only %zu components", index, VMFMixture::NumKernels::value);
}

static std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>* ScalarArray_fromList(bp::list list)
{
    std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value> array;
    std::fill(array.begin(), array.end(), VMFMixture::KernelType::ScalarType{0.0f});

    const size_t numElements = bp::len(list);

    if (numElements > VMFMixture::MaxK::value)
        SLog(EWarn, "cannot insert %zu elements into array of length %zu. cutting off excess values", numElements, VMFMixture::NumKernels::value);

    for (size_t i=0; i < numElements && i < VMFMixture::MaxK::value; ++i)
    {
        bp::extract<Float> value(list[i]);
        if (value.check())
        {
            array[i/VMFMixture::ScalarType::Width::value].insert(i%VMFMixture::ScalarType::Width::value, value());
        }
        else
        {
            SLog(EWarn, "skipping invalid value at index %zu during array construction", i);
        }
    }

    return new std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>{array};
}

static std::string ScalarArray_repr(const std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>& a)
{
    std::ostringstream oss;
    oss << '[';

    for (size_t i=0; i < VMFMixture::MaxK::value; ++i)
    {
        if (i > 0)
            oss << ", ";
        oss << a[i/VMFMixture::KernelType::ScalarType::Width::value][i%VMFMixture::KernelType::ScalarType::Width::value];
    }
    oss << ']';

    return oss.str();
}


static Matrix2x2 *Matrix2x2_fromList(bp::list list) {
    if (bp::len(list) == 4) {
        Float buf[2][2];
        for (int i=0; i<2; ++i) {
            bp::list subList = bp::extract<bp::list>(list[i]);
            if (bp::len(subList) != 2)
                SLog(EError, "Matrix2x2 list constructor: invalid argument");
            for (int j=0; j<2; ++j)
                buf[i][j] = bp::extract<Float>(subList[j]);
        }
        return new Matrix2x2(buf);
    } else if (bp::len(list) == 4) {
        Float buf[4];
        for (int i=0; i<4; ++i)
            buf[i] = bp::extract<Float>(list[i]);
        return new Matrix2x2(buf);
    } else {
        SLog(EError, "Matrix2x2 list constructor: invalid argument");
        return NULL;
    }
}

static void Matrix2x2_setItem(Matrix2x2 *matrix, bp::tuple tuple, Float value) {
    if (bp::len(tuple) != 2)
        SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
    int i = bp::extract<int>(tuple[0]);
    int j = bp::extract<int>(tuple[1]);

    if (i < 0 || j < 0 || i >= 4 || j >= 4)
        SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

    matrix->operator()(i, j) = value;
}

static Float Matrix2x2_getItem(Matrix2x2 *matrix, bp::tuple tuple) {
    if (bp::len(tuple) != 2)
        SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
    int i = bp::extract<int>(tuple[0]);
    int j = bp::extract<int>(tuple[1]);

    if (i < 0 || j < 0 || i >= 2 || j >= 2)
        SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

    return matrix->operator()(i, j);
}

static Matrix2x2 Matrix2x2_transpose(Matrix2x2 *matrix) {
    Matrix2x2 result;
    matrix->transpose(result);
    return result;
}

static Matrix2x2 Matrix2x2_invert(Matrix2x2 *matrix) {
    Matrix2x2 result;
    matrix->invert(result);
    return result;
}

static bp::tuple Matrix2x2_symEig(Matrix2x2 *matrix) {
    Matrix2x2 Q;
    Float d[2];
    matrix->symEig(Q, d);

    bp::list list;
    for (int i=0; i<2; ++i)
        list.append(d[i]);

    return bp::make_tuple(Q, list);
}

static bp::tuple Matrix2x2_lu(Matrix2x2 *matrix) {
    Matrix2x2 LU;
    int piv[2];
    int pivsign;
    matrix->lu(LU, piv, pivsign);

    bp::list list;
    for (int i=0; i<2; ++i)
        list.append(piv[i]);

    return bp::make_tuple(LU, list, pivsign);
}

static Vector2 Matrix2x2_cholSolve(Matrix2x2 *matrix, Vector B) {
    typedef Matrix<2, 1, Float> Matrix2x1;
    Vector2 X;
    Matrix2x1 & aliasedB MTS_MAY_ALIAS = reinterpret_cast<Matrix2x1 &>(B);
    Matrix2x1 & aliasedX MTS_MAY_ALIAS = reinterpret_cast<Matrix2x1 &>(X);
    matrix->cholSolve<1>(aliasedB, aliasedX);
    return X;
}

static Vector2 Matrix2x2_luSolve(Matrix2x2 *matrix, Vector B, bp::list pivList) {
    typedef Matrix<2, 1, Float> Matrix2x1;
    Vector2 X;
    int piv[2];

    if (bp::len(pivList) != 2)
        SLog(EError, "Matrix2x2 list constructor: invalid argument");
    for (int i=0; i<2; ++i)
        piv[i] = bp::extract<Float>(pivList[i]);

    Matrix2x1 & aliasedB MTS_MAY_ALIAS = reinterpret_cast<Matrix2x1 &>(B);
    Matrix2x1 & aliasedX MTS_MAY_ALIAS = reinterpret_cast<Matrix2x1 &>(X);
    matrix->luSolve<1>(aliasedB, aliasedX, piv);
    return X;
}

/*
static void scalarSetItem(VMFMixture::KernelType::ScalarType& scalar, uint32_t index, float value) {
    if (index >= VMFMixture::KernelType::ScalarType::Width::value)
        SLog(EError, "Index %zu is out of bounds!", index);

    scalar.insert(index, value);
}

static float scalarGetItem(const VMFMixture::KernelType::ScalarType& scalar, uint32_t index) {
    if (index >= VMFMixture::KernelType::ScalarType::Width::value)
        SLog(EError, "Index %zu is out of bounds!", index);

    return scalar[index];
}
*/

static std::vector<PathGuidingSampleData>& appendToSampleData(std::vector<PathGuidingSampleData>& samples, const std::vector<PathGuidingSampleData>& range)
{
    samples.insert(samples.end(), range.begin(), range.end());
    return samples;
}

void export_guiding(){
	bp::object guidingModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mitsuba.guiding"))));
	bp::scope().attr("guiding") = guidingModule;
	PyObject *oldScope = bp::detail::current_scope;

	BP_SETSCOPE(guidingModule);
	guidingModule.attr("__path__") = "mitsuba.guiding";

    bp::class_<std::vector<float>>("VectorOfFloat")
            .def(bp::vector_indexing_suite<std::vector<float>>());

    bp::class_<std::vector<PathGuidingSampleData>>("VectorOfPathGuidingSampleData")
            .def(bp::vector_indexing_suite<std::vector<PathGuidingSampleData>>())
            .def("__iadd__", &appendToSampleData, BP_RETURN_VALUE);

    bp::class_<PathGuidingSampleData>("PathGuidingSampleData")
            .def(bp::init<>())
            .def(bp::init<Point, Vector, float, float, float>())
            .def_readwrite("position",  &PathGuidingSampleData::position)
            .def_readwrite("direction", &PathGuidingSampleData::direction)
            .def_readwrite("weight",    &PathGuidingSampleData::weight)
            .def_readwrite("pdf",       &PathGuidingSampleData::pdf)
            .def_readwrite("distance",  &PathGuidingSampleData::distance)
            .def("__repr__",  &PathGuidingSampleData::toString);

    bp::class_<lightpmm::VMMFactoryProperties>("VMMFactoryProperties")
            .def(bp::init<>())
            .def_readwrite("minItr", &lightpmm::VMMFactoryProperties::minItr)
            .def_readwrite("maxItr", &lightpmm::VMMFactoryProperties::maxItr)
            .def_readwrite("relLogLikelihoodThreshold", &lightpmm::VMMFactoryProperties::relLogLikelihoodThreshold)
            .def_readwrite("numInitialComponents", &lightpmm::VMMFactoryProperties::numInitialComponents)
            .def_readwrite("initKappa", &lightpmm::VMMFactoryProperties::initKappa)
            .def_readwrite("maxKappa", &lightpmm::VMMFactoryProperties::maxKappa)
            .def_readwrite("vPrior", &lightpmm::VMMFactoryProperties::vPrior)
            .def_readwrite("rPrior", &lightpmm::VMMFactoryProperties::rPrior)
            .def_readwrite("rPriorWeight", &lightpmm::VMMFactoryProperties::rPriorWeight)
            .def("__repr__", &lightpmm::VMMFactoryProperties::toString);

    using VMMFactory = lightpmm::VMMFactory<VMFMixture>;

    bp::class_<VMMFactory>("VMMFactory")
            .def(bp::init<>())
            .def(bp::init<const lightpmm::VMMFactoryProperties&>())
            .def_readwrite("properties", &VMMFactory::m_properties)
            .def("initialize",  &VMMFactory::initialize)
            .def("fit",         &VMMFactory::template fit<std::vector<PathGuidingSampleData>::const_iterator>)
            .def("maskedFit",   &VMMFactory::template maskedFit<std::vector<PathGuidingSampleData>::const_iterator>)
            .def("updateFit",   &VMMFactory::template updateFit<std::vector<PathGuidingSampleData>::const_iterator>)
            .def("serialize",   &VMMFactory::serialize)
            .def("deserialize", &VMMFactory::deserialize)
            .def("__repr__",    &VMMFactory::toString);

    bp::class_<VMFMixture>("VMFMixture")
            .def(bp::init<>())
            .def(bp::init<const VMFMixture&>())
            .def_readonly("MaxK", VMFMixture::MaxK::value)
            .add_property("K",   &VMFMixture::getK, &VMFMixture::setK)
            .def_readwrite("sampleWeight",    &VMFMixture::m_sampleWeight)
            .def_readwrite("numSamples",      &VMFMixture::m_numSamples)
            .def_readwrite("totalNumSamples", &VMFMixture::m_totalNumSamples)
            .def_readwrite("numEMIterations", &VMFMixture::m_numEMIterations)
            .def("sample",          &VMFMixture::sample, BP_RETURN_VALUE)
            .def("pdf",             &VMFMixture::pdf, BP_RETURN_VALUE)
            .def("getKappa",        &vmfMixtureGetKappa, BP_RETURN_VALUE)
            .def("setKappa",        &vmfMixtureSetKappa)
            .def("getMu",           &vmfMixtureGetMu, BP_RETURN_VALUE)
            .def("setMu",           &vmfMixtureSetMu)
            .def("getWeight",       &vmfMixtureGetWeight, BP_RETURN_VALUE)
            .def("setWeight",       &vmfMixtureSetWeight)
            .def("pdfComponent",    &vmfMixturePDFComponent, BP_RETURN_VALUE)
            .def("normWeights",     &VMFMixture::normWeights)
            .def("mergeComponents", &VMFMixture::mergeComponents)
            .def("removeWeightPrior", &VMFMixture::removeWeightPrior)
            .def("applyWeightPrior", &VMFMixture::applyWeightPrior)
            .def("product",         &VMFMixture::product,BP_RETURN_VALUE)
            .def("convolve",        &VMFMixture::convolve)
            .add_property("valid",  &VMFMixture::valid)
            .def("__repr__",        &VMFMixture::toString);

    bp::class_<VMFMixture::KernelType>("VMFKernel")
            .def(bp::init<>())
            .def(bp::init<float, float, Vector3>())
            //TODO: pdf needs Scalar type exposed to Python
            .def("pdf", &VMFMixture::KernelType::pdf)
            .def("sample", &VMFMixture::KernelType::sample)
            .add_property("valid", &VMFMixture::KernelType::valid)
            .def("__repr__", &VMFMixture::KernelType::toString);

    using SufficientStats = VMFMixture::KernelType::SufficientStats;

    bp::class_<SufficientStats>("SufficientStats")
            .def(bp::init<>())
            .def_readwrite("muTimesAvgCosineTimesSumWeight", &SufficientStats::muTimesAvgCosineTimesSumWeight)
            .def_readwrite("sumWeight", &SufficientStats::sumWeight)
            .def("accumulateSample", &SufficientStats::accumulateSample<PathGuidingSampleData>)
            .def("__repr__", &SufficientStats::toString);

    //eposing SSE vectors to python is buggy.
    /*
    using Scalar = VMFMixture::KernelType::ScalarType;

    bp::class_<Scalar>("Scalar")
            .def(bp::init<>())
            .def(bp::init<float>())
            .def("__setitem__", &scalarSetItem)
            .def("__getitem__", &scalarGetItem)
            .def("__repr__", &Scalar::toString);

    using VecOfScalar = typename lightpmm::Vec<Scalar>;

    bp::class_<VecOfScalar>("VecOfScalar")
            .def(bp::init<>())
            .def(bp::init<Vector>())
            .def_readwrite("x", &VecOfScalar::x)
            .def_readwrite("y", &VecOfScalar::y)
            .def_readwrite("z", &VecOfScalar::z)
            .add_property("length", &VecOfScalar::length)
            .def(bp::self + bp::self)
            .def(bp::self += bp::self)
            .def(bp::self - bp::self)
            .def(bp::self -= bp::self)
            .def(bp::self * bp::self)
            .def(bp::self *= bp::self)
            .def(bp::self * Scalar())
            .def(bp::self *= Scalar())
            .def(bp::self / bp::self)
            .def(bp::self /= bp::self)
            .def(bp::self / Scalar())
            .def(bp::self /= Scalar())
            .def("__repr__", &VecOfScalar::toString);
    */

    bp::implicitly_convertible<float, lightpmm::Scalar4>();
    bp::implicitly_convertible<Vector, lightpmm::Vec<lightpmm::Scalar4>>();
    bp::implicitly_convertible<float, lightpmm::Scalar8>();
    bp::implicitly_convertible<Vector, lightpmm::Vec<lightpmm::Scalar8>>();

    bp::class_<VMFMixture::SoftAssignmentWeights>("SoftAssignmentWeights", bp::no_init)
            .def(bp::init<VMFMixture, Vector>())
            .def("__getitem__",         &VMFMixture::SoftAssignmentWeights::getSoftAssignment)
            .def("getSoftAssignments",  &VMFMixture::SoftAssignmentWeights::getSoftAssignments)
            .add_property("valid",      &VMFMixture::SoftAssignmentWeights::valid)
            .add_property("mixturePDF", &VMFMixture::SoftAssignmentWeights::getMixturePDF)
            .def("__repr__",            &VMFMixture::SoftAssignmentWeights::toString);

    bp::class_<GuidingFieldType>("GuidingField")
            .def(bp::init<const Properties&>())
            .def(bp::init<Stream*>())
            .def("serialize", &GuidingFieldType::serialize)
            //TODO: getVolumeGuidingDistribution, configureGuidedBSDF
            .def("initField", &GuidingFieldType::initField)
            .def("buildField", static_cast<void (GuidingFieldType::*)(const AABB&, std::vector<PathGuidingSampleData>&)>(&GuidingFieldType::buildField))
            .def("updateField", static_cast<void (GuidingFieldType::*)(std::vector<PathGuidingSampleData>&, std::vector<Point>&, const RenderJob*)>(&GuidingFieldType::updateField))
            .def("addTrainingIteration", &GuidingFieldType::addTrainingIteration)
            .def_readwrite("guidingTree", &GuidingFieldType::m_guidingTree)
            .def_readwrite("factory", &GuidingFieldType::m_factory)
            .def_readwrite("iteration", &GuidingFieldType::m_iteration)
            .def_readwrite("totalSPP", &GuidingFieldType::m_totalSPP)
            .def("__repr__", &GuidingFieldType::toString);

    bp::class_<GuidingFieldType::GuidingTreeType>("GuidingTree")
            .def("getNumRegions", &GuidingFieldType::GuidingTreeType::getNumRegions)
            .def("getRegionIndex", &GuidingFieldType::GuidingTreeType::getRegionIndex)
            .def("getRegion", static_cast<GuidingFieldType::GuidingRegionType&(GuidingFieldType::GuidingTreeType::*)(Point)>(&GuidingFieldType::GuidingTreeType::getRegion), BP_RETURN_INTREF)
            .def("getRegion", static_cast<GuidingFieldType::GuidingRegionType&(GuidingFieldType::GuidingTreeType::*)(size_t)>(&GuidingFieldType::GuidingTreeType::getRegion), BP_RETURN_INTREF)
            .def("__repr__", &GuidingFieldType::GuidingTreeType::toString);

    bp::class_<GuidingFieldType::GuidingRegionType>("GuidingRegion")
            .def(bp::init<>())
            .def(bp::init<const GuidingFieldType::GuidingRegionType&>())
            .def_readwrite("distribution", &GuidingFieldType::GuidingRegionType::distribution)
            .def_readwrite("bound", &GuidingFieldType::GuidingRegionType::bound)
            .def_readwrite("statistics", &GuidingFieldType::GuidingRegionType::statistics)
            .def_readwrite("statsSinceLastSplit", &GuidingFieldType::GuidingRegionType::statsSinceLastSplit)
            .def_readwrite("statsUntilLastSplit", &GuidingFieldType::GuidingRegionType::statsUntilLastSplit)
            .def_readwrite("dStart", &GuidingFieldType::GuidingRegionType::dStart)
            .def_readwrite("nSamples", &GuidingFieldType::GuidingRegionType::nSamples)
            .def_readwrite("valid", &GuidingFieldType::GuidingRegionType::valid)
            .def("__repr__", &GuidingFieldType::GuidingRegionType::toString);

    using PathGuidingMixtureStatsType = PathGuidingMixtureStats<VMFMixture>;

    bp::class_<PathGuidingMixtureStatsType>("PathGuidingMixtureStats")
            .def_readwrite("flux", &PathGuidingMixtureStatsType::flux)
            .def_readwrite("lastFitStats", &PathGuidingMixtureStatsType::lastFitStats)
            .def_readwrite("incrementalDistance", &PathGuidingMixtureStatsType::incrementalDistance)
            .def_readwrite("incrementalPearsonChiSquared", &PathGuidingMixtureStatsType::incrementalPearsonChiSquared)
            .def_readwrite("incrementalCovariance2D", &PathGuidingMixtureStatsType::incrementalCovariance2D)
            .def(bp::self *= float())
            .def("clear", &PathGuidingMixtureStatsType::clear)
            .def("__repr__", &PathGuidingMixtureStatsType::toString);

    bp::class_<std::array<typename VMFMixture::KernelType::ScalarType::BooleanType, VMFMixture::NumKernels::value>>("VMFMixtureComponentMask");

    using GuidingFieldFactoryType = GuidingFieldType::GuidingFieldFactoryType;

    bp::class_<GuidingFieldFactoryType>("GuidingFieldFactory")
            .def(bp::init<>())
            .def(bp::init<const Properties&>())
            .def_readwrite("distributionFactory", &GuidingFieldFactoryType::m_distributionFactory)
            .def_readwrite("statisticsFactory", &GuidingFieldFactoryType::m_statisticsFactory)
            .def("fit", &GuidingFieldFactoryType::fit<std::vector<PathGuidingSampleData>>)
            .def("updateFit", &GuidingFieldFactoryType::updateFit<std::vector<PathGuidingSampleData>>)
            .def("__repr__", &GuidingFieldFactoryType::toString);

    using GuidingDistributionFactoryType = typename GuidingFieldFactoryType::DistributionFactoryType;

    bp::class_<GuidingDistributionFactoryType>("GuidingDistributionFactory")
            .def(bp::init<>())
            .def(bp::init<const Properties&>())
            .def_readwrite("vmmFactory", &GuidingDistributionFactoryType::m_vmmFactory)
            .def("serialize", &GuidingDistributionFactoryType::serialize)
            .def("initialize", &GuidingDistributionFactoryType::initialize)
            .def("fit", &GuidingDistributionFactoryType::fit<std::vector<PathGuidingSampleData>>)
            .def("updateFit", &GuidingDistributionFactoryType::updateFit<std::vector<PathGuidingSampleData>>)
            .def("__repr__", &GuidingDistributionFactoryType::toString);

    using GuidingStatisticsFactoryType = typename GuidingFieldFactoryType::StatisticsFactoryType;

    bp::class_<GuidingStatisticsFactoryType>("GuidingStatisticsFactory")
            .def(bp::init<>())
            .def(bp::init<const Properties&>())
            .def_readwrite("useSplitAndMerge", &GuidingStatisticsFactoryType::m_splitAndMerge)
            .def_readwrite("splitMinDivergence", &GuidingStatisticsFactoryType::m_splitMinDivergence)
            .def_readwrite("mergeMaxDivergence", &GuidingStatisticsFactoryType::m_mergeMaxDivergence)
            .def_readwrite("minSamplesForSplitting", &GuidingStatisticsFactoryType::m_minSamplesForSplitting)
            .def_readwrite("minSamplesForPostSplitFitting", &GuidingStatisticsFactoryType::m_minSamplesForPostSplitFitting)
            .def("serialize", &GuidingStatisticsFactoryType::serialize)
            .def("preFit", &GuidingStatisticsFactoryType::preFit<std::vector<PathGuidingSampleData>>)
            .def("postFit", &GuidingStatisticsFactoryType::postFit<std::vector<PathGuidingSampleData>>)
            .def("preSample", &GuidingStatisticsFactoryType::preSample)
            .def("mergeComponents", &GuidingStatisticsFactoryType::mergeComponents).staticmethod("mergeComponents")
            .def("mergeOne", &GuidingStatisticsFactoryType::mergeOne)
            .def("mergeAll", &GuidingStatisticsFactoryType::mergeAll)
            .def("splitComponentUsingPCA", &GuidingStatisticsFactoryType::splitComponentUsingPCA).staticmethod("splitComponentUsingPCA")
            .def("splitOne", &GuidingStatisticsFactoryType::splitOne<std::vector<PathGuidingSampleData>>)
            .def("splitAll", &GuidingStatisticsFactoryType::splitAll<std::vector<PathGuidingSampleData>>)
            .def("repositionSample", &GuidingStatisticsFactoryType::repositionSample<PathGuidingSampleData>).staticmethod("repositionSample")
            .def("computePearsonChiSquaredMergeMetric", &pathGuidingMixtureStatsComputePearsonChiSquaredMergeMetric)
            .def("computeMaskFromComponentIndexVector", &computeMaskFromComponentIndexVector).staticmethod("computeMaskFromComponentIndexVector")
            .def("__repr__", &GuidingStatisticsFactoryType::toString);

    bp::class_<Guiding::SampleStats>("SampleStats")
            .def_readwrite("mean", &Guiding::SampleStats::mean)
            .def_readwrite("unnormalizedVariance", &Guiding::SampleStats::unnormalizedVariance)
            .def_readonly("variance", &Guiding::SampleStats::computeVariance)
            .def_readwrite("numSamples", &Guiding::SampleStats::numSamples)
            .def_readwrite("numZeroValuedSamples", &Guiding::SampleStats::numZeroValuedSamples)
            .def_readwrite("sampleBound", &Guiding::SampleStats::sampleBound)
            .def("__repr__", &Guiding::SampleStats::toString);

    bp::class_<PPG::STree>("PPG_STree")
            .def("findDtree", &PPG::STree::findDTree, BP_RETURN_VALUE)
            .def("__repr__", &PPG::STree::toString);

    bp::class_<PPG::DTree>("PPG_DTree")
            .def("pdf", &PPG::DTree::pdf)
            .add_property("aabb", bp::make_function(&PPG::DTree::aabb, BP_RETURN_VALUE));

    bp::class_<PathGuidingUtilities>("PathGuidingUtilities")
            .def_readwrite("renderSize", &PathGuidingUtilities::renderSize)
            .def("loadEnvmapForSampling", &PathGuidingUtilities::loadEnvmapForSampling).staticmethod("loadEnvmapForSampling")
            .def("sampleEnvMap", &PathGuidingUtilities::sampleEnvMap).staticmethod("sampleEnvMap")
            .def("renderEnvmapPDF", &PathGuidingUtilities::renderEnvmapPDF).staticmethod("renderEnvmapPDF")
            .def("getFirstSmoothSurfaceInteraction", &PathGuidingUtilities::getFirstSmoothSurfaceInteraction, BP_RETURN_VALUE).staticmethod("getFirstSmoothSurfaceInteraction")
            .def("renderSphericalView", &PathGuidingUtilities::renderSphericalView).staticmethod("renderSphericalView")
            .def("renderSphericalViewDistance", &PathGuidingUtilities::renderSphericalViewDistance).staticmethod("renderSphericalViewDistance")
            .def("renderBSDF", &PathGuidingUtilities::renderBSDF).staticmethod("renderBSDF")
            .def("normalizeSphericalView", &PathGuidingUtilities::normalizeSphericalView).staticmethod("normalizeSphericalView")
            .def("exportSamplesAsCSV", &PathGuidingUtilities::exportSamplesAsCSV).staticmethod("exportSamplesAsCSV")
            .def("exportSamplesAsOBJ", &PathGuidingUtilities::exportSamplesAsOBJ).staticmethod("exportSamplesAsOBJ")
            .def("renderSampleBitmap", &PathGuidingUtilities::renderSampleBitmap).staticmethod("renderSampleBitmap")
            .def("renderSamplePDFBitmap", &PathGuidingUtilities::renderSamplePDFBitmap).staticmethod("renderSamplePDFBitmap")
            .def("renderSampleDistanceBitmap", &PathGuidingUtilities::renderSampleDistanceBitmap).staticmethod("renderSampleDistanceBitmap")
            .def("loadSamples", &PathGuidingUtilities::loadSamples).staticmethod("loadSamples")
            .def("loadSampleRange", &PathGuidingUtilities::loadSampleRange).staticmethod("loadSampleRange")
            .def("reweightSamplesWithComponentSoftAssignment", &PathGuidingUtilities::reweightSamplesWithComponentSoftAssignment).staticmethod("reweightSamplesWithComponentSoftAssignment")
            .def("reweightSamplesWithRelativePDF", &PathGuidingUtilities::reweightSamplesWithRelativePDF).staticmethod("reweightSamplesWithRelativePDF")
            .def("exportVMMAsCSV", &PathGuidingUtilities::exportVMMAsCSV).staticmethod("exportVMMAsCSV")
            .def("renderVMMPDF", &PathGuidingUtilities::renderVMMPDF).staticmethod("renderVMMPDF")
            .def("renderVMMDistance", &PathGuidingUtilities::renderVMMDistance).staticmethod("renderVMMDistance")
            .def("renderVMFPDF", &PathGuidingUtilities::renderVMFPDF).staticmethod("renderVMFPDF")
            .def("renderVMMPDFColored", &PathGuidingUtilitiesRenderVMMPDFColored).staticmethod("renderVMMPDFColored")
            .def("computeLogLikelihood", &PathGuidingUtilities::computeLogLikelihood).staticmethod("computeLogLikelihood")
            .def("computePearsonChiSquaredDivergence", &PathGuidingUtilities::computePearsonChiSquaredDivergence).staticmethod("computePearsonChiSquaredDivergence")
            .def("exportGuidingTreeASCII", &PathGuidingUtilities::exportGuidingTreeASCII).staticmethod("exportGuidingTreeASCII")
            .def("exportGuidingTreeOBJ", &PathGuidingUtilities::exportGuidingTreeOBJ).staticmethod("exportGuidingTreeOBJ")
            .def("exportGuidingTreeSampleBoundsOBJ", &PathGuidingUtilities::exportGuidingTreeSampleBoundsOBJ).staticmethod("exportGuidingTreeSampleBoundsOBJ")
            .def("exportGuidingTreeCentersOBJ", &PathGuidingUtilities::exportGuidingTreeCentersOBJ).staticmethod("exportGuidingTreeCentersOBJ")
            .def("loadPPGSDTree", &PathGuidingUtilities::loadPPGSDTree).staticmethod("loadPPGSDTree")
            .def("renderDTreePDF", &PathGuidingUtilities::renderDTreePDF).staticmethod("renderDTreePDF")
            .def("exportSTreeOBJ", &PathGuidingUtilities::exportSTreeOBJ).staticmethod("exportSTreeOBJ");

    bp::class_<IncrementalPearsonChiSquared<VMFMixture>>("IncrementalPearsonChiSquared")
            .def_readwrite("divergencePlusOneTimesIntegralSqr", &IncrementalPearsonChiSquared<VMFMixture>::divergencePlusOneTimesIntegralSqr)
            .def_readwrite("numSamples", &IncrementalPearsonChiSquared<VMFMixture>::numSamples)
            .def("updateDivergence", &IncrementalPearsonChiSquared<VMFMixture>::updateDivergence<std::vector<PathGuidingSampleData>>)
            .def("updateDivergenceMasked", &IncrementalPearsonChiSquared<VMFMixture>::updateDivergenceMasked<std::vector<PathGuidingSampleData>>)
            .def("computeDivergence", &incrementalPearsonChiSquaredComputeDivergence)
            .def("merge", &IncrementalPearsonChiSquared<VMFMixture>::merge)
            .def("split", &IncrementalPearsonChiSquared<VMFMixture>::split)
            .def(bp::self *= float())
            .def("clear", &IncrementalPearsonChiSquared<VMFMixture>::clear)
            .def("__repr__", &IncrementalPearsonChiSquared<VMFMixture>::toString);

    bp::class_<IncrementalCovariance2D<VMFMixture>>("IncrementalCovariance2D")
            .def_readwrite("varianceAndCovariance", &IncrementalCovariance2D<VMFMixture>::varianceAndCovariance)
            .def_readwrite("sumWeights", &IncrementalCovariance2D<VMFMixture>::sumWeights)
            .def("updateStatistics", &IncrementalCovariance2D<VMFMixture>::updateStatistics<std::vector<PathGuidingSampleData>>)
            .def("computeCovarianceMatrix", &IncrementalCovariance2D<VMFMixture>::computeCovarianceMatrix)
            .def("merge", &IncrementalCovariance2D<VMFMixture>::merge)
            .def("split", &IncrementalCovariance2D<VMFMixture>::split)
            .def(bp::self *= float())
            .def("clear", &IncrementalCovariance2D<VMFMixture>::clear)
            .def("__repr__", &IncrementalCovariance2D<VMFMixture>::toString);

    bp::class_<IncrementalDistance<VMFMixture>>("IncrementalDistance")
            .def_readwrite("distances", &IncrementalDistance<VMFMixture>::distances)
            .def_readwrite("sumWeights", &IncrementalDistance<VMFMixture>::sumWeights)
            .def("updateDistances", &IncrementalDistance<VMFMixture>::updateDistances<std::vector<PathGuidingSampleData>>)
            .def("reposition", &IncrementalDistance<VMFMixture>::reposition)
            .def("repositionDistribution", &IncrementalDistance<VMFMixture>::repositionDistribution)
            .def("merge", &IncrementalDistance<VMFMixture>::merge)
            .def("split", &IncrementalDistance<VMFMixture>::split)
            .def(bp::self *= float())
            .def("clear", &IncrementalDistance<VMFMixture>::clear)
            .def("__repr__", &IncrementalDistance<VMFMixture>::toString);

    bp::class_<std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>>("ScalarArray")
            .def(bp::init<>())
            .def("__init__", bp::make_constructor(ScalarArray_fromList))
            .def("__getitem__", &ScalarArray_getItem)
            .def("__setitem__", &ScalarArray_setItem)
            .def("__repr__", &ScalarArray_repr);

    //bp::implicitly_convertible<bp::list, std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>>();

    bp::class_<Matrix2x2>("Matrix2x2", bp::init<>())
            .def(bp::init<Float>())
            .def(bp::init<Stream *>())
            .def(bp::init<Matrix2x2>())
            .def(bp::init<Vector2, Vector2>())
            .def("__init__", bp::make_constructor(Matrix2x2_fromList))
            .def(bp::init<Stream *>())
            .def("__setitem__", &Matrix2x2_setItem)
            .def("__getitem__", &Matrix2x2_getItem)
            .def("setIdentity", &Matrix2x2::setIdentity)
            .def("setZero", &Matrix2x2::setZero)
            .def("isZero", &Matrix2x2::isZero)
            .def("isIdentity", &Matrix2x2::isIdentity)
            .def("trace", &Matrix2x2::trace)
            .def("frob", &Matrix2x2::frob)
            .def("det", &Matrix2x2::det)
            .def("cholDet", &Matrix2x2::cholDet)
            .def("luDet", &Matrix2x2::luDet)
            .def("chol", &Matrix2x2::chol)
            .def("symEig", &Matrix2x2_symEig)
            .def("lu", &Matrix2x2_lu)
            .def("cholSolve", &Matrix2x2_cholSolve)
            .def("luSolve", &Matrix2x2_luSolve)
            .def("transpose", &Matrix2x2_transpose)
            .def("invert", &Matrix2x2_invert)
            .def("serialize", &Matrix2x2::serialize)
            .def(bp::self != bp::self)
            .def(bp::self == bp::self)
            .def(-bp::self)
            .def(bp::self + bp::self)
            .def(bp::self += bp::self)
            .def(bp::self + Float())
            .def(bp::self += Float())
            .def(bp::self - bp::self)
            .def(bp::self -= bp::self)
            .def(bp::self - Float())
            .def(bp::self -= Float())
            .def(bp::self * Float())
            .def(Float() * bp::self)
            .def(bp::self *= Float())
            .def(bp::self * bp::self)
            .def(bp::self *= bp::self)
            .def(bp::self / Float())
            .def(bp::self /= Float())
            .def("__repr__", &Matrix2x2::toString);

	bp::detail::current_scope = oldScope;
}


