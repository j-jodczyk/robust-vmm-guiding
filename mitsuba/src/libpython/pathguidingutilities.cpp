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

#include "pathguidingutilities.h"

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/imageblock.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/renderqueue.h>

#include <pmm/DirectionalData.h>
#include <mitsuba/guiding/BSPTree.h>

#include <stack>
#include <vector>
#include <fstream>
#include <limits>

GUIDING_NAMESPACE_BEGIN

Vector2i PathGuidingUtilities::renderSize(512, 256);

ref<Emitter> PathGuidingUtilities::loadEnvmapForSampling(const std::string &filename)
{
    const Transform envmapTransform = Transform{
            Matrix4x4{0.0, 0.0, 1.0, 0.0,
                      1.0, 0.0, 0.0, 0.0,
                      0.0, 1.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 1.0}};

    Properties envmapProperties{"envmap"};
    envmapProperties.setString("filename", filename);
    envmapProperties.setTransform("toWorld", envmapTransform);

    ref<Emitter> envmap = static_cast<Emitter*> (PluginManager::getInstance()->createObject(MTS_CLASS(Emitter), envmapProperties));
    envmap->configure();

    //generate some geometry to set the envmap's bounding sphere to something non-zero
    //(this is needed for EnvironmentMap::sampleDirect)
    Properties sphereProperties{"sphere"};
    sphereProperties.setPoint("center", Point{0.0f});
    sphereProperties.setFloat("radius", 1.0f);
    ref<Shape> sphere = static_cast<Shape*> (PluginManager::getInstance()->createObject(MTS_CLASS(Shape), sphereProperties));
    ref<Scene> scene = new Scene();
    scene->addChild(sphere);
    scene->initialize();
    scene->configure();

    envmap->createShape(scene);

    return envmap;
}

std::vector<PathGuidingSampleData> PathGuidingUtilities::sampleEnvMap(const Emitter *envMap, unsigned int numSamples, bool useImportanceSampling, const Point p, const Point2i seed)
{
    std::vector<PathGuidingSampleData> samples;

    if (!envMap)
        return samples;

    ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("deterministic")));
    sampler->generate(seed);

    samples.reserve(numSamples);

    DirectSamplingRecord dRec(p, 0.0f);

    for (unsigned int i=0; i<numSamples; i++)
    {
        if (useImportanceSampling)
        {
            // importance sample the env map
            Spectrum liDivPDF = envMap->sampleDirect(dRec, sampler->next2D());

            samples.emplace_back(p, dRec.d, PATHGUIDING_SPECTRUM_TO_FLOAT(liDivPDF), dRec.pdf, std::numeric_limits<float>::infinity());
        }
        else
        {
            // sample the env map randomly
            Intersection its;
            its.p = p;
            const Vector d = warp::squareToUniformSphere(sampler->next2D());
            Spectrum liDivPDF = envMap->eval(its, d)/warp::squareToUniformSpherePdf();

            samples.emplace_back(its.p, -d, PATHGUIDING_SPECTRUM_TO_FLOAT(liDivPDF), warp::squareToUniformSpherePdf(), std::numeric_limits<float>::infinity());
        }

        sampler->advance();
    }

    return samples;
}

ref<Bitmap> PathGuidingUtilities::renderEnvmapPDF(const Emitter *envMap)
{
    ref<Bitmap> pdfBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    pdfBitmap->clear();

    if (!envMap)
        return pdfBitmap;

    DirectSamplingRecord dRec;
    dRec.measure = ESolidAngle;

    Point2i p;
    for (p.y=0; p.y<renderSize.y; p.y++)
    {
        for (p.x=0; p.x<renderSize.x; p.x++)
        {
            dRec.d = sphericalPlanePointToDirVector(p);
            pdfBitmap->setPixel(p, Spectrum(envMap->pdfDirect(dRec)));
        }
    }

    return pdfBitmap;
}

Intersection PathGuidingUtilities::getFirstSmoothSurfaceInteraction(const Scene *scene, Point2 pixel)
{
    SAssert(scene != nullptr);

    ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("independent")));
    sampler->generate(Point2i(4,2));

    Intersection its;
    Ray ray;

    scene->getSensor()->sampleRay(ray, pixel, sampler->next2D(), sampler->next1D());

    if (scene->rayIntersect(ray, its))
    {
        int depth = 1;
        const int maxDepth = 32;

        while (depth < maxDepth && (its.getBSDF()->getType() & BSDF::ESmooth) == 0)
        {
            BSDFSamplingRecord bRec(its, sampler, ERadiance);
            its.getBSDF()->sample(bRec, sampler->next2D());
            const Vector wo = its.toWorld(bRec.wo);
            ray = Ray(its.p, wo, ray.time);
            if (!scene->rayIntersect(ray, its))
                break;
            depth++;
        }
    }

    if (!its.isValid())
        SLog(EDebug, "no surface interaction at pixel location %s", pixel.toString().c_str());

    return its;
}

ref<Bitmap> PathGuidingUtilities::renderSphericalView(const Scene *scene, const Point p, const std::string& integratorName, int numSamples)
{
    ref<PluginManager> pluginManager = PluginManager::getInstance();

    Properties rFilterProps{"box"};
    ref<ReconstructionFilter> rFilter = static_cast<ReconstructionFilter*>(pluginManager->createObject(MTS_CLASS(ReconstructionFilter), rFilterProps));
    rFilter->configure();

    Properties sensorProps{"spherical"};
    sensorProps.setTransform("toWorld", Transform::lookAt(p, Point{p+Vector{1,0,0}}, Vector{0,0,1}));
    ref<Sensor> sensor = static_cast<Sensor*>(pluginManager->createObject(MTS_CLASS(Sensor), sensorProps));

    Properties filmProps{"hdrfilm"};
    filmProps.setInteger("width",  renderSize.x);
    filmProps.setInteger("height", renderSize.y);
    filmProps.setBoolean("banner", false);
    ref<Film> film = static_cast<Film*>(pluginManager->createObject(MTS_CLASS(Film), filmProps));
    film->clear();
    film->addChild(rFilter);
    film->setDestinationFile("", 0);
    film->configure();

    Properties samplerProps{"independent"};
    samplerProps.setInteger("sampleCount", numSamples);
    ref<Sampler> sampler = static_cast<Sampler*>(pluginManager->createObject(MTS_CLASS(Sampler), samplerProps));
    sampler->configure();

    sensor->addChild(sampler);
    sensor->addChild(film);
    sensor->configure();

    Properties integratorProps{integratorName};
    //disable features not compatible with a spherical camera
    integratorProps.setBoolean("lensPerturbation", false);
    integratorProps.setBoolean("multiChainPerturbation", false);
    integratorProps.setBoolean("causticPerturbation", false);
    integratorProps.setBoolean("bidirectionalMutation", true);
    integratorProps.setBoolean("manifoldPerturbation", true);
    if (integratorName == "pathguiding")
        integratorProps.setInteger("trainingSamples", std::max(0, numSamples-4));
    ref<Integrator> integrator = static_cast<Integrator*>(pluginManager->createObject(MTS_CLASS(Integrator), integratorProps));
    integrator->configure();

    ref<Scene> sceneCopy = new Scene{const_cast<Scene*>(scene)};
    sceneCopy->setSensor(sensor);
    sceneCopy->setIntegrator(integrator);
    sceneCopy->setSampler(sampler);
    sceneCopy->configure();

    ref<RenderQueue> renderQueue = new RenderQueue{};

    ref<RenderJob> job = new RenderJob{"sphericalView", sceneCopy, renderQueue, -1, -1, -1, false};
    job->start();

    renderQueue->waitLeft(0);
    renderQueue->join();

    ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    bitmap->clear();

    film->develop(Point2i{0,0}, renderSize, Point2i{0,0}, bitmap);

    return bitmap;
}

ref<Bitmap> PathGuidingUtilities::renderSphericalViewDistance(const Scene *scene, const Point p, int numSamples)
{
    ref<PluginManager> pluginManager = PluginManager::getInstance();

    Properties rFilterProps{"box"};
    ref<ReconstructionFilter> rFilter = static_cast<ReconstructionFilter*>(pluginManager->createObject(MTS_CLASS(ReconstructionFilter), rFilterProps));
    rFilter->configure();

    Properties sensorProps{"spherical"};
    sensorProps.setTransform("toWorld", Transform::lookAt(p, Point{p+Vector{1,0,0}}, Vector{0,0,1}));
    ref<Sensor> sensor = static_cast<Sensor*>(pluginManager->createObject(MTS_CLASS(Sensor), sensorProps));

    Properties filmProps{"hdrfilm"};
    filmProps.setInteger("width",  renderSize.x);
    filmProps.setInteger("height", renderSize.y);
    filmProps.setBoolean("banner", false);
    ref<Film> film = static_cast<Film*>(pluginManager->createObject(MTS_CLASS(Film), filmProps));
    film->clear();
    film->addChild(rFilter);
    film->setDestinationFile("", 0);
    film->configure();

    Properties samplerProps{"independent"};
    samplerProps.setInteger("sampleCount", numSamples);
    ref<Sampler> sampler = static_cast<Sampler*>(pluginManager->createObject(MTS_CLASS(Sampler), samplerProps));
    sampler->configure();

    sensor->addChild(sampler);
    sensor->addChild(film);
    sensor->configure();

    Properties integratorProps{"field"};
    integratorProps.setString("field", "distance");
    integratorProps.setFloat("undefined", 0.0f);

    ref<Integrator> integrator = static_cast<Integrator*>(pluginManager->createObject(MTS_CLASS(Integrator), integratorProps));
    integrator->configure();

    ref<Scene> sceneCopy = new Scene{const_cast<Scene*>(scene)};
    sceneCopy->setSensor(sensor);
    sceneCopy->setIntegrator(integrator);
    sceneCopy->setSampler(sampler);
    sceneCopy->configure();

    ref<RenderQueue> renderQueue = new RenderQueue{};

    ref<RenderJob> job = new RenderJob{"sphericalView", sceneCopy, renderQueue, -1, -1, -1, false};
    job->start();

    renderQueue->waitLeft(0);
    renderQueue->join();

    ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    bitmap->clear();

    film->develop(Point2i{0,0}, renderSize, Point2i{0,0}, bitmap);

    return bitmap;
}

ref<Bitmap> PathGuidingUtilities::renderBSDF(const Intersection &its)
{
    ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    Point2i p;
    for (p.y=0; p.y<renderSize.y; p.y++)
    {
        for (p.x=0; p.x<renderSize.x; p.x++)
        {
            BSDFSamplingRecord bRec(its, its.toLocal(sphericalPlanePointToDirVector(p)));
            bitmap->setPixel(p, Spectrum(its.getBSDF()->eval(bRec)));
        }
    }

    return bitmap;
}

ref<Bitmap> PathGuidingUtilities::normalizeSphericalView(const ref<Bitmap> bitmap)
{
    ref<Bitmap> normalizedBitmap = bitmap->clone();
    const Vector2i size = normalizedBitmap->getSize();

    float sum{0.0f};

    Point2i p;
    for (p.y=0; p.y<size.y; p.y++)
    {
        float rowSum{0.0f};

        for (p.x=0; p.x<size.x; p.x++)
        {
            rowSum += PATHGUIDING_SPECTRUM_TO_FLOAT(normalizedBitmap->getPixel(p));
        }

        const float sinTheta = sin((static_cast<float>(p.y)+0.5f)*M_PI/static_cast<float>(size.y));
        sum += rowSum*sinTheta;
    }

    normalizedBitmap->scale(static_cast<float>(size.x*size.y)/(4.0f*M_PI*sum));

    return normalizedBitmap;
}

void PathGuidingUtilities::exportSamplesAsCSV(const std::string &filename, const std::vector<PathGuidingSampleData>& samples)
{
    std::ofstream os(filename.c_str());

    for (const PathGuidingSampleData& sample : samples)
    {
        os << sample.position.x  << ';' << sample.position.y  << ';' << sample.position.z << ';'
           << sample.direction.x << ';' << sample.direction.y << ';' << sample.direction.z << ';'
           << sample.weight << ';'
           << sample.pdf << ';'
           << sample.distance
           << '\n';
    }

    os.close();

    SLog(EInfo, "%zu samples written to %s", samples.size(), filename.c_str());
}

void PathGuidingUtilities::exportSamplesAsOBJ(const std::string &filename, const std::vector<PathGuidingSampleData>& samples)
{
    std::ofstream os(filename.c_str());
    os << "o Samples\n";

    for (const PathGuidingSampleData& sample : samples)
    {
        const Point& p = sample.position;
        os << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
    }

    os.close();

    SLog(EInfo, "%zu samples written to %s", samples.size(), filename.c_str());
}

ref<Bitmap> PathGuidingUtilities::renderSampleBitmap(const std::vector<PathGuidingSampleData>& samples)
{
    ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    bitmap->clear();

    //accumulate samples
    for (const PathGuidingSampleData& sample : samples)
    {
        const Vector dir = sample.direction;
        const Point2i p = dirVectorToSphericalPlanePoint(dir);
        const Spectrum radiance(sample.weight*INV_TWOPI);
        bitmap->setPixel(p, bitmap->getPixel(p)+radiance);
    }

    //normalize to one sample per pixel on one hemisphere
    float normalization = static_cast<float>(renderSize.x*renderSize.y)*0.5f/static_cast<float>(samples.size());

    bitmap->scale(normalization);

    return bitmap;
}

ref<Bitmap> PathGuidingUtilities::renderSamplePDFBitmap(const std::vector<PathGuidingSampleData>& samples)
{
    ref<Bitmap> bitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, renderSize);
    bitmap->clear();

    std::vector<size_t> sampleCountPerPixel(renderSize.x*renderSize.y, 0);

    for (const PathGuidingSampleData& sample : samples)
    {
        const Vector dir = sample.direction;
        const Point2i p = dirVectorToSphericalPlanePoint(dir);
        const size_t index = p.y*renderSize.x+p.x;

        const float previousSamples = sampleCountPerPixel[index];
        ++sampleCountPerPixel[index];

        const Spectrum pdf(sample.pdf);
        bitmap->setPixel(p, ((bitmap->getPixel(p)*previousSamples+pdf)/(previousSamples+1.0f)));
    }

    bitmap->scale(1.0f/(bitmap->average().average()*static_cast<float>(renderSize.x*renderSize.y)));

    return bitmap;
}

ref<Bitmap> PathGuidingUtilities::renderSampleDistanceBitmap(const std::vector<PathGuidingSampleData>& samples)
{
    ref<Bitmap> bitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, renderSize);
    bitmap->clear();

    std::vector<float> pixelWeight(renderSize.x*renderSize.y, 0.0f);

    //accumulate samples
    for (const PathGuidingSampleData& sample : samples)
    {
        const Vector dir = sample.direction;
        const Point2i p = dirVectorToSphericalPlanePoint(dir);
        const Spectrum depth(sample.distance);

        bitmap->setPixel(p, bitmap->getPixel(p)+sample.weight*depth);
        pixelWeight[p.x+p.y*renderSize.x] += sample.weight;
    }

    Point2i p;
    for (p.y=0; p.y<renderSize.y; p.y++)
        for (p.x=0; p.x<renderSize.x; p.x++)
            if (pixelWeight[p.x+p.y*renderSize.x] > 0.0f)
                bitmap->setPixel(p, bitmap->getPixel(p)/pixelWeight[p.x+p.y*renderSize.x]);

    return bitmap;
}

std::vector<PathGuidingSampleData> PathGuidingUtilities::loadSamples(const std::string &filename)
{
    static_assert(std::is_trivially_copyable<PathGuidingSampleData>::value, "cannot trivially deserialize sample data");

    std::vector<PathGuidingSampleData> sampleData;

    if (fs::exists(filename))
    {
        ref<FileStream> sampleFileStream = new FileStream(filename, FileStream::EReadOnly);

        size_t availableSamples = sampleFileStream->readSize();
        sampleData.resize(availableSamples);
        sampleFileStream->read(sampleData.data(), availableSamples*sizeof(PathGuidingSampleData));

        sampleFileStream->close();
    }
    else
        SLog(EWarn, "Serialized sample data could not be opened! (%s)", filename.c_str());

    return sampleData;
}

std::vector<PathGuidingSampleData> PathGuidingUtilities::loadSampleRange(const std::string &filename, size_t startIndex, size_t numSamples)
{
    static_assert(std::is_trivially_copyable<PathGuidingSampleData>::value, "cannot trivially deserialize sample data");

    std::vector<PathGuidingSampleData> sampleData;

    if (fs::exists(filename))
    {
        ref<FileStream> sampleFileStream = new FileStream(filename, FileStream::EReadOnly);

        size_t availableSamples = sampleFileStream->readSize();

        if (startIndex+numSamples >= startIndex && startIndex+numSamples <= availableSamples)
        {
            size_t startingPos = sampleFileStream->getPos();
            startingPos += startIndex*sizeof(PathGuidingSampleData);
            sampleFileStream->seek(startingPos);

            sampleData.resize(numSamples);
            sampleFileStream->read(sampleData.data(), numSamples*sizeof(PathGuidingSampleData));
        }
        else
        {
            SLog(EWarn, "requested sample range out of available data range (%zu-%zu requested, %zu available)", startIndex, startIndex+numSamples, availableSamples);
        }

        sampleFileStream->close();
    }
    else
        SLog(EWarn, "Serialized sample data could not be opened! (%s)", filename.c_str());

    return sampleData;
}

std::vector<PathGuidingSampleData> PathGuidingUtilities::reweightSamplesWithComponentSoftAssignment(const VMFMixture &vmm, size_t component, const std::vector<PathGuidingSampleData>& samples)
{
    std::vector<PathGuidingSampleData> reweightedSamples{samples.begin(), samples.end()};
    for (PathGuidingSampleData& sample : reweightedSamples)
    {
        const VMFMixture::SoftAssignmentWeights softAssignments {vmm, sample.direction};
        sample.weight *= softAssignments.getSoftAssignment(component);
    }

    return reweightedSamples;
}

std::vector<PathGuidingSampleData> PathGuidingUtilities::reweightSamplesWithRelativePDF(const VMFMixture &vmm, const std::vector<PathGuidingSampleData>& samples)
{
    const float avgWeight = std::accumulate(samples.begin(), samples.end(), 0.0f, [](float sum, PathGuidingSampleData sample){return sum+sample.weight;})/static_cast<float>(samples.size());

    std::vector<PathGuidingSampleData> reweightedSamples{samples.begin(), samples.end()};
    for (PathGuidingSampleData& sample : reweightedSamples)
    {
        sample.weight *= sample.pdf/(vmm.pdf(sample.direction)*avgWeight);
    }

    return reweightedSamples;
}

void PathGuidingUtilities::exportVMMAsCSV(const std::string &filename, const VMFMixture &vmm)
{
    std::ofstream os(filename.c_str());

    for (unsigned int i=0; i<vmm.getK(); i++)
    {
        const Vector& mu = vmm.m_comps[i/VMFMixture::KernelType::ScalarType::Width::value].getMu(i%VMFMixture::KernelType::ScalarType::Width::value);

        os << mu[0] << ';' << mu[1] << ';' << mu[2] << ';'
           << vmm.m_comps[i/VMFMixture::KernelType::ScalarType::Width::value].getKappa(i%VMFMixture::KernelType::ScalarType::Width::value) << ';'
           << vmm.m_comps[i/VMFMixture::KernelType::ScalarType::Width::value].getWeight(i%VMFMixture::KernelType::ScalarType::Width::value) << '\n';
    }

    os.close();

    SLog(EInfo, "vMF mixture written to %s", filename.c_str());
}

ref<Bitmap> PathGuidingUtilities::renderVMMPDF(const VMFMixture &vmm)
{
    ref<Bitmap> pdfBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    for (int y=0; y<renderSize.y; y++)
    {
        for (int x=0; x<renderSize.x; x++)
        {
            const Point2i p{x,y};
            const Vector dir = sphericalPlanePointToDirVector(p);

            pdfBitmap->setPixel(p, Spectrum(vmm.pdf(dir)));
        }
    }

    Spectrum vMFCenterHighlight(0.0f);
    vMFCenterHighlight[0] = 1.0f;
    Spectrum vMFCenterData(0.0f);

    for (unsigned int i=0; i<vmm.getK(); i++)
    {
        vMFCenterData[1] = vmm.m_comps[i/VMFMixture::KernelType::ScalarType::Width::value].getWeight(i%VMFMixture::KernelType::ScalarType::Width::value);
        vMFCenterData[2] = vmm.m_comps[i/VMFMixture::KernelType::ScalarType::Width::value].getKappa(i%VMFMixture::KernelType::ScalarType::Width::value);

        const Vector dir = vmm.m_comps[i/VMFMixture::KernelType::ScalarType::Width::value].getMu(i%VMFMixture::KernelType::ScalarType::Width::value);

        const Point2i p = dirVectorToSphericalPlanePoint(dir);

        pdfBitmap->setPixel(p, pdfBitmap->getPixel(p)*vMFCenterHighlight+vMFCenterData);
    }

    return pdfBitmap;
}

ref<Bitmap> PathGuidingUtilities::renderVMMDistance(const VMFMixture &vmm, const std::array<VMFMixture::KernelType::ScalarType, VMFMixture::NumKernels::value>& distances)
{
    ref<Bitmap> distanceBitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, renderSize);
    for (int y=0; y<renderSize.y; y++)
    {
        for (int x=0; x<renderSize.x; x++)
        {
            const Point2i p{x,y};
            const Vector dir = sphericalPlanePointToDirVector(p);

            VMFMixture::SoftAssignmentWeights softAssignments {vmm, dir};

            float pixelValue{0.0f};
            for (uint32_t k=0; k<(vmm.getK()+VMFMixture::KernelType::ScalarType::Width::value-1)/VMFMixture::KernelType::ScalarType::Width::value; ++k)
            {
                VMFMixture::KernelType::ScalarType::BooleanType activeComponents = VMFMixture::KernelType::ScalarType::BooleanType{k < vmm.getK()/VMFMixture::KernelType::ScalarType::Width::value} || !!VMFMixture::KernelType::ScalarType{1.0}.cutoff(vmm.getK()%VMFMixture::KernelType::ScalarType::Width::value);
                pixelValue += lightpmm::sum(lightpmm::ifthen(activeComponents, distances[k], VMFMixture::KernelType::ScalarType{0.0f})*softAssignments.getSoftAssignments(k));
            }

            distanceBitmap->setPixel(p, Spectrum{pixelValue});
        }
    }

    return distanceBitmap;
}

ref<Bitmap> PathGuidingUtilities::renderVMFPDF(const VMFMixture &vmm, size_t component)
{
    ref<Bitmap> pdfBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    if (component >= vmm.getK())
    {
        SLog(EWarn, "rendering vMF PDF of out of range component");
        pdfBitmap->clear();
        return pdfBitmap;
    }

    Point2i p;
    for (p.y=0; p.y<renderSize.y; p.y++)
        for (p.x=0; p.x<renderSize.x; p.x++)
            pdfBitmap->setPixel(p, Spectrum{vmm.pdfK(component, sphericalPlanePointToDirVector(p))});

    return pdfBitmap;
}

ref<Bitmap> PathGuidingUtilities::renderVMMPDFColored(const VMFMixture &vmm, const std::array<Spectrum, VMFMixture::MaxK::value> &colors)
{
    ref<Bitmap> pdfBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    for (int x=0; x<renderSize.x; x++)
    {
        for (int y=0; y<renderSize.y; y++)
        {
            const Point2i p{x,y};
            const Vector dir = sphericalPlanePointToDirVector(p);

            Spectrum pixelValue(0.0f);
            for (size_t i=0; i<vmm.getK(); ++i)
                pixelValue += colors[i]*vmm.pdfK(i, dir)*vmm.weightK(i);

            pdfBitmap->setPixel(p, pixelValue);
        }
    }

    Spectrum vMFCenterHighlight(0.0f);
    vMFCenterHighlight[0] = 1.0f;
    Spectrum vMFCenterData(0.0f);

    for (unsigned int i=0; i<vmm.getK(); i++)
    {
        vMFCenterData[1] = vmm.getComponent(i/VMFMixture::KernelType::ScalarType::Width::value).getWeight(i%VMFMixture::KernelType::ScalarType::Width::value);
        vMFCenterData[2] = vmm.getComponent(i/VMFMixture::KernelType::ScalarType::Width::value).getKappa(i%VMFMixture::KernelType::ScalarType::Width::value);

        const Vector dir = vmm.getComponent(i/VMFMixture::KernelType::ScalarType::Width::value).getMu(i%VMFMixture::KernelType::ScalarType::Width::value);

        const Point2i p = dirVectorToSphericalPlanePoint(dir);

        pdfBitmap->setPixel(p, pdfBitmap->getPixel(p)*vMFCenterHighlight+vMFCenterData);
    }

    return pdfBitmap;
}

float PathGuidingUtilities::computeLogLikelihood(const VMFMixture &distribution, const std::vector<PathGuidingSampleData>& samples)
{
    float logLikelihood{0.0f};
    float totalWeight{0.0f};

    for (const PathGuidingSampleData& sample : samples)
    {
        totalWeight   += sample.weight;
        logLikelihood += sample.weight*log(distribution.pdf(sample.direction));
    }

    return logLikelihood/totalWeight;
}

float PathGuidingUtilities::computePearsonChiSquaredDivergence(const VMFMixture &distribution, const std::vector<PathGuidingSampleData>& samples)
{
    float divergence{0.0f};
    float totalWeight{0.0f};

    for (const PathGuidingSampleData& sample : samples)
    {
        totalWeight += sample.weight;
        divergence  += sample.weight*sample.weight*sample.pdf/distribution.pdf(sample.direction);
    }

    divergence *= static_cast<float>(samples.size());
    divergence /= totalWeight*totalWeight;
    divergence -= 1.0f;

    return divergence;
}

void PathGuidingUtilities::exportGuidingTreeASCII(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename)
{
    std::ofstream oss(filename.c_str());
    oss << guidingTree.getBounds().toString() << endl;

    const Guiding::BSPTree& tree = guidingTree.getTree();
    const size_t numRegions = tree.size();

    struct StackItem
    {
        uint32_t id;
        uint32_t depth;
        StackItem() = default;
        explicit StackItem(uint32_t id, uint32_t depth) : id{id}, depth{depth} {}
    };

    std::stack<StackItem> nodesToVisit {{StackItem{0,0}}};

    char buffer[12];
    char blanks[64];
    std::fill_n(blanks, sizeof(blanks)-1, ' ');
    blanks[sizeof(blanks)-1] = 0;

    char stars[64];
    std::fill_n(stars, sizeof(stars)-1, '*');
    stars[sizeof(stars)-1] = 0;

    while (!nodesToVisit.empty())
    {
        StackItem current = nodesToVisit.top();
        nodesToVisit.pop();

        int buffLen = snprintf(buffer, sizeof(buffer), "%d ", current.id);
        if (buffLen < 0)
        {
            buffer[0] = 0;
            buffLen = 0;
        }
        else if (buffLen >= static_cast<int>(sizeof(buffer)))
        {
            buffLen = sizeof(buffer)-1;
            buffer[sizeof(buffer)-1] = 0;
        }

        const uint32_t numStars = std::min(32U, current.depth);

        // id and depth
        oss << blanks+sizeof(blanks)-12+buffLen << buffer;
        oss << stars+sizeof(stars)-1-numStars << blanks+sizeof(blanks)-34+numStars;

        const Guiding::BSPTree::Node& currentNode = tree[current.id];

        if (currentNode.isInnerNode())
        {
            // split dimension and split position
            oss << ((currentNode.splitDim == 0) ? "x " : (currentNode.splitDim == 1) ? "y " : "z ");
            oss << currentNode.splitPos;
            nodesToVisit.emplace(currentNode.childIdx+1, current.depth+1);
            nodesToVisit.emplace(currentNode.childIdx, current.depth+1);
        }
        else {
            oss << "-> " << currentNode.childIdx;
        }

        oss << '\n';
    }

    oss.close();

    SLog(EInfo, "GuidingField spatial structure of %zu regions written to %s", numRegions, filename.c_str());
}

void PathGuidingUtilities::exportGuidingTreeOBJ(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename)
{
    std::ofstream oss(filename.c_str());
    oss << "o GuidingFieldBounds" << endl;

    const size_t numRegions = guidingTree.getNumRegions();

    for (size_t i=0; i<numRegions; i++)
    {
        const AABB& aabb = guidingTree.getRegion(i).bound;
        for (unsigned int i=0; i<8; i++)
        {
            const Point corner = aabb.getCorner(i);

            oss << "v " << corner.x << " " << corner.y << " " << corner.z << endl;
        }
    }

    // generate box faces
    for (size_t i=0; i<numRegions; i++)
    {
        oss << "f " << 8*i+0+1 << " " << 8*i+1+1 << " " << 8*i+3+1 << " " << 8*i+2+1 << endl;
        oss << "f " << 8*i+1+1 << " " << 8*i+5+1 << " " << 8*i+7+1 << " " << 8*i+3+1 << endl;
        oss << "f " << 8*i+5+1 << " " << 8*i+4+1 << " " << 8*i+6+1 << " " << 8*i+7+1 << endl;
        oss << "f " << 8*i+4+1 << " " << 8*i+0+1 << " " << 8*i+2+1 << " " << 8*i+6+1 << endl;
        oss << "f " << 8*i+2+1 << " " << 8*i+3+1 << " " << 8*i+7+1 << " " << 8*i+6+1 << endl;
        oss << "f " << 8*i+0+1 << " " << 8*i+1+1 << " " << 8*i+5+1 << " " << 8*i+4+1 << endl;
    }
    oss.close();

    SLog(EInfo, "GuidingField bounds of %zu regions written to %s", numRegions, filename.c_str());
}

void PathGuidingUtilities::exportGuidingTreeSampleBoundsOBJ(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename)
{
    std::ofstream oss(filename.c_str());
    oss << "o GuidingFieldSampleBounds" << endl;

    const size_t numRegions = guidingTree.getNumRegions();

    for (size_t i=0; i<numRegions; i++)
    {
        const AABB& aabb = guidingTree.getRegion(i).statsSinceLastSplit.sampleBound;
        for (unsigned int i=0; i<8; i++)
        {
            const Point corner = aabb.getCorner(i);

            oss << "v " << corner.x << " " << corner.y << " " << corner.z << endl;
        }
    }

    // generate box faces
    for (size_t i=0; i<numRegions; i++)
    {
        oss << "f " << 8*i+0+1 << " " << 8*i+1+1 << " " << 8*i+3+1 << " " << 8*i+2+1 << endl;
        oss << "f " << 8*i+1+1 << " " << 8*i+5+1 << " " << 8*i+7+1 << " " << 8*i+3+1 << endl;
        oss << "f " << 8*i+5+1 << " " << 8*i+4+1 << " " << 8*i+6+1 << " " << 8*i+7+1 << endl;
        oss << "f " << 8*i+4+1 << " " << 8*i+0+1 << " " << 8*i+2+1 << " " << 8*i+6+1 << endl;
        oss << "f " << 8*i+2+1 << " " << 8*i+3+1 << " " << 8*i+7+1 << " " << 8*i+6+1 << endl;
        oss << "f " << 8*i+0+1 << " " << 8*i+1+1 << " " << 8*i+5+1 << " " << 8*i+4+1 << endl;
    }
    oss.close();

    SLog(EInfo, "GuidingField bounds of %zu regions written to %s", numRegions, filename.c_str());
}

void PathGuidingUtilities::exportGuidingTreeCentersOBJ(const GuidingFieldType::GuidingTreeType& guidingTree, const std::string &filename)
{
    std::ofstream oss(filename.c_str());
    oss << "o GuidingFieldCenters" << endl;

    const size_t numRegions = guidingTree.getNumRegions();

    for (size_t i=0; i<numRegions; i++)
    {
        const AABB& aabb = guidingTree.getRegion(i).bound;
        const Point center = aabb.getCenter();
        oss << "v " << center.x << " " << center.y << " " << center.z << endl;
    }

    oss.close();

    SLog(EInfo, "GuidingField centers of %zu regions written to %s", numRegions, filename.c_str());
}

PPG::STree PathGuidingUtilities::loadPPGSDTree(const std::string &filename)
{
    PPG::STree sTree;
    PPG::BlobReader reader(filename);

    if (reader.isValid())
    {
        SLog(EInfo, "reading SDTree from file %s", filename.c_str());

        //unused
        {
            Matrix4x4 camera;
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    reader >> camera(i, j);
                }
            }
        }

        PPG::DTree dTree;
        while (dTree.read(reader))
        {
            sTree.aabb.expandBy(dTree.aabb());
            sTree.dTrees.emplace_back(std::move(dTree));
            dTree = PPG::DTree{};
        }

        SLog(EInfo, sTree.toString().c_str());
    }
    else
        SLog(EWarn, "failed to read SDTree from file %s", filename.c_str());

    return sTree;
}

ref<Bitmap> PathGuidingUtilities::renderDTreePDF(const PPG::DTree &dTree)
{
    ref<Bitmap> pdfBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat32, renderSize);
    Point2i p;
    for (p.y=0; p.y<renderSize.y; p.y++)
        for (p.x=0; p.x<renderSize.x; p.x++)
            pdfBitmap->setPixel(p, Spectrum(dTree.pdf(sphericalPlanePointToDirVector(p))));

    return pdfBitmap;
}

void PathGuidingUtilities::exportSTreeOBJ(const PPG::STree &sTree, const std::string &filename)
{
    std::ofstream oss(filename.c_str());
    oss << "o STreeBounds" << endl;

    const size_t numRegions = sTree.dTrees.size();

    for (size_t i=0; i<numRegions; i++)
    {
        const AABB& aabb = sTree.dTrees.at(i).aabb();
        for (unsigned int i=0; i<8; i++)
        {
            const Point corner = aabb.getCorner(i);

            oss << "v " << corner.x << " " << corner.y << " " << corner.z << endl;
        }
    }

    // generate box faces
    for (size_t i=0; i<numRegions; i++)
    {
        oss << "f " << 8*i+0+1 << " " << 8*i+1+1 << " " << 8*i+3+1 << " " << 8*i+2+1 << endl;
        oss << "f " << 8*i+1+1 << " " << 8*i+5+1 << " " << 8*i+7+1 << " " << 8*i+3+1 << endl;
        oss << "f " << 8*i+5+1 << " " << 8*i+4+1 << " " << 8*i+6+1 << " " << 8*i+7+1 << endl;
        oss << "f " << 8*i+4+1 << " " << 8*i+0+1 << " " << 8*i+2+1 << " " << 8*i+6+1 << endl;
        oss << "f " << 8*i+2+1 << " " << 8*i+3+1 << " " << 8*i+7+1 << " " << 8*i+6+1 << endl;
        oss << "f " << 8*i+0+1 << " " << 8*i+1+1 << " " << 8*i+5+1 << " " << 8*i+4+1 << endl;
    }
    oss.close();

    SLog(EInfo, "STree bounds of %zu DTrees written to %s", numRegions, filename.c_str());
}

GUIDING_NAMESPACE_END
