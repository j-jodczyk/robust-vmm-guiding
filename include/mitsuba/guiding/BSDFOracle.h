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

#pragma once
 
#include "guiding.h"
#include <pmm/VMFKernel.h>

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/scene.h>

#include <unordered_map>

GUIDING_NAMESPACE_BEGIN

//#define DIFFUSE_KAPPA 2.133
#define DIFFUSE_KAPPA 2.18853f

template<typename PMM>
struct BSDFRepresentation{
    typedef std::integral_constant<int, 4> MaxLobes;
    typedef std::integral_constant<int, 4096> ThetaDiscretizationSteps;

    static constexpr Float m_minTheta = 0.0f;
    static constexpr Float m_maxTheta = M_PI/2.0f - (1e-3f);
    static constexpr Float m_thetaStepsize = (m_maxTheta-m_minTheta)/Float(ThetaDiscretizationSteps::value-1);

    const BSDF* m_bsdf;

    Float m_roughThetasDiffuse[ThetaDiscretizationSteps::value];
    Float m_roughThetasGlossy[ThetaDiscretizationSteps::value];

    PMM m_pmm;
    PMM m_fittedPMMs[ThetaDiscretizationSteps::value];

    BSDFRepresentation() = default;

    BSDFRepresentation(const BSDF* bsdf)
        : m_bsdf(bsdf)
    {
        Intersection its;
        its.uv = Point2(0.5, 0.5);
        its.dudx = its.dudy = 0.5;
        its.dvdx = its.dvdy = 0.5;
        its.wi = Vector(0.0, 0.0, 1.0);
        BSDFSamplingRecord bRec(its, Vector(0.0f, 0.0f, 1.0f));

        switch (m_bsdf->getModel())
        {
        case BSDF::EBSDFModel::EMPhong:
        {
            const Float glossySamplingRate = bsdf->getGlossySamplingRate(bRec);
            const Float roughness = bsdf->getRoughness(its,0);
            const Float exponent = (2.0/(roughness*roughness))-2.0;
            m_pmm.setK(2);
            m_pmm.m_comps[0].setWeight(0,glossySamplingRate);
            m_pmm.m_comps[0].setKappa(0,exponent);
            m_pmm.m_comps[0].setMu(0,Vector3(0.0,0.0,1.0));

            m_pmm.m_comps[0].setWeight(1,1.0-glossySamplingRate);
            m_pmm.m_comps[0].setKappa(1,DIFFUSE_KAPPA);
            m_pmm.m_comps[0].setMu(1,Vector3(0.0,0.0,1.0));
            break;
        }

        case BSDF::EBSDFModel::EMRoughConductor:
            fitGlossyBSDFLobes(bsdf);
            break;

        case BSDF::EBSDFModel::EMRoughPlastic:
        case BSDF::EBSDFModel::EMTorrance:
            fitDiffuseGlossyBSDFLobes(bsdf);
            break;

        case BSDF::EBSDFModel::EMOther:
            SLog(EWarn, "Approximating BSDF %s as diffuse.", bsdf->getClass()->getName().c_str());
        case BSDF::EBSDFModel::EMPlastic:
        case BSDF::EBSDFModel::EMSmoothDiffuse:
            m_pmm.setK(1);
            m_pmm.m_comps[0].setWeight(0,1.0);
            m_pmm.m_comps[0].setKappa(0,DIFFUSE_KAPPA);
            m_pmm.m_comps[0].setMu(0,Vector3(0.0,0.0,1.0));
            break;

        case BSDF::EBSDFModel::EMConductor:
        case BSDF::EBSDFModel::EMDielectric:
        case BSDF::EBSDFModel::EMRoughDielectric:
        default:
            m_pmm.setK(0);
            break;
        }
    }

    PMM getBSDFPMM(const Intersection &its, const Normal n, const Vector wi) const
    {
        PMM samplingPMM = m_pmm;

        Float thetaSign = math::signum(Frame::cosTheta(its.wi));

        Vector wo = reflect(wi, n);

        Point2 thetaPhi = toSphericalCoordinates(its.wi);
        thetaPhi[0] *= thetaSign;
        int thetaIdx = std::max(0,std::min(int(thetaPhi[0]/m_thetaStepsize), ThetaDiscretizationSteps::value-1));

        switch (m_bsdf->getModel())
        {
        case BSDF::EBSDFModel::EMPhong:
            samplingPMM.m_comps[0].setMu(0,wo);
            samplingPMM.m_comps[0].setMu(1,n);
            break;

        case BSDF::EBSDFModel::EMRoughConductor:
            samplingPMM = m_fittedPMMs[thetaIdx];
            samplingPMM.m_comps[0].setMu(0,its.toWorld(sphericalDirection(thetaSign*m_roughThetasGlossy[thetaIdx], thetaPhi[1]+M_PI)));
            break;

        case BSDF::EBSDFModel::EMRoughPlastic:
        case BSDF::EBSDFModel::EMTorrance:
            samplingPMM = m_fittedPMMs[thetaIdx];
            samplingPMM.m_comps[0].setMu(0,its.toWorld(sphericalDirection(m_roughThetasGlossy[thetaIdx]*thetaSign, thetaPhi[1]+M_PI)));
            //samplingPMM.m_comps[0].setMu(1,its.toWorld(sphericalDirection(m_roughThetasDiffuse[thetaIdx]*thetaSign, thetaPhi[1]+M_PI)));
            samplingPMM.m_comps[0].setMu(1,n);
            break;

        case BSDF::EBSDFModel::EMOther:
        case BSDF::EBSDFModel::EMPlastic:
        case BSDF::EBSDFModel::EMSmoothDiffuse:
            samplingPMM.m_comps[0].setMu(0,n);
            break;

        case BSDF::EBSDFModel::EMConductor:
        case BSDF::EBSDFModel::EMDielectric:
        case BSDF::EBSDFModel::EMRoughDielectric:
        default:
            SLog(EError, "Product sampling is not supported for BSDF %s.", m_bsdf->getClass()->getName().c_str());
            break;
        }

        return samplingPMM;
    }

    void vmfWeightedMLE(const int N, const Vector3f *wos, const Float *weights, Float &theta, Float &kappa)
    {
        Vector3f vecR(0.0,0.0,0.0);
        Float sumW = 0.0;
        for (int i=0; i<N;i++){
            vecR += wos[i]*weights[i]; 
            sumW += weights[i];
        }

        vecR /= sumW;
        Float meanCosine = vecR.length();
        Vector3f mean = vecR/meanCosine;

        theta = std::acos(Frame::cosTheta(mean));

        const Float maxMeanCosine = lightpmm::kappaToMeanCosine(32768.0f);

        if (meanCosine > maxMeanCosine)
        {
            //SLog(EWarn, "clipped meanCosine %f to %f", meanCosine, maxMeanCosine);
            meanCosine = maxMeanCosine;
        }

        kappa = 0.0f;
        if (meanCosine>0.0f){
            kappa = ((meanCosine*3.0f)-(meanCosine*meanCosine*meanCosine))/(1.0f-meanCosine*meanCosine);
        }
    }

    void fitGlossyBSDFLobes(const BSDF* bsdf)
    {
        const int N = 2048*1;
        
        Properties props("independent");
        props.setInteger("sampleCount", 10*N);

#ifdef MTS_OPENMP
#pragma omp parallel
#endif
        {
            Vector3f wos[N];
            Float weights[N];

            ref<Sampler> sampler = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
            sampler->configure();

            Intersection itsBSDF;

#ifdef MTS_OPENMP
#pragma omp for
#endif
            for(int t = 0; t<ThetaDiscretizationSteps::value; t++)
            {
                sampler->generate(Point2i(0, 0));
                Float thetaIn = t*m_thetaStepsize;

                itsBSDF.wi = sphericalDirection(thetaIn,0.0);
                BSDFSamplingRecord bRec(itsBSDF, sampler);

                for (int n=0; n<N; ++n)
                {
                    Float bsdfWeight;
                    do
                    {
                        bsdfWeight = bsdf->sample(bRec, sampler->next2D()).average();
                        sampler->advance();
                    } while (bsdfWeight == 0.0);

                    wos[n] = bRec.wo;
                    weights[n] = bsdfWeight;
                }

                Float thetaOut = 0.0;
                Float kappa = DIFFUSE_KAPPA;
                vmfWeightedMLE(N, wos, weights, thetaOut, kappa);
                m_roughThetasGlossy[t] = thetaOut;
                //m_roughKappasGlossy[t] = kappa;

                Vector3 muOut = sphericalDirection(thetaOut, 0.0f);
                m_fittedPMMs[t].setK(1);
                m_fittedPMMs[t].m_comps[0].setWeight(0,1.0f);
                m_fittedPMMs[t].m_comps[0].setMu(0,muOut);
                m_fittedPMMs[t].m_comps[0].setKappa(0,kappa);
            }
        }
    }

    void fitDiffuseGlossyBSDFLobes(const BSDF* bsdf)
    {
        const int N = 2048*1;
        
        Properties props("independent");
        props.setInteger("sampleCount", 10*N);

#ifdef MTS_OPENMP
#pragma omp parallel
#endif
        {
            Vector3f wos[N];
            Float weights[N];

            ref<Sampler> sampler = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
            sampler->configure();

            Intersection itsBSDF;

#ifdef MTS_OPENMP
#pragma omp for
#endif
            for(int t = 0; t< ThetaDiscretizationSteps::value;t++)
            {
                sampler->generate(Point2i(0, 0));
                Float thetaIn = t*m_thetaStepsize;

                Float diffuseAlbedo = 0.0f;
                Float glossyAlbedo  = 0.0f;

                itsBSDF.wi = sphericalDirection(thetaIn,0.0);

                BSDFSamplingRecord bRecDiffuse(itsBSDF, sampler);
                bRecDiffuse.typeMask = BSDF::EDiffuseReflection;


                Float thetaOutDiffuse = 0.0;
                Float kappaDiffuse = DIFFUSE_KAPPA;
                Float diffuseReflectance = bsdf->getDiffuseReflectance(itsBSDF).average();

                if (diffuseReflectance > 0) {
                    size_t nDiffuseSamples = 0;

                    for (int n=0; n<N; ++n)
                    {
                        Float bsdfWeight;
                        do
                        {
                            bsdfWeight = bsdf->sample(bRecDiffuse, sampler->next2D()).average();
                            sampler->advance();
                            ++nDiffuseSamples;
                        } while (bsdfWeight == 0.0);

                        diffuseAlbedo += bsdfWeight;
                        wos[n] = bRecDiffuse.wo;
                        weights[n] = bsdfWeight;
                    }

                    diffuseAlbedo /= Float(nDiffuseSamples);

                    vmfWeightedMLE(N, wos, weights, thetaOutDiffuse, kappaDiffuse);
                }
                m_roughThetasDiffuse[t] = thetaOutDiffuse;
                kappaDiffuse = DIFFUSE_KAPPA;
                //m_roughKappasDiffuse[t] = kappaDiffuse;


                sampler->generate(Point2i(0, 0));

                BSDFSamplingRecord bRecGlossy(itsBSDF, sampler);
                bRecGlossy.typeMask = BSDF::EGlossyReflection;

                Float thetaOutGlossy = 0.0;
                Float kappaGlossy = DIFFUSE_KAPPA;

                Float specularReflectance = bsdf->getSpecularReflectance(itsBSDF).average();
                if (specularReflectance > 0)
                {
                    size_t nGlossySamples = 0;

                    for (int n=0; n<N; ++n)
                    {
                        Float bsdfWeight;
                        do
                        {
                            bsdfWeight = bsdf->sample(bRecGlossy, sampler->next2D()).average();
                            sampler->advance();
                            ++nGlossySamples;
                        } while (bsdfWeight == 0.0);

                        glossyAlbedo += bsdfWeight;
                        wos[n] = bRecGlossy.wo;
                        weights[n] = bsdfWeight;
                    }

                    glossyAlbedo /= Float(nGlossySamples);

                    vmfWeightedMLE(N, wos, weights, thetaOutGlossy, kappaGlossy);
                }
                SAssert(std::isfinite(thetaOutGlossy));
                SAssert(!std::isnan(thetaOutGlossy));
                m_roughThetasGlossy[t] = thetaOutGlossy;

                Float glossyWeight = glossyAlbedo/(glossyAlbedo+diffuseAlbedo);

                Vector3 muOutDiffuse = sphericalDirection(thetaOutDiffuse, 0.0f);
                Vector3 muOutGlossy = sphericalDirection(thetaOutGlossy, 0.0f);
                m_fittedPMMs[t].setK(2);
                m_fittedPMMs[t].m_comps[0].setWeight(0,glossyWeight);
                m_fittedPMMs[t].m_comps[0].setMu(0,muOutGlossy);
                m_fittedPMMs[t].m_comps[0].setKappa(0,kappaGlossy);

                m_fittedPMMs[t].m_comps[0].setWeight(1,1.0f-glossyWeight);
                m_fittedPMMs[t].m_comps[0].setMu(1,muOutDiffuse);
                m_fittedPMMs[t].m_comps[0].setKappa(1,kappaDiffuse);
            }
        }
    }

    inline Vector3f reflect(Vector3f wi, Vector3f n) const {
        return n*(2.0*dot(wi,n))-wi;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "BSDFRepresentation[" << endl
            << "  pmm = " <<  m_pmm.toString() << "," << endl
            << "  model = " << m_bsdf->getModel() << "," << endl
            << "  bsdf = " << indent(m_bsdf->toString()) << endl
            << "]";
        return oss.str();
    }
};

template<typename PMM>
class BSDFOracleVMF {
public:
    typedef BSDFRepresentation<PMM> BSDFRepresentationType;
    
    static void initBSDFRepresentations(const Scene* scene)
    {
#ifdef MTS_OPENMP
        ref<Scheduler> scheduler = Scheduler::getInstance();
        size_t nCores = scheduler->getCoreCount();

        Thread::initializeOpenMP(nCores);
#endif

        const ref_vector<ConfigurableObject>& objects = scene->getReferencedObjects();
        for(const ref<ConfigurableObject>& object : objects)
            if (const BSDF* bsdf = dynamic_cast<const BSDF*>(object.get()))
                m_bsdfRepresentaions.emplace(bsdf, bsdf);
    }

    static const BSDFRepresentationType *getBSDFRepresentation(const BSDF* bsdf)
    {
        Guiding_Assert(bsdf);
        const auto& bsdfRepresentation = m_bsdfRepresentaions.find(bsdf);
        if (EXPECT_NOT_TAKEN(bsdfRepresentation == m_bsdfRepresentaions.end()))
        {
            SLog(EWarn,"BSDF product representation not found: %s",bsdf->toString().c_str());
            SLog(EError,"make sure that all used BSDFs are registered as immediate children of the <scene> node and then referenced where needed.");
        }
        return &(*bsdfRepresentation).second;
    }

private:
    static std::unordered_map<const BSDF*, BSDFRepresentationType> m_bsdfRepresentaions;
};


GUIDING_NAMESPACE_END
