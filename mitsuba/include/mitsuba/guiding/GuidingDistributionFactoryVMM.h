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

#ifndef GUIDINGDISTRIBUTIONFACTORYVMM_H
#define GUIDINGDISTRIBUTIONFACTORYVMM_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/math.h>

#include <pmm/pmm.h>
#include <pmm/ParametricMixtureModel.h>
#include <pmm/VMFKernel.h>
#include <pmm/VMMFactory.h>

#include <mitsuba/guiding/guiding.h>
#include <mitsuba/guiding/GuidingFieldFactoryFwd.h>
#include <mitsuba/guiding/GuidedBSDF.h>
#ifdef GUIDING_USE_BSDF_PRODUCT
#include <mitsuba/guiding/BSDFOracle.h>
#endif

#include <string>
#include <sstream>

GUIDING_NAMESPACE_BEGIN

template<typename TScalar, uint32_t NKernels>
class GuidingDistributionFactory<lightpmm::<lightpmm::VMFKernel<TScalar>, NKernels>>
{
private:
    using TDistribution = lightpmm::ParametricMixtureModel<lightpmm::VMFKernel<TScalar>, NKernels>;
    using TDistributionFactory = lightpmm::VMMFactory<TDistribution>;
    using TVMFKernel = lightpmm::VMFKernel<TScalar>;
    using TVMM = lightpmm::ParametricMixtureModel<lightpmm::VMFKernel<TScalar>, NKernels>;
#ifdef GUIDING_USE_BSDF_PRODUCT
    //instantiation of the template without explicit instantiation of its static member variable result in a linker error
    //therefore, this is disabled by default
    using TBSDFVMM = lightpmm::ParametricMixtureModel<lightpmm::VMFKernel<TScalar>, 1>;
    using TBSDFOracle = Guiding::BSDFOracleVMF<TBSDFVMM>;
#endif
public:

    GuidingDistributionFactory(const Properties& props=Properties{})
    {
        lightpmm::VMMFactoryProperties vmmFactoryProps;
        vmmFactoryProps.minItr                    = props.getSize( "vmmFactory.minItr",                    vmmFactoryProps.minItr);
        vmmFactoryProps.maxItr                    = props.getSize( "vmmFactory.maxItr",                    vmmFactoryProps.maxItr);
        vmmFactoryProps.relLogLikelihoodThreshold = props.getFloat("vmmFactory.relLogLikelihoodThreshold", vmmFactoryProps.relLogLikelihoodThreshold);
        vmmFactoryProps.numInitialComponents      = props.getSize( "vmmFactory.numInitialComponents",      vmmFactoryProps.numInitialComponents);
        vmmFactoryProps.initKappa                 = props.getFloat("vmmFactory.initKappa",                 vmmFactoryProps.initKappa);
        vmmFactoryProps.maxKappa                  = props.getFloat("vmmFactory.maxKappa",                  vmmFactoryProps.maxKappa);
        vmmFactoryProps.vPrior                    = props.getFloat("vmmFactory.vPrior",                    vmmFactoryProps.vPrior);
        vmmFactoryProps.rPrior                    = props.getFloat("vmmFactory.rPrior",                    vmmFactoryProps.rPrior);
        vmmFactoryProps.rPriorWeight              = props.getFloat("vmmFactory.rPriorWeight",              vmmFactoryProps.rPriorWeight);

        m_vmmFactory = TDistributionFactory{vmmFactoryProps};

        m_decayOnSpatialSplit = props.getFloat("decayOnSpatialSplit", 0.25f);
    }

    GuidingDistributionFactory(Stream* stream)
    {
        size_t factoryStreamSize = stream->readSize();
        char factoryStreamString[factoryStreamSize];
        stream->read(factoryStreamString, factoryStreamSize);
        std::istringstream factorySerializationStream{std::string{factoryStreamString, factoryStreamSize}};
        m_vmmFactory.deserialize(factorySerializationStream);
        m_decayOnSpatialSplit = stream->readSingle();
    }

    void serialize(Stream* stream) const
    {
        std::ostringstream factorySerializationStream;
        m_vmmFactory.serialize(factorySerializationStream);
        const std::string factoryStreamString = factorySerializationStream.str();
        const size_t factoryStreamSize = factoryStreamString.size();
        stream->writeSize(factoryStreamSize);
        stream->write(factoryStreamString.data(), factoryStreamSize);
        stream->writeSingle(m_decayOnSpatialSplit);
    }

    void initialize(TDistribution& distribution) const {
        m_vmmFactory.initialize(distribution);
    }

    ///to be called on spatial splits (only once per tree update, even for repeated splits)
    void onSpatialSplit(TDistribution& distribution) const
    {
        distribution.m_numSamples   *= m_decayOnSpatialSplit;
        distribution.m_sampleWeight *= m_decayOnSpatialSplit;
    }

    template<typename TSampleRange>
    void fit(TDistribution& distribution, TSampleRange& samples) const {
        m_vmmFactory.fit(samples.begin(), samples.end(), distribution, false);
    }

    template<typename TSampleRange>
    void updateFit(TDistribution& distribution, TSampleRange& samples) const {
        m_vmmFactory.updateFit(samples.begin(), samples.end(), distribution);
    }

    template<typename TStatistics>
    void prepareGuidedBSDFForSampling(GuidedBSDF<TDistribution, TStatistics>& gBSDF, const Intersection& its, const Vector wiWorld) const
    {
        if (EXPECT_NOT_TAKEN(gBSDF.m_bsdf->hasComponent(BSDF::ETransmission&~BSDF::ENull)))
        {
            gBSDF.m_guidingData.clear();
            return;
        }

        gBSDF.m_guidingData.m_distributions[0] = gBSDF.m_guidingData.m_liDistribution;
        gBSDF.m_guidingData.m_weights[0] = 1.0f;
        gBSDF.m_guidingData.m_numDistributions = 1;
        gBSDF.m_guidingData.m_productIntegral = 1.0f;

        const Vector surfaceNormalTowardsWi = its.shFrame.n*math::signum(Frame::cosTheta(its.wi));

        if (gBSDF.m_useBSDFProduct)
        {
#ifdef GUIDING_USE_BSDF_PRODUCT
            const typename TBSDFOracle::BSDFRepresentationType* bsdfRep = TBSDFOracle::getBSDFRepresentation(gBSDF.m_bsdf);

            const TBSDFVMM bsdfPMM = bsdfRep->getBSDFPMM(its, surfaceNormalTowardsWi, wiWorld);

            gBSDF.m_guidingData.m_productIntegral = 0.0f;
            unsigned int bsdfLobeCount = 0;
            for(unsigned int i = 0; i < bsdfPMM.getK();i++){
                if(bsdfPMM.m_comps[0].getWeight(i) > 0.0f){
                    gBSDF.m_guidingData.m_distributions[bsdfLobeCount] = gBSDF.m_guidingData.m_liDistribution;
                    const TVMFKernel bsdfLobe = bsdfPMM.m_comps[0].extract(i);

                    Float productWeight = gBSDF.m_guidingData.m_distributions[bsdfLobeCount].product(bsdfLobe);
                    gBSDF.m_guidingData.m_weights[bsdfLobeCount] = productWeight;
                    gBSDF.m_guidingData.m_productIntegral += productWeight;
                    bsdfLobeCount++;
                }
            }
            for(unsigned int i = 0; i < bsdfLobeCount;i++){
                gBSDF.m_guidingData.m_weights[i] /= gBSDF.m_guidingData.m_productIntegral;
            }
            gBSDF.m_guidingData.m_numDistributions = bsdfLobeCount;
#else
            SLog(EError, "bsdf-product support was not enabled during compilation");
#endif
        }
        else if (gBSDF.m_useCosineProduct)
        {
            const TVMFKernel vmfCosineLobe(1.0f, 2.18853f, surfaceNormalTowardsWi);
            gBSDF.m_guidingData.m_productIntegral = gBSDF.m_guidingData.m_distributions[0].product(vmfCosineLobe);
        }
    }

    FINLINE bool checkValidity(const TDistribution& distribution) const {
#ifdef GUIDING_NO_VALIDITY_CHECKS
        return true;
#else
        return distribution.valid(m_vmmFactory.computeMinComponentWeight(distribution.getK()));
#endif
    }

    std::string toString() const
    {
        std::ostringstream oss;
        oss.setf(std::ios_base::boolalpha);

        oss << "GuidingDistributionFactory<VMF Mixture> [\n"
            << "  vmmFactory = " << m_vmmFactory.toString() << '\n'
            << "  decayOnSpatialSplits = " << m_decayOnSpatialSplit << '\n'
            << ']';
        return oss.str();
    }

    TDistributionFactory m_vmmFactory;
    float m_decayOnSpatialSplit {0.25f};
};

GUIDING_NAMESPACE_END

#endif // GUIDINGDISTRIBUTIONFACTORYVMM_H
