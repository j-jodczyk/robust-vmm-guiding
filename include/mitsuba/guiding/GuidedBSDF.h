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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/platform.h>

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/shape.h>

#include "guiding.h"
#include "GuidingRegion.h"

GUIDING_NAMESPACE_BEGIN

template<typename TDistribution, typename TStatistics>
class GuidedBSDF : public BSDF {
public:
    using TGuidingRegion = GuidingRegion<TDistribution, TStatistics>;

    GuidedBSDF(Float bsdfSamplingProbability, bool useCosineProduct=true, bool useBSDFProduct=false)
        : BSDF(Properties()),
          m_bsdfSamplingProbability(bsdfSamplingProbability),
          m_useCosineProduct(useCosineProduct),
          m_useBSDFProduct(useBSDFProduct)
    {
	}

    void configure() override
    {
        Guiding_Assert(m_bsdf);

        // check if we can use guiding
        const bool useGuiding = m_guidingData.valid() && bsdfTypeAllowsGuiding();

        m_glossyComponentId = -1;
        m_diffuseComponentId = -1;

        m_combinedType = m_bsdf->getType();
        m_components.resize(m_bsdf->getComponentCount());
        for (int i=0; i<m_bsdf->getComponentCount(); ++i)
        {
            m_components[i] = m_bsdf->getType(i);

            if (m_components[i]&BSDF::EDiffuse)
                m_diffuseComponentId = i;
            else if (m_components[i]&BSDF::EGlossy)
                m_glossyComponentId = i;
        }

        if (useGuiding)
            m_combinedType |= BSDF::EGuiding;
    }

    Spectrum eval(const BSDFSamplingRecord & bRec, EMeasure measure = ESolidAngle) const override
    {
        return m_bsdf->eval(bRec, measure);
	}

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const override
    {
        if (!canUseGuiding()) {
            return m_bsdf->sample(bRec, sample);
        }

        Float pdf;
        return this->sample(bRec, pdf, sample);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const override
    {
        if (!canUseGuiding()) {
            //use only the BSDF for sampling
            return m_bsdf->sample(bRec, pdf, sample);
        }

        const bool hasDelta = m_bsdf->getType()&bRec.typeMask&BSDF::EDelta;
        const Float deltaSamplingRate = hasDelta ? m_bsdf->getDeltaSamplingRate(bRec) : 0.0f;
        const unsigned int bRecTypemask = bRec.typeMask;

        if (hasDelta && bRec.sampler->next1D() < deltaSamplingRate) {
            //sample the delta component
            bRec.typeMask &= BSDF::EDelta;

            Spectrum throughput = m_bsdf->sample(bRec, pdf, sample);

            pdf *= deltaSamplingRate;
            throughput /= deltaSamplingRate;

            bRec.typeMask = bRecTypemask;

            return throughput;
        }

        //from here on, sample only smooth components
        bRec.typeMask &= BSDF::ESmooth;

        // select between BSDF sampling and guiding distribution sampling
        const bool sampleBSDF = bRec.sampler->next1D() < m_bsdfSamplingProbability;

        pdf = 0.0f;
        Spectrum throughput = Spectrum(0.f);
        Float bsdfPDF, guidingPDF;

        if (!sampleBSDF)
        {
            // sample the guiding distribution
            Vector3 wo = m_guidingData.sample(sample);
            Guiding_Assert(Guiding::isValid(wo));
            bRec.wo = bRec.its.toLocal(wo);

            throughput = m_bsdf->eval(bRec);
            //for invalid samples, we can terminate early
            if (throughput.isZero())
            {
                bRec.sampledType = BSDF::EGuiding;
                return Spectrum{0.0f};
            }

            bsdfPDF = m_bsdf->pdf(bRec);
            guidingPDF = m_guidingData.pdf(wo);

            // the change in index of refraction depends upon which side has been sampled
            bRec.eta = Frame::cosTheta(bRec.wi) > 0.0f ? m_bsdf->getEta() : 1.0f/m_bsdf->getEta();

            // set the sampled component and type
            if (m_combinedType&BSDF::EGlossy && m_combinedType&BSDF::EDiffuse)
            {
                const float glossySamplingRate = m_bsdf->getGlossySamplingRate(bRec);
                bRec.typeMask &= BSDF::EGlossy;
                const Float glossyPDF = m_bsdf->pdf(bRec);
                if (glossySamplingRate*glossyPDF >= bsdfPDF*0.5f)
                    bRec.sampledComponent = m_glossyComponentId;
                else
                    bRec.sampledComponent = m_diffuseComponentId;
            }
            else if (m_combinedType&BSDF::EGlossy)
                bRec.sampledComponent = m_glossyComponentId;
            else if (m_combinedType&BSDF::EDiffuse)
                bRec.sampledComponent = m_diffuseComponentId;
            //this should never happen
            else
                SLog(EError, "BSDF has no glossy or diffuse component, yet guiding has been used for sampling.");
            bRec.sampledType = getType(bRec.sampledComponent)|BSDF::EGuiding;
        }
        else
        {
            // sample the BSDF
            throughput = m_bsdf->sample(bRec, bsdfPDF, sample);
            //for invalid samples, we can terminate early
            if (throughput.isZero())
            {
                bRec.sampledType = m_bsdf->getType();
                return Spectrum{0.0f};
            }

            Guiding_Assert(Guiding::isValid(bRec.wo));
            // correct the sample result, since it is already
            // multiplied by the inverse pdf
            throughput *= bsdfPDF;

            guidingPDF = m_guidingData.pdf(bRec.its.toWorld(bRec.wo));
        }

        //compute the one-sample MIS PDF
        pdf = guidingPDF+(bsdfPDF-guidingPDF)*m_bsdfSamplingProbability;

        if (hasDelta)
            pdf *= 1.0f-deltaSamplingRate;

        bRec.typeMask = bRecTypemask;

        if (EXPECT_NOT_TAKEN(!Guiding::isValid(pdf)))
            SLog(EWarn, "invalid PDF evaluated for direction %s in GuidedBSDF\n%s", bRec.its.toWorld(bRec.wo).toString().c_str(), toString().c_str());

        Guiding_Assert(throughput.isValid());
        Guiding_Assert(Guiding::isValid(pdf));

        return throughput/pdf;
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure = ESolidAngle) const override {
        if (!canUseGuiding()) {
            //use only the BSDF for sampling
            return m_bsdf->pdf(bRec);
        }
        BSDFSamplingRecord bRecCopy(bRec);

        // start with evaluating the pdf for smooth reflectance
        bRecCopy.typeMask &= BSDF::ESmooth;

        const Float bsdfPDFSmooth = m_bsdf->pdf(bRecCopy, measure);
        const Vector woWorld = bRec.its.toWorld(bRec.wo);
        const Float guidedPDF = m_guidingData.pdf(woWorld);
        Float pdf = guidedPDF+(bsdfPDFSmooth-guidedPDF)*m_bsdfSamplingProbability;

        // if there is a delta component, consider its contribution
        if (m_bsdf->getType()&bRec.typeMask&BSDF::EDelta) {
            const Float deltaReflectanceRate = m_bsdf->getDeltaSamplingRate(bRec);

            bRecCopy.typeMask = bRec.typeMask&BSDF::EDelta;
            const Float bsdfPDFDelta = m_bsdf->pdf(bRecCopy, measure);

            pdf += (bsdfPDFDelta-pdf)*deltaReflectanceRate;
        }

        return pdf;
    }

    FINLINE bool canUseGuiding() const {
        return (m_combinedType&BSDF::EGuiding);
    }

	/// Return a string representation
    std::string toString() const override {
		std::ostringstream oss;
        oss << "GuidedBSDF[\n"
            << "  bsdf: " << (m_bsdf ? m_bsdf->toString() : "nullptr") << '\n'
            << "  region: " << (m_region ? m_region->toString() : "nullptr") << '\n'
            << "  guidingData: " << m_guidingData.toString() << '\n'
            << "  bsdfSamplingProbability: " << m_bsdfSamplingProbability << '\n'
            << "  useCosineProduct: " << m_useCosineProduct << '\n'
            << "  useBSDFProduct: " << m_useBSDFProduct << '\n'
            << ']';
		return oss.str();
	}

    Spectrum estimateReflectedRadiance(const BSDFSamplingRecord &bRec, const Intersection &its) const
    {
        Guiding_Assert(m_region);
        Guiding_Assert(m_guidingData.valid());

        const bool hasGlossy  = m_bsdf->hasComponent(BSDF::EGlossy);
        const bool hasDiffuse = m_bsdf->hasComponent(BSDF::EDiffuse);
        Spectrum averageReflectance{0.0f};
        if (hasGlossy && hasDiffuse)
        {
            Float glossySamplingRate = m_bsdf->getGlossySamplingRate(bRec);
            averageReflectance = (1.0f-glossySamplingRate)*m_bsdf->getDiffuseReflectance(its)+glossySamplingRate*m_bsdf->getSpecularReflectance(its);
        }
        else if (hasGlossy)
            averageReflectance = m_bsdf->getSpecularReflectance(its);
        else if (hasDiffuse)
            averageReflectance = m_bsdf->getDiffuseReflectance(its);

        //ignoring potential contribution from delta component
        const bool hasDelta = m_bsdf->hasComponent(BSDF::EDelta);
        if (hasDelta)
            averageReflectance *= 1.0f-m_bsdf->getDeltaSamplingRate(bRec);

        const Spectrum irradiance = Spectrum{m_region->statistics.flux*m_guidingData.m_productIntegral};
        const Spectrum reflectedRadiance = averageReflectance * irradiance;

        return reflectedRadiance;
    }

    FINLINE Spectrum estimateIncidentRadiance(const Vector& wiWorld) const
    {
        Guiding_Assert(m_region);
        Guiding_Assert(m_guidingData.valid());
        if (!m_guidingData.valid())
            return Spectrum{0.0f};
        return m_region->statistics.flux*m_guidingData.m_liDistribution.pdf(wiWorld);
    }

    FINLINE Float getRadius() const {
        Guiding_Assert(m_region);
        return m_region->radius;
    }

    FINLINE const AABB& getRegionBounds() const {
        Guiding_Assert(m_region);
        return m_region->bound;
    }

    FINLINE const AABB& getSampleBounds() const {
        Guiding_Assert(m_region);
        return m_region->statsSinceLastSplit.sampleBound;
    }

    FINLINE uint32_t getNumMixtureComponents() const {
        if (!m_guidingData.m_numDistributions)
            return 0;
        return m_guidingData.m_distributions[0].getK();
    }

    //TODO: should implement Mitsuba's RTTI interface
    //MTS_DECLARE_CLASS()

    FINLINE bool bsdfTypeAllowsGuiding() const {
        const int type = m_bsdf->getType();
        //TODO: should also exclude materials with very low roughness (e.g. < 0.01)
        return !(isBSDFPureSpecular(type) || isBSDFTransmissionExcludingNull(type));
        //return !isBSDFPureSpecular(type);
        //return !(isBSDFPureSpecular(type) || isBSDFTransmissionExcludingNull(type));
    }

	/// Returns true if bsdf is purely specular (i.e. only delta transmission and/or reflection)
    FINLINE static constexpr bool isBSDFPureSpecular(const int type) {
        return (type & BSDF::EAll & ~BSDF::EDelta) == 0;
	}

    /// Returns true if the bsdf contains a transmission component which is not a Null pseudo-intersection
    FINLINE static constexpr bool isBSDFTransmissionExcludingNull(const int type) {
        return (type & BSDF::EAll & BSDF::ETransmission & ~BSDF::ENull);
	}

private:

    template<typename, typename> friend class GuidingField;
    template<typename, typename> friend class GuidingFieldFactory;
    template<typename> friend class GuidingDistributionFactory;

    struct GuidingData
    {
        typedef std::integral_constant<size_t, 2> MaxNumProductDistributions;

        /// the region's Li distribution with applied parallax-compensation
        TDistribution m_liDistribution;

        /// product guiding distribution (may be invalid)
        std::array<TDistribution, MaxNumProductDistributions::value> m_distributions;
        /// distribution sampling weights, sum up to 1.0f
        std::array<Float, MaxNumProductDistributions::value> m_weights;
        /// when 0 use the non-product distribution instead
        uint32_t m_numDistributions;
        /// guiding cosine/BSDF product integral (= irradiance/flux, for cosine)
        Float m_productIntegral;

        FINLINE Vector sample(const Point2 sample) const
        {
            Guiding_Assert(m_numDistributions > 0);

            Float weight {0.0f};
            uint32_t i=0;
            for (; i<m_numDistributions-1; ++i)
            {
                const Float nextWeight = weight+m_weights[i];

                if (nextWeight > sample.x)
                    break;

                weight = nextWeight;
            }

            return m_distributions[i].sample(Point2{(sample.x-weight)/m_weights[i], sample.y});
        }

        FINLINE Float pdf(const Vector dir) const
        {
            Guiding_Assert(m_numDistributions > 0);

            Float pdf {0.0f};
            for (uint32_t i=0; i<m_numDistributions; ++i)
                pdf += m_weights[i]*m_distributions[i].pdf(dir);

            return pdf;
        }

        FINLINE bool valid() const {
            return m_numDistributions > 0;
        }

        FINLINE void clear() {
            m_numDistributions = 0;
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "GuidingData [\n";
            for (uint32_t i=0; i<m_numDistributions; ++i)
            {
                oss << '[' << i << "]: " << m_distributions[i].toString() << '\n'
                    << "weight: " << m_weights[i] << '\n';
            }
            oss << "product: " << m_productIntegral << '\n'
                << ']';
            return oss.str();
        }

    } m_guidingData;

    /// bsdf at the current intersection
    const BSDF* m_bsdf;
    int m_glossyComponentId;
    int m_diffuseComponentId;

    /// guiding region at the current intersection
    const TGuidingRegion* m_region;

    /// the probability of sampling the bsdf instead of the
    /// guiding distribution if both are viable options
    Float m_bsdfSamplingProbability;

    bool m_useCosineProduct;
    bool m_useBSDFProduct;
};


GUIDING_NAMESPACE_END
