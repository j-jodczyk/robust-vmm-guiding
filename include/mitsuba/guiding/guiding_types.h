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

#ifndef _GUIDING_CONFIG_H_
#define _GUIDING_CONFIG_H_

#include <mitsuba/mitsuba.h>

#include <pmm/pmm.h>
#include <pmm/ParametricMixtureModel.h>
#include <pmm/VMFKernel.h>
#include <pmm/DirectionalData.h>

#include <mitsuba/guiding/guiding.h>
#include <mitsuba/guiding/GuidingField.h>
#include <mitsuba/guiding/GuidingFieldFactory.h>
#include <mitsuba/guiding/pathguidingmixturestats.h>
#include <mitsuba/guiding/pathguidingmixturestatsfactory.h>
#include <mitsuba/guiding/BSDFOracle.h>

GUIDING_NAMESPACE_BEGIN

// Sample data and mixture configuration.
// This can be changed to use Scalar8 for AVX, instead of Scalar4 for SSE2
// But you also need to set the -mavx compiler flag. Otherwise, it is still using SSE2 internally.
typedef lightpmm::DirectionalData PathGuidingSampleData;
typedef lightpmm::ParametricMixtureModel<lightpmm::VMFKernel<lightpmm::Scalar4>, (32+lightpmm::Scalar4::Width::value-1)/lightpmm::Scalar4::Width::value> VMFMixture;
typedef GuidingField<VMFMixture, PathGuidingMixtureStats<VMFMixture>> GuidingFieldType;
typedef GuidingFieldType::GuidingTreeType GuidingTreeType;
typedef lightpmm::ParametricMixtureModel<lightpmm::VMFKernel<lightpmm::Scalar4>, 1> BSDFVMFMixture;
typedef BSDFOracleVMF<BSDFVMFMixture> BSDFOracle;

GUIDING_NAMESPACE_END

#endif
