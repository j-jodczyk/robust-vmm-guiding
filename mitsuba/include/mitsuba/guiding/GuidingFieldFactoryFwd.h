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

#ifndef GUIDINGFIELDFACTORYFWD_H
#define GUIDINGFIELDFACTORYFWD_H

#include "guiding.h"

GUIDING_NAMESPACE_BEGIN

template<typename TDistribution>
class GuidingDistributionFactory;

template<typename TDistribution, typename TStatistics>
class GuidingStatsFactory;

template<typename TDistribution, typename TStatistics>
class GuidingFieldFactory;

GUIDING_NAMESPACE_END

#endif // GUIDINGFIELDFACTORYFWD_H
