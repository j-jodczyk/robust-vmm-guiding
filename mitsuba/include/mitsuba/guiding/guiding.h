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

#ifndef _GUIDING_H_
#define _GUIDING_H_

#define GUIDING_NO_ASSERTS
#define GUIDING_NO_VALIDITY_CHECKS

// enable this for more verbose logging with additional statistics
//#define GUIDING_DETAILED_STATISTICS
//#define GUIDING_VALIDATE_INTERMEDIATE_MIXTURE_STATES

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/logger.h>

#define GUIDING_USE_PARALLEL
#define GUIDING_USE_BSDF_PRODUCT


#define GUIDING_NAMESPACE_BEGIN MTS_NAMESPACE_BEGIN namespace Guiding{

#define GUIDING_NAMESPACE_END } MTS_NAMESPACE_END

#ifndef _MSC_VER
#define __forceinline inline __attribute__((always_inline))
#endif

#define GUIDING_INLINE __forceinline
//#define GUIDING_INLINE

#ifndef GUIDING_NO_ASSERTS
#define Guiding_Assert(cond) SAssert(cond);
#define Guiding_Assert_MSG(cond, msg) SAssertEx(cond, msg);
#else
#define Guiding_Assert(cond)
#define Guiding_Assert_MSG(cond, msg)
#endif

GUIDING_NAMESPACE_BEGIN

GUIDING_INLINE bool isValid(const Float &v){
    return (std::isfinite(v) && !std::isnan(v));
}

GUIDING_INLINE bool isFinite(const Vector &vec){
    return (std::isfinite(vec[0]) && std::isfinite(vec[1]) && std::isfinite(vec[2]));
}

GUIDING_INLINE bool isNan(const Vector &vec){
    return (std::isnan(vec[0]) || std::isnan(vec[1]) || std::isnan(vec[2]));
}

GUIDING_INLINE bool isValid(const Vector &vec){
    return (isFinite(vec) && !isNan(vec));
}

enum GuidingFlags
{
    ESplatted = 1<<0 // point does not represent any real scene intersection point
};


FINLINE Float PATHGUIDING_SPECTRUM_TO_FLOAT(Spectrum spectrum)
{
    return spectrum.average();
    //return spectrum.max();
}

GUIDING_NAMESPACE_END


#endif
