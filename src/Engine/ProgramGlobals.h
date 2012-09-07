/*
Copyright (c) 2012, UT-Battelle, LLC
All rights reserved

[MPS++, Version 0.1]
[by K. Al-Hassanieh, Oak Ridge National Laboratory]
[by J. Rincon, Oak Ridge National Laboratory]
[by G.A., Oak Ridge National Laboratory]

See full open source LICENSE under file LICENSE
in the root directory of this distribution.

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/

/** \ingroup MPSPP */
/*@{*/

/*! \file ProgramGlobals.h
 *
 *
 *
 */
#ifndef PROGRAM_LIMITS_H
#define PROGRAM_LIMITS_H
#include "Matrix.h"
#include <vector>
#include "CrsMatrix.h"

namespace Mpspp {
struct ProgramGlobals {
	//		static size_t const MaxNumberOfSites = 300; // max number of sites that a model can use
	static size_t const MaxLanczosSteps = 1000000; // max number of internal Lanczos steps
	static size_t const LanczosSteps = 200; // max number of external Lanczos steps
	static double const LanczosTolerance; // tolerance of the Lanczos Algorithm
	//		enum {FERMION,BOSON};

	enum {TO_THE_RIGHT,TO_THE_LEFT};

	template<typename ComplexOrRealType>
	class Real {
	public:
		typedef ComplexOrRealType Type;
	};

	template<typename RealType>
	class Real<std::complex<RealType> > {
	public:
		typedef RealType Type;
	};

	template<typename SomeFieldType>
	class Matrix {
	public:
#ifndef USE_MATH_UTILS
		typedef PsimagLite::Matrix<SomeFieldType> Type;
#else
		typedef Utils::Matrix<SomeFieldType> Type;
#endif
	};

	template<typename SomeFieldType>
	class Vector {
	public:
#ifndef USE_MATH_UTILS
		typedef std::vector<SomeFieldType> Type;
#else
		typedef Utils::Vector<SomeFieldType> Type;
#endif
	};

	template<typename SomeFieldType>
	class CrsMatrix {
	public:
#ifndef USE_MATH_UTILS
		typedef PsimagLite::CrsMatrix<SomeFieldType> Type;
#else
		typedef Utils::CrsMatrix<SomeFieldType> Type;
#endif
	};
}; // ProgramGlobals

	double const ProgramGlobals::LanczosTolerance = 1e-12;
}; // namespace Mpspp
/*@}*/
#endif

