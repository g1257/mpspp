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

#ifndef STATE_PREDICTOR_H
#define STATE_PREDICTOR_H

#include "ProgramGlobals.h"
#include "RandomForTests.h"

namespace Mpspp {

template<typename RealType,typename VectorType>
class StatePredictor {

public:

	StatePredictor()
	    : rng_(3433117),energy_(0),symmetrySector_(0)
	{}

	void createRandomVector(VectorType& y,SizeType offset,SizeType final)
	{
		typedef typename VectorType::value_type ComplexOrRealType;
		ComplexOrRealType tmp;
		typename PsimagLite::Real<ComplexOrRealType>::Type atmp=0;
		y.resize(final-offset);
		for (SizeType i=offset;i<final;i++) {
			myRandomT(tmp);
			y[i]=tmp;
			atmp += PsimagLite::real(y[i]*PsimagLite::conj(y[i]));
		}

		atmp = 1.0 / sqrt (atmp);
		for (SizeType i=offset;i<final;i++) y[i] *= atmp;
		vectorSaved_=y;
	}

	void push(const RealType& e,const VectorType& v,SizeType symmetrySector)
	{
		energy_ = e;
		vectorSaved_ = v;
		symmetrySector_ = symmetrySector;
	}

	const RealType& energy() const { return energy_; }

	const VectorType& vector() const { return vectorSaved_; }

	SizeType symmSector() const { return symmetrySector_; }

private:

	void myRandomT(std::complex<RealType>& value)
	{
		value = std::complex<RealType>(rng_() - 0.5, rng_() - 0.5);
	}

	void myRandomT(RealType& value)
	{
		value = rng_() - 0.5;
	}

	PsimagLite::RandomForTests<RealType> rng_;
	RealType energy_;
	VectorType vectorSaved_;
	SizeType symmetrySector_;
}; // StatePredictor

} // namespace Mpspp

/*@}*/
#endif // STATE_PREDICTOR_H

