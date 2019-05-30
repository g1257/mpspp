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

#ifndef MPO_LOCAL_H
#define MPO_LOCAL_H

#include "ProgramGlobals.h"
#include "MpsLocal.h"
#include "MpoFactor.h"
#include "SymmetryHelper.h"

namespace Mpspp {

/* PSIDOC MpsLocal
		Essentially same as MpoFactor but with little bit
		more information.
		*/

template<typename ComplexOrRealType_,typename SymmetryLocalType>
class MpoLocal {

public:

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef MpsLocal<ComplexOrRealType,SymmetryLocalType> MpsLocalType;
	typedef SymmetryHelper<SymmetryLocalType> SymmetryHelperType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef MpoFactor<RealType,ComplexOrRealType> MpoFactorType;
	typedef typename PsimagLite::Vector<MpoFactorType>::Type VectorType;

	MpoLocal(size_t nsites) : data_(nsites,MpoFactorType(0,0)) {}

	const MpoFactorType& operator()(size_t site) const
	{
		assert(site<data_.size());
		return data_[site];
	}

	MpoFactorType& operator()(size_t site)
	{
		assert(site<data_.size());
		return data_[site];
	}

	size_t size() const { return data_.size(); }

private:

	VectorType data_;

}; // MpoLocal

} // namespace Mpspp

/*@}*/
#endif // MPO_LOCAL_H

