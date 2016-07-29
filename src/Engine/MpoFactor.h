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

#ifndef MPO_FACTOR_H
#define MPO_FACTOR_H

#include "ProgramGlobals.h"
#include "Operator.h"

namespace Mpspp {

/* PSIDOC MpoFactor: MPO (Matrix Product Operator)
		MpoFactor is the templated class which represents the
		matricies "W" needed to build the hamiltonian. An
		Example of this can be found on eq. 184-188 of
		Schollowock-2011. This matricies "W" are similar
		to the added central site in the DMRG language and are
		needed for the definition of the full Hamiltonian.
		*/

template<typename RealType,typename ComplexOrRealType>
class MpoFactor {

	typedef MpoFactor<RealType,ComplexOrRealType> ThisType;
	typedef Operator<ComplexOrRealType> OperatorType_;
	typedef typename OperatorType_::SparseMatrixType SparseMatrixType;
	typedef PsimagLite::Matrix<OperatorType_> MatrixType;

	static const int MAX_SITES = ProgramGlobals::MAX_SITES;

public:

	typedef OperatorType_ OperatorType;

	MpoFactor(size_t wdim1,size_t wdim2)
		: data_(wdim1,wdim2) {}

	const OperatorType& operator()(size_t i,size_t j) const
	{
		assert(i<n_row() && j<n_col());
		return data_(i,j);
	}

	void setTo(const OperatorType& op)
	{
		data_.setTo(op);
	}

	OperatorType& operator()(size_t i,size_t j)
	{
		assert(i<n_row() && j<n_col());
		return data_(i,j);
	}

	size_t n_row() const { return data_.n_row(); }

	size_t n_col() const { return data_.n_col(); }

	bool operator==(const ThisType& other) const
	{
		return (data_ == other.data_);
	}

private:

	PsimagLite::Matrix<OperatorType> data_;

}; // MpoFactor

} // namespace Mpspp

/*@}*/
#endif // MPO_FACTOR_H

