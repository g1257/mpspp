/*
Copyright (c) 2012-2013, UT-Battelle, LLC
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

#ifndef OPERATOR_H
#define OPERATOR_H

#include "ProgramGlobals.h"
#include "MpsLocal.h"
#include "MpoFactor.h"

namespace Mpspp {

template<typename ComplexOrRealType>
class Operator {

	typedef Operator<ComplexOrRealType> ThisType;

public:

	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef std::pair<SparseMatrixType,int> PairType;

	explicit Operator(const SparseMatrixType& m, int f)
		:data_(m),fermionSign_(f)
	{}

	explicit Operator(int  = 0)
	{}

	ThisType& operator=(const PairType& pair)
	{
		data_=pair.first;
		fermionSign_=pair.second;
		return *this;
	}

	ThisType& operator=(const SparseMatrixType& m)
	{
		data_=m;
		fermionSign_=1;
		return *this;
	}

	const SparseMatrixType& matrix() const
	{
		return data_;
	}

	const int fermionSign() const
	{
		return fermionSign_;
	}

	bool operator==(const ThisType& op) const
	{
		return (data_ == op.data_ && fermionSign_ == op.fermionSign_);
	}

private:

	SparseMatrixType data_;
	int fermionSign_;

}; // Operator

} // namespace Mpspp

/*@}*/
#endif // OPERATOR_H

