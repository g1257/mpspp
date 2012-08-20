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

#ifndef MODEL_BASE_H
#define MODEL_BASE_H
#include "MatrixProductOperator.h"
#include "ModelHelper.h"
#include "ReflectionSymmetryEmpty.h"

namespace Mpspp {

template<typename ParametersSolverType_,
		 typename InputValidatorType_,
		 typename GeometryType_,
		 typename ConcurrencyType_>
class ModelBase {

public:

	typedef ParametersSolverType_ ParametersSolverType;
	typedef InputValidatorType_ InputValidatorType;
	typedef GeometryType_ GeometryType;
	typedef ConcurrencyType_ ConcurrencyType;
	typedef MatrixProductOperator MatrixProductOperatorType;
	typedef typename ParametersSolverType_::RealType RealType;
	typedef ModelHelper<RealType,RealType> ModelHelperType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef ReflectionSymmetryEmpty<SparseMatrixType> ReflectionSymmetryType;
	typedef typename ProgramGlobals::Vector<RealType>::Type VectorType;

	virtual const MatrixProductOperatorType& hamiltonian(size_t site) const=0;

	virtual const ParametersSolverType& solverParams() const=0;

	virtual void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
	{

//		data_.resize(total,total);
//		for (size_t i=0;i<total;i++) {
//			data_.setRow(i,counter);
//			PairType iLa2 = symm.super().unpack(i+offset);
//			PairType a1sigma2 = symm.left().upack(iLa2.first);
//			size_t sigma2 = a1sigma2.second;
//			size_t a1 = a1sigma2.first;
//			for (size_t b1=0;b1<hamiltonian.rank();b1++) {
//				const SparseMatrixType& cLm = cL(b1);
//				for (size_t k1=cLm.getRowPtr(a1);k1<cLm.getRowPtr(a1+1);k1++) {
//					size_t a1prime = cLm.getCol(k1);
//					for (size_t b2=0;b2<hamiltonian.rank();b2++) {
//						const MatrixType& wm = hamiltonian(b1,b2);
//						const SparseMatrixType& cLm = cL(b2);
//						for (size_t k2=cRm.getRowPtr(a2);k2<cRm.getRowPtr(a2+1);k2++) {
//							size_t a2prime = cLm.getCol(k2);
//							for (size_t sigma2prime = 0;sigma2prime<hilbertSize;sigma2prime++) {
//								size_t iLprime = packLeft.pack(a1prime,sigma2prime);
//								size_t iprime = packSuper.pack(iLprime,a2prime);
//								tmpVector[iprime] += cLm.getValue(k1) * wm(sigma2,sigma2prime) * cRm.getValue(k2);
//							}
//						}
//					}
//				}
//			}
//			counter += addThisRow(tmpVector);
//		}
//		data_.setRow(total,counter);
	}

	virtual void matrixVectorProduct(VectorType& x,const VectorType& y,const ModelHelperType& modelHelper) const=0;
}; // ModelBase

} // namespace Mpspp

/*@}*/
#endif // MODEL_BASE_H

