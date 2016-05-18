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
#include "ContractedLocal.h"

namespace Mpspp {

template<typename ParametersSolverType_,
		 typename InputValidatorType_,
		 typename SymmetryLocalType_,
		 typename GeometryType_,
		 typename ConcurrencyType_>
class ModelBase {

public:

	typedef ParametersSolverType_ ParametersSolverType;
	typedef typename ParametersSolverType::ComplexOrRealType ComplexOrRealType;
	typedef InputValidatorType_ InputValidatorType;
	typedef SymmetryLocalType_ SymmetryLocalType;
	typedef GeometryType_ GeometryType;
	typedef ConcurrencyType_ ConcurrencyType;
	typedef MatrixProductOperator<ComplexOrRealType,SymmetryLocalType> MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MpsLocalType MpsLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ProgramGlobals::Vector<RealType>::Type VectorType;
	typedef ContractedLocal<MatrixProductOperatorType> ContractedLocalType;
	typedef ModelHelper<ContractedLocalType> ModelHelperType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef ReflectionSymmetryEmpty<SparseMatrixType> ReflectionSymmetryType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type MatrixType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;

	virtual ~ModelBase() {}

	virtual const MatrixProductOperatorType& hamiltonian() const=0;

	virtual const ParametersSolverType& solverParams() const=0;

	virtual const GeometryType& geometry() const=0;

	virtual void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
	{
		modelHelper.fullHamiltonian(matrix);
	}

	virtual void matrixVectorProduct(VectorType& x,const VectorType& y,const ModelHelperType& modelHelper) const
	{
		modelHelper.matrixVectorProduct(x,y);
	}

	virtual void getOneSite(std::vector<size_t>& quantumNumbers,size_t site) const=0;

}; // ModelBase

} // namespace Mpspp

/*@}*/
#endif // MODEL_BASE_H

