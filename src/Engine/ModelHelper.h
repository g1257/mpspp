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

#ifndef MODEL_HELPER_H
#define MODEL_HELPER_H
#include "ProgramGlobals.h"

namespace Mpspp {

template<typename LeftRightSuperType>
class ModelHelper {

public:

	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename LeftRightSuperType::ComplexOrRealType ComplexOrRealType;
	typedef typename LeftRightSuperType::ContractedPartType ContractedPartType;
	typedef typename ContractedPartType::ContractedFactorType ContractedFactorType;
	typedef typename LeftRightSuperType::MatrixProductStateType MatrixProductStateType;
	typedef typename MatrixProductStateType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename ContractedPartType::SparseMatrixType SparseMatrixType;
	typedef typename ProgramGlobals::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename LeftRightSuperType::MatrixProductOperatorType MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;

	ModelHelper(const LeftRightSuperType& lrs,
				size_t symmetrySector,
				size_t currentSite,
				size_t direction,
				const MpoFactorType& hamiltonian)
	: lrs_(lrs),
	  symmetrySector_(symmetrySector),
	  currentSite_(currentSite),
	  direction_(direction),
	  hamiltonian_(hamiltonian),
	  symmetry_(lrs_.A().symmetry()(currentSite_))
	{}

	size_t size() const
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need size(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	size_t symmetrySector() const { return symmetrySector_; }

	size_t hilbertSize() const { return hamiltonian_(0,0).row(); }

	const LeftRightSuperType& lrs() const { return lrs_; }

	const MpoFactorType& hamiltonian() const { return hamiltonian_; }

	const SymmetryFactorType& symmetry() const { return symmetry_; }

	const ContractedFactorType& contractedFactorLeft() const { return lrs_.contractedLeft()(currentSite_); }

	const ContractedFactorType& contractedFactorRight() const { return lrs_.contractedRight()(currentSite_); }

private:

	const LeftRightSuperType& lrs_;
	size_t symmetrySector_;
	size_t currentSite_;
	size_t direction_;
	size_t hilbertSize_;
	const MpoFactorType& hamiltonian_;
	const SymmetryFactorType& symmetry_;

}; // ModelHelper

} // namespace Mpspp

/*@}*/
#endif // MODEL_HELPER_H

