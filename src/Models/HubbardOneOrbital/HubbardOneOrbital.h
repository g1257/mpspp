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

#ifndef MODEL_HUBBARD_ONE_ORBITAL_H
#define MODEL_HUBBARD_ONE_ORBITAL_H

#include "ModelBase.h"
#include "ParametersModelHubbard.h"
#include "ProgramGlobals.h"

/** Hamiltonian of the Hubbard Model:

  For a chain:
  left corner: [U_0 nup ndown, t_{01} c^\dagger,I]
  right corner: [I,c,U_{n-1} nup ndown]
  on middle site s: [ I                0          0 ]
					[ c                0          0 ]
					[ U_{s} nup ndown t_{s,s+1} c, I]

  For a ladder: ?

*/

namespace Mpspp {

template<typename ParametersSolverType,
         typename InputValidatorType,
         typename SymmetryLocalType,
         typename GeometryType>
class HubbardOneOrbital : public ModelBase<ParametersSolverType,
        InputValidatorType,
        SymmetryLocalType,
        GeometryType> {

	typedef ModelBase<ParametersSolverType,
	InputValidatorType,
	SymmetryLocalType,
	GeometryType> ModelBaseType;

	typedef typename ModelBaseType::MpoLocalType MpoLocalType;
	typedef typename MpoLocalType::MpoFactorType MpoFactorType;
	typedef typename MpoFactorType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairForOperatorType;
	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::ComplexOrRealType ComplexOrRealType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename MpoLocalType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::VectorIntegerType VectorIntegerType;
	typedef ParametersModelHubbard<RealType> ParametersModelType;

	static const int MAX_SITES = ProgramGlobals::MAX_SITES;

public:

	enum {SPIN_UP,SPIN_DOWN};

	HubbardOneOrbital(const ParametersSolverType& solverParams,
	                  InputValidatorType& io,
	                  const GeometryType& geometry)
	    : solverParams_(solverParams),
	      io_(io),
	      geometry_(geometry),
	      hilbert_(4),
	      mp_(io),
	      hamiltonian_(geometry_.numberOfSites())
	{
		// FIXME: CONNECT WITH THE GEOMETRY HERE!!
		ComplexOrRealType tiip1up = -1.0;
		ComplexOrRealType tiip1down = 0.0;
		SizeType n = hamiltonian_.size();
		SizeType wdim = 6;

		SparseMatrixType identity(hilbert_,hilbert_);
		identity.makeDiagonal(hilbert_,1.0);
		identity.toDense();
		SparseMatrixType cup(hilbert_,hilbert_);
		fillDestructionMatrix(cup,SPIN_UP);
		SparseMatrixType cdaggerUp(hilbert_,hilbert_);
		transposeConjugate(cdaggerUp,cup);

		SparseMatrixType cdown(hilbert_,hilbert_);
		fillDestructionMatrix(cdown,SPIN_DOWN);
		SparseMatrixType cdaggerDown(hilbert_,hilbert_);
		transposeConjugate(cdaggerDown,cdown);

		//		SparseMatrixType nupndown = (cdaggerUp*cup) * (cdaggerDown*cdown);

		MpoFactorType mleft(1,wdim);
		RealType mysign = -1.0;
		//		mleft(0,0) = PairForOperatorType(mp_.hubbardU[0]* nupndown,1);
		mleft(0,1) = PairForOperatorType(tiip1up * cdaggerUp,-1);
		mleft(0,2) = PairForOperatorType((mysign * tiip1up) * cup,-1);
		mleft(0,3) = PairForOperatorType(tiip1down * cdaggerDown,-1);
		mleft(0,4) = PairForOperatorType((mysign * tiip1down) * cdown,-1);
		mleft(0,5) = PairForOperatorType(identity,1);
		hamiltonian_(0)=mleft;

		for (SizeType i=1;i<n-1;i++) {
			MpoFactorType m(wdim,wdim);
			m(0,0) = PairForOperatorType(identity,1);
			m(1,0) = PairForOperatorType(cup,-1);
			m(2,0) = PairForOperatorType(cdaggerUp,-1);
			m(3,0) = PairForOperatorType(cdown,-1);
			m(4,0) = PairForOperatorType(cdaggerDown,-1);
			//			m(5,0) = PairForOperatorType(mp_.hubbardU[i]* nupndown,1);
			m(5,1) = PairForOperatorType(tiip1up * cdaggerUp,-1);
			m(5,2) = PairForOperatorType((mysign * tiip1up) * cup,-1);
			m(5,3) = PairForOperatorType(tiip1down * cdaggerDown,-1);
			m(5,4) = PairForOperatorType((mysign * tiip1down) * cdown,-1);
			m(5,5) = PairForOperatorType(identity,1);
			hamiltonian_(i)=m;
		}

		MpoFactorType mright(wdim,1);
		//		mright(5,0) = PairForOperatorType(mp_.hubbardU[n-1]* nupndown,1);
		mright(4,0) = PairForOperatorType(cdaggerDown,-1);
		mright(3,0) = PairForOperatorType(cdown,-1);
		mright(2,0) = PairForOperatorType(cdaggerUp,-1);
		mright(1,0) = PairForOperatorType(cup,-1);
		mright(0,0) = PairForOperatorType(identity,1);
		hamiltonian_(n-1)=mright;
	}

	virtual const MpoLocalType& hamiltonian() const
	{
		return hamiltonian_;
	}

	virtual const ParametersSolverType& solverParams() const { return solverParams_; }

	virtual const GeometryType& geometry() const { return geometry_; }

	virtual void getOneSite(VectorIntegerType& quantumNumbers,SizeType site) const
	{
		quantumNumbers.push_back(0);
		quantumNumbers.push_back(1);
		quantumNumbers.push_back(MAX_SITES);
		quantumNumbers.push_back(1+MAX_SITES);
		//		VectorIntegerType partition_;
		//		PsimagLite::Vector<SizeType>::Type permutation_;
		//		PsimagLite::Vector<SizeType>::Type permutationInverse_;
	}

private:

	void fillDestructionMatrix(SparseMatrixType& cm,SizeType spin) const
	{
		MatrixType m(4,4);
		if (spin==SPIN_UP) {
			m(0,1) = 1;
			m(2,3) = 1;
		} else {
			m(0,2) = 1;
			m(1,3) = -1;
		}
		fullMatrixToCrsMatrix(cm,m);
	}

	const ParametersSolverType& solverParams_;
	InputValidatorType& io_;
	const GeometryType& geometry_;
	SizeType hilbert_;
	ParametersModelType mp_;
	MpoLocalType hamiltonian_;

}; // HubbardOneOrbital

} // namespace Mpspp

/*@}*/
#endif // MODEL_HUBBARD_ONE_ORBITAL_H

