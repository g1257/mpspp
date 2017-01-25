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

#ifndef HEISENBERG_SPIN_ONE_HALF
#define HEISENBERG_SPIN_ONE_HALF

#include "ModelBase.h"
#include "ParametersHeisenbergSpinOneHalf.h"
#include "ProgramGlobals.h"
#include "MpoLocal.h"

/** Hamiltonian of the Heisenberg Model:

  For a chain:
  left corner:
  right corner:
  on middle site s:

  For a ladder: ?

*/

namespace Mpspp {

template<typename ParametersSolverType,
         typename InputValidatorType,
         typename SymmetryLocalType,
         typename GeometryType>
class HeisenbergSpinOneHalf : public ModelBase<ParametersSolverType,
        InputValidatorType,
        SymmetryLocalType,
        GeometryType> {

    typedef ModelBase<ParametersSolverType,
    InputValidatorType,
    SymmetryLocalType,
    GeometryType> ModelBaseType;

    typedef typename ModelBaseType::MpoLocalType MpoLocalType;
    typedef typename MpoLocalType::MpoFactorType MpoFactorType;
    typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
    typedef typename ModelBaseType::ModelHelperType ModelHelperType;
    typedef typename ModelBaseType::VectorType VectorType;
    typedef typename ModelBaseType::ComplexOrRealType ComplexOrRealType;
    typedef typename ParametersSolverType::RealType RealType;
    typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
    typedef typename MpoLocalType::MpsLocalType MpsLocalType;
    typedef typename MpsLocalType::VectorIntegerType VectorIntegerType;

    typedef ParametersHeisenbergSpinOneHalf<RealType> ParametersModelType;
    typedef typename MpoFactorType::OperatorType OperatorType;

    static const int MAX_SITES = ProgramGlobals::MAX_SITES;

public:

	HeisenbergSpinOneHalf(const ParametersSolverType& solverParams,
	                      InputValidatorType& io,
	                      const GeometryType& geometry)
	    : solverParams_(solverParams),
	      io_(io),
	      geometry_(geometry),
	      hilbert_(2),
	      mp_(io),
	      hamiltonian_(geometry_.numberOfSites())
	{
		// FIXME: CONNECT WITH THE GEOMETRY HERE!!
		RealType J = 1.0;
		RealType Jz = 1.0;
		RealType Jover2 = 0.5*J;
		SizeType n = hamiltonian_.size();
		SizeType wdim = 5;

		SparseMatrixType identity(hilbert_,hilbert_);
		identity.makeDiagonal(hilbert_,1.0);
		SparseMatrixType zero(hilbert_,hilbert_);
		zero.makeDiagonal(hilbert_,0.0);

		OperatorType zeroop(zero,1);

		SparseMatrixType splus(hilbert_,hilbert_);
		fillSplusMatrix(splus);
		SparseMatrixType sminus(hilbert_,hilbert_);
		transposeConjugate(sminus,splus);

		SparseMatrixType sz(hilbert_,hilbert_);
		fillSzMatrix(sz);


		MpoFactorType mleft(1,wdim);
		mleft(0,0) = zero;
		mleft(0,1) = Jover2*sminus;
		mleft(0,2) = Jover2*splus;
		mleft(0,3) = Jz*sz;
		mleft(0,4) = identity;
		hamiltonian_(0)=mleft;

		for (SizeType i=1;i<n-1;i++) {
			MpoFactorType m(wdim,wdim);
			m.setTo(zeroop);
			m(0,0) = identity;
			m(1,0) = splus;
			m(2,0) = sminus;
			m(3,0) = sz;

			m(4,1) = Jover2*sminus;
			m(4,2) = Jover2*splus;
			m(4,3) = Jz*sz;
			m(4,4) = identity;
			hamiltonian_(i)=m;
			if (i > 1)
				assert(hamiltonian_(i) == hamiltonian_(1));
		}

		MpoFactorType mright(wdim,1);
		mright(4,0) = zero;
		mright(3,0) = sz;
		mright(2,0) = sminus;
		mright(1,0) = splus;
		mright(0,0) = identity;
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
        if (solverParams_.options.find("nolocalsymm") != PsimagLite::String::npos) {
			quantumNumbers.resize(2,0);
			return;
		}

		quantumNumbers.push_back(0);
		quantumNumbers.push_back(1);
	}

private:

	void fillSplusMatrix(SparseMatrixType& cm) const
	{
		MatrixType m(hilbert_,hilbert_);
		m(1,0) = 1;
		fullMatrixToCrsMatrix(cm,m);
	}

	void fillSzMatrix(SparseMatrixType& cm) const
	{
		MatrixType m(hilbert_,hilbert_);
		m(0,0) = 0.5;
		m(1,1) = -0.5;
		fullMatrixToCrsMatrix(cm,m);
	}

	const ParametersSolverType& solverParams_;
	InputValidatorType& io_;
	const GeometryType& geometry_;
	SizeType hilbert_;
	ParametersModelType mp_;
	MpoLocalType hamiltonian_;

}; // HeisenbergSpinOneHalf

} // namespace Mpspp

/*@}*/
#endif // HEISENBERG_SPIN_ONE_HALF

