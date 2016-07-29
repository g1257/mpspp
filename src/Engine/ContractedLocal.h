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

#ifndef CONTRACTED_LOCAL_H
#define CONTRACTED_LOCAL_H

#include "ProgramGlobals.h"
#include "ContractedFactor.h"

namespace Mpspp {

/* PSIDOC ContractedLocal
		ContractedLocal contains the two ContractedFactors L_ and R_, which
		can be seen as the wavefunction of left section of the lattice and
		right section of the lattice. They are calculated using a recursive
		relation found in equation 192 and 193 of Schollwock-2011. This
		templated class is used to direct the algorithm to ContractedFactor
		where the "guts" of move left/right of the finite/infite sweeps
		exists.
		*/

template<typename MatrixProductOperatorType_>
class ContractedLocal {

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT,
		  TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT,
		  PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:

	typedef MatrixProductOperatorType_ MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MatrixProductOperatorType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::MpsFactorType MpsFactorType;
	typedef typename MpsFactorType::VectorIntegerType VectorIntegerType;
	typedef typename MpsLocalType::SymmetryLocalType SymmetryLocalType;
	typedef typename MpsLocalType::ComplexOrRealType ComplexOrRealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef ContractedFactor<MatrixProductOperatorType> ContractedFactorType;
	typedef typename MatrixProductOperatorType::SymmetryHelperType SymmetryHelperType;

	ContractedLocal(const MpsLocalType& abState,const MatrixProductOperatorType& h)
		: abState_(abState),
		  h_(h),
		  R_(abState.sites(),ProgramGlobals::PART_RIGHT),
		  L_(abState.sites(),ProgramGlobals::PART_LEFT)
	{}

	void initialGuess(SizeType currentSite,
					  const SymmetryHelperType& symmHelper,
					  SizeType nsites)
	{
		ContractedFactorType* ptr = (currentSite == 0) ? 0 : &R_[currentSite];
		R_[currentSite+1].build(abState_.B(nsites-1-currentSite),
								h_(nsites-1-currentSite),
								ptr,
								symmHelper,
								currentSite);
	}

	//! From As (or Bs) and Ws reconstruct *this
	void move(SizeType currentSite,
			  SizeType direction,
			  const SymmetryHelperType& symmHelper)
	{
		if (direction==TO_THE_RIGHT) {
			moveLeft(currentSite,abState_,symmHelper);
		} else {
			moveRight(currentSite,abState_,symmHelper);
		}
	}

	template<typename SomeTruncationType>
	void truncate(SizeType site,
				  SizeType part,
				  SizeType cutoff,
				  SizeType nsites,
				  const SomeTruncationType& trunc)
	{
		if (part==ProgramGlobals::PART_LEFT) {
			if (site+1>=L_.size()) return;
			L_[site+1].truncate(cutoff,trunc);
		} else {
			SizeType siteToSet = nsites - site;
			R_[siteToSet].truncate(cutoff,trunc);
		}
	}

	const ContractedFactorType& operator()(SizeType currentSite,
										   SizeType leftOrRight) const
	{
		if (leftOrRight == PART_LEFT) {
			assert(currentSite<L_.size());
		} else {
			assert(currentSite<R_.size());
		}

		return (leftOrRight == PART_LEFT) ? L_[currentSite] : R_[currentSite];
	}

private:

	void moveLeft(SizeType currentSite,
				  const MpsLocalType& abState,
				  const SymmetryHelperType& symm)
	{
        if (currentSite+1==L_.size()) return;
		assert(currentSite+1<L_.size());

		L_[currentSite+1].move(abState.A(currentSite),
							   h_(currentSite),
							   L_[currentSite],
							   symm,
							   currentSite);
		std::cout<<"set L_["<<(currentSite+1)<<"]="<<L_[currentSite+1].row()<<"\n";
	}

	void moveRight(SizeType currentSite,
				   const MpsLocalType& abState,
				   const SymmetryHelperType& symm)
	{
		SizeType nsites = symm.symmLocal()(currentSite).super().block().size();
		SizeType siteToSet = nsites - currentSite;
		if (siteToSet==R_.size()) return;
		assert(siteToSet<R_.size());
		R_[siteToSet].move(abState.B(siteToSet-1),
						   h_(currentSite),
						   R_[siteToSet-1],
				symm,
				currentSite);
		std::cout<<"set R_["<<siteToSet<<"]="<<R_[siteToSet].row()<<"\n";
	}

	const MpsLocalType& abState_;
	const MatrixProductOperatorType& h_;
	typename PsimagLite::Vector<ContractedFactorType>::Type R_;
	typename PsimagLite::Vector<ContractedFactorType>::Type L_;

}; // ContractedLocal

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_LOCAL_H
