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

/*! \file InternalProductStored.h
 *
 *  A class to encapsulate the product x+=Hy, where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef InternalProductStored_HEADER_H
#define InternalProductStored_HEADER_H

#include <vector>
#include "ProgressIndicator.h"

/* PSIDOC FermionSign
		Used to output the Fermionic Sign which comes from the
		anticommutation relations of fermions. It is generally needed
		for the non-interacting part of the Hamiltonian.
		*/

namespace Mpspp {
template<typename T,typename ModelType>
class InternalProductStored {
public:
	typedef T HamiltonianElementType;
	typedef HamiltonianElementType value_type;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::RealType RealType;
	typedef PsimagLite::Matrix<RealType> MatrixType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	//typedef typename SparseMatrixType::value_type SparseElementType;

	InternalProductStored(ModelType const *model,
	                      ModelHelperType const *modelHelper,
	                      const ReflectionSymmetryType* rs=0)
	    : matrixStored_(2),pointer_(0),progress_("InternalProductStored")
	{
		model_ = model;
		modelHelper_=modelHelper;

		if (!rs) {
			matrixStored_[0].clear();
			model->fullHamiltonian(matrixStored_[0],*modelHelper);
			if (model_->solverParams().options.find("debugmatrix") !=
			        PsimagLite::String::npos) {
				MatrixType fullm;
				crsMatrixToFullMatrix(fullm,matrixStored_[0]);
				if (PsimagLite::isZero(fullm)) std::cerr<<"Matrix is zero\n";
				if (fullm.n_row()>40) {
					printNonZero(fullm,std::cerr);
				} else {
					std::cout<<fullm;
					//					printFullMatrix(fullm,"matrix",1);
				}

			}

			assert(isHermitian(matrixStored_[0],true));
			PsimagLite::OstringStream msg;
			msg<<"fullHamiltonian has rank="<<matrixStored_[0].row();
			msg<<" nonzeros="<<matrixStored_[0].nonZero();
			progress_.printline(msg,std::cout);
			return;
		}
		SparseMatrixType matrix2;
		model->fullHamiltonian(matrix2,*modelHelper);
		rs->transform(matrixStored_[0],matrixStored_[1],matrix2);
		PsimagLite::OstringStream msg;
		msg<<" sector="<<matrixStored_[0].row()<<" and sector="<<matrixStored_[1].row();
		progress_.printline(msg,std::cout);
	}

	SizeType rank() const { return matrixStored_[pointer_].row(); }

	template<typename SomeVectorType>
	void matrixVectorProduct(SomeVectorType &x, SomeVectorType const &y) const
	{
		matrixStored_[pointer_].matrixVectorProduct(x,y);
	}

//	HamiltonianElementType operator()(SizeType i,SizeType j) const
//	{
//		return matrixStored_[pointer_](i,j);
//	}

private:
	ModelType const *model_;
	ModelHelperType const *modelHelper_;
	typename PsimagLite::Vector<SparseMatrixType>::Type matrixStored_;
	SizeType pointer_;
	PsimagLite::ProgressIndicator progress_;
}; // class InternalProductStored
} // namespace Mpspp

/*@}*/
#endif
