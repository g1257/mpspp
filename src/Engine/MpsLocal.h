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

#ifndef MPS_LOCAL_H
#define MPS_LOCAL_H

#include "ProgramGlobals.h"
#include "MpsFactor.h"

namespace Mpspp {

template<typename ComplexOrRealType_,typename SymmetryLocalType_>
class MpsLocal {

	// FIXME: IDEA: PULL SYMMETRY OUT, PASS THROUGH FUNCTIONS

	/* PSIDOC MpsLocal
			Contains the Vector of MpsFactor, which are A_ and B_
			operators obtained from the SVD of the ground-state of the
			full system. Here, the index of the vector corresponds to
			the "site". It also intializes the size of A_ and B_ -
			containing the "triangular" size structure of each, and later
			randomly set elements of the operator.
			*/

public:

	typedef SymmetryLocalType_ SymmetryLocalType;
	typedef typename SymmetryLocalType::IoInputType IoInputType;
	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef MpsFactor<ComplexOrRealType,SymmetryFactorType> MpsFactorType;
	typedef typename MpsFactorType::VectorType VectorType;
	typedef typename MpsFactorType::SparseMatrixType SparseMatrixType;
	typedef typename MpsFactorType::MatrixType MatrixType;
	typedef typename MpsFactorType::VectorRealType VectorRealType;
	typedef typename MpsFactorType::VectorIntegerType VectorIntegerType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	MpsLocal(SizeType nsites)
	    : nsites_(nsites),center_(0)
	{}

	~MpsLocal()
	{
		for (SizeType i=0;i<A_.size();i++)
			if (A_[i]) delete A_[i];
		for (SizeType i=0;i<B_.size();i++)
			if (B_[i]) delete B_[i];
	}

	void initialGuess(SizeType currentSite,const SymmetryLocalType&,SizeType)
		{
			center_ = currentSite;
			SizeType d = 2;
			SizeType middle = static_cast<SizeType>(nsites_/2);
			SizeType x = nsites_ - currentSite;
			SizeType n = (currentSite<middle) ? pow(d,currentSite+1) : pow(d,x);

			MpsFactorType* mpsFactor = new MpsFactorType(MpsFactorType::TYPE_B);
			std::cout << "currentSite=" << currentSite << ", n=" << n << "\n";
			mpsFactor->setRandom(currentSite,n);
			B_.push_back(mpsFactor);
		}

	//! Returns the number of sites
	SizeType sites() const
	{
		return nsites_;
	}

	//! tmpVec[i] --> M^\sigma2 _ {a1,a2}
	template<typename SomeTruncationType>
	void move(SomeTruncationType& truncation,
	          SizeType currentSite,
	          const VectorType& v,
	          SizeType direction,
	          SizeType symmetrySector,
	          const SymmetryFactorType& symm)
	{
		center_ = currentSite;
		if (direction==ProgramGlobals::TO_THE_RIGHT) {
			if (currentSite==A_.size()) {
				MpsFactorType* mpsFactor = new MpsFactorType(MpsFactorType::TYPE_A);
				SizeType n = symm.left().size();
				mpsFactor->setRandom(currentSite,n);
				A_.push_back(mpsFactor);
			}
			assert(currentSite<A_.size());
            A_[currentSite]->move(truncation,v,symmetrySector,symm);
			std::cout<<"moved A["<<currentSite<<"].row= ";
			std::cout<<A_[currentSite]->operator()().row()<<"\n";
		} else {
			SizeType nsites = symm.super().block().size();
			SizeType siteToSet = nsites - 1 - currentSite;
			assert(siteToSet<nsites);

			if (siteToSet==B_.size()) {
				MpsFactorType* mpsFactor = new MpsFactorType(MpsFactorType::TYPE_B);
				SizeType n = symm.right().size();
				mpsFactor->setRandom(currentSite,n);
				B_.push_back(mpsFactor);
			}
            B_[siteToSet]->move(truncation,v,symmetrySector,symm);
			std::cout<<"moved B["<<(siteToSet)<<"].row= ";
			std::cout<<B_[siteToSet]->operator()().row()<<"\n";
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
			A_[site]->truncate(cutoff,trunc);
		} else {
			SizeType siteToSet = nsites-site-1;
			B_[siteToSet]->truncate(cutoff,trunc);
		}
	}

	const MpsFactorType& A(SizeType site) const
	{
		assert(site<A_.size());
		return *(A_[site]);
	}

	const MpsFactorType& B(SizeType site) const
	{
		assert(site<B_.size());
		return *(B_[site]);
	}

private:

	// copy ctor:
	MpsLocal(const MpsLocal& other);

	// assignment
	MpsLocal& operator=(const MpsLocal& other);

	//	const SymmetryLocalType& symm_;
	SizeType nsites_;
	SizeType center_;
	typename PsimagLite::Vector<MpsFactorType*>::Type B_;
	typename PsimagLite::Vector<MpsFactorType*>::Type A_;
}; // MpsLocal

} // namespace Mpspp

/*@}*/
#endif // MPS_LOCAL_H


/*

//	void computeIntermediate(MatrixType& matrixNew,
//	                         const MatrixType& matrixOld,
//	                         SizeType center,
//	                         const SymmetryLocalType& symmLocal) const
//	{
//		const SparseMatrixType& Bmatrix = B_[center]->operator()();
//		SparseMatrixType Btranspose;
//		transposeConjugate(Btranspose,Bmatrix);

//		matrixNew.resize(Bmatrix.row(),Bmatrix.row());
//		matrixNew.setTo(0.0);

//		if (matrixOld.n_row()==0) {
//			assert(center+1==B_.size());
//			for (SizeType x=0;x<Bmatrix.row();x++) {
//				for (int k=Bmatrix.getRowPtr(x);k<Bmatrix.getRowPtr(x+1);k++) {
//					SizeType sigma = Bmatrix.getCol(k);
//					for (int k2=Btranspose.getRowPtr(sigma);
//					     k2<Btranspose.getRowPtr(sigma+1);
//					     k2++) {
//						SizeType y = Btranspose.getCol(k2);
//						matrixNew(x,y) += Bmatrix.getValue(k)*
//						        std::conj(Btranspose.getValue(k2));
//					}
//				}
//			}

//			return;
//		}

//		const SymmetryFactorType& symm = symmLocal(center);
//		for (SizeType anm2=0;anm2<Bmatrix.row();anm2++) {
//			for (int k=Bmatrix.getRowPtr(anm2);k<Bmatrix.getRowPtr(anm2+1);k++) {
//				PairType p = symm.right().unpack(Bmatrix.getCol(k));
//				SizeType sigmanm2=p.first;
//				SizeType anm1=p.second;
//				for (SizeType apnm1=0;apnm1<matrixOld.n_row();apnm1++) {
//					SizeType x = symm.right().pack(sigmanm2,apnm1);
//					for (int k2=Btranspose.getRowPtr(x);k2<Btranspose.getRowPtr(x+1);k2++) {
//						SizeType apnm2 = Btranspose.getCol(k2);
//						matrixNew(anm2,apnm2) += matrixOld(anm1,apnm1)*
//						        Bmatrix.getValue(k)*Btranspose.getValue(k2);
//					} // k2
//				} // apnm1
//			} // k
//		} // anm2
//	}

//	RealType normA(const SymmetryLocalType& symm) const
//	{
//		MatrixType tmpOld;

//		for (SizeType i=0;i<A_.size();i++) {
//			SizeType center = i;
//			MatrixType tmpNew;
//			computeIntermediateA(tmpNew,tmpOld,center,symm);
//			tmpOld = tmpNew;
//		}

//		ComplexOrRealType sum = 0;
//		for (SizeType i=0;i<tmpOld.n_row();i++) {
//			for (SizeType j=0;j<tmpOld.n_col();j++) {
//				sum += tmpOld(i,j);
//			}
//		}
//		return std::real(sum);
//	}

//	void computeIntermediateA(MatrixType& matrixNew,
//	                          const MatrixType& matrixOld,
//	                          SizeType center,
//	                          const SymmetryLocalType& symmLocal) const
//	{
//		const SymmetryFactorType& symm = symmLocal(center+1);
//		SizeType leftSize = symm.left().split();
//		SizeType leftSize2 = symm.left().size();

//		SparseMatrixType Atranspose;
//		transposeConjugate(Atranspose,A_[center]->operator()());
//		matrixNew.resize(leftSize*leftSize,leftSize2*leftSize2);
//		matrixNew.setTo(0.0);
//		assert(Atranspose.row()==symm.left().size());
//		for (SizeType a1=0;a1<Atranspose.row();a1++) {
//			for (int k1=Atranspose.getRowPtr(a1);k1<Atranspose.getRowPtr(a1+1);k1++) {
//				PairType a0sigma1=symm.left().unpack(Atranspose.getCol(k1));
//				SizeType a0 = a0sigma1.first;
//				SizeType sigma1= a0sigma1.second;
//				for (SizeType a1p=0;a1p<Atranspose.row();a1p++) {
//					for (int k2=Atranspose.getRowPtr(a1p);
//					     k2<Atranspose.getRowPtr(a1p+1);
//					     k2++) {

//						PairType a0sigma1p=symm.left().unpack(Atranspose.getCol(k2));
//						SizeType a0p = a0sigma1p.first;
//						SizeType sigma1p= a0sigma1p.second;
//						if (sigma1!=sigma1p) continue;
//						matrixNew(a0+a0p*leftSize,a1+a1p*leftSize2) +=
//						        Atranspose.getValue(k1)*std::conj(Atranspose.getValue(k2));
//					} // k2
//				} // a1p
//			} // k1
//		} // a1

//		if (matrixOld.n_row()==0) return;
//		MatrixType m = matrixOld * matrixNew;
//		matrixNew = m;
//	}

*/

