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

template<typename ContractedPartType>
class ModelHelper {

	enum {PART_RIGHT = ProgramGlobals::PART_RIGHT,
		  PART_LEFT = ProgramGlobals::PART_LEFT};

public:

	typedef typename ContractedPartType::ContractedFactorType ContractedFactorType;
	typedef typename ContractedPartType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename MpsLocalType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef typename ContractedPartType::SparseMatrixType SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ContractedPartType::MatrixProductOperatorType MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MpoFactorType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairForOperatorType;
	typedef typename MatrixProductOperatorType::SymmetryHelperType SymmetryHelperType;

	ModelHelper(const ContractedPartType& contractedPart,
				SizeType symmetrySector,
				SizeType currentSite,
				SizeType direction,
				const MpoFactorType& hamiltonian,
				const SymmetryHelperType& symmetry,
				SizeType siteForSymm)
		: contractedPart_(contractedPart),
		  symmetrySector_(symmetrySector),
		  currentSite_(currentSite),
		  part_(direction),
		  hamiltonian_(hamiltonian),
		  symmetry_(symmetry),
		  siteForSymm_(siteForSymm)
	{}

	SizeType size() const
	{
		return symmetry_.symmLocal()(siteForSymm_).super().partitionSize(symmetrySector_);
	}

	//! Eq. (201) but very modified
	void matrixVectorProduct(VectorType& x,const VectorType& y) const
	{
		SizeType nsites = symmetry_.symmLocal().size();
		const SymmetryFactorType& symm = symmetry_.symmLocal()(siteForSymm_);
		SizeType offset = symm.super().partitionOffset(symmetrySector_);
		SizeType total = symm.super().partitionSize(symmetrySector_);

		SizeType leftIndex = (part_ == PART_RIGHT) ? currentSite_ : currentSite_+1;
		SizeType rightIndex = (part_ == PART_RIGHT) ? currentSite_ : nsites - currentSite_ -1;

		const ContractedFactorType& cL = contractedPart_(leftIndex,
														 ProgramGlobals::PART_LEFT);
		const ContractedFactorType& cR = contractedPart_(rightIndex,
														 ProgramGlobals::PART_RIGHT);
		for (SizeType blm1=0;blm1<cL.size();blm1++) {
			const SparseMatrixType& l1 = cL(blm1);
			for (SizeType bl=0;bl<cR.size();bl++) {
				const OperatorType& wOp =  hamiltonian_(blm1,bl);
				const SparseMatrixType& w = wOp.matrix();
				if (w.row()==0) continue;
				//				SparseMatrixType w;
				//				transposeConjugate(w,w1);
				const SparseMatrixType& r1 = cR(bl);
				for (SizeType i=0;i<total;i++) {
					PairType ab = symm.super().unpack(i+offset);
					SizeType alm1=0;
					SizeType sigmaL=0;
					SizeType alB=0;
					SizeType electronsLeft = symmetry_.electronsFromQn(symm.left().
																	   qn(ab.first));
					if (part_==PART_RIGHT) {
						PairType tmpPair1 = symm.left().unpack(ab.first);
						alm1=tmpPair1.first;
						sigmaL=tmpPair1.second;
						alB = ab.second;
						SizeType electronsBlock = symmetry_.electronsFromState(sigmaL);
						assert(electronsBlock<=electronsLeft);
						electronsLeft -= electronsBlock;
					} else {
						PairType tmpPair1 = symm.right().unpack(ab.second);
						alB=tmpPair1.second;
						sigmaL=tmpPair1.first;
						alm1=ab.first;
					}

					RealType fermionSign = (electronsLeft & 1) ? wOp.fermionSign() : 1.0;

					for (int k1=l1.getRowPtr(alm1);k1<l1.getRowPtr(alm1+1);k1++) {
						SizeType alm1p=l1.getCol(k1);
						for (int kw=w.getRowPtr(sigmaL);kw<w.getRowPtr(sigmaL+1);kw++) {
							SizeType sigmaLp=w.getCol(kw);
							for (int k2=r1.getRowPtr(alB);k2<r1.getRowPtr(alB+1);k2++) {
								SizeType alBp = r1.getCol(k2);
								SizeType j = 0;
								if (part_==PART_RIGHT) {
									SizeType tmp1 = symm.left().pack(alm1p,sigmaLp);
									j = symm.super().pack(tmp1,alBp);
								} else {
									SizeType tmp1 = symm.right().pack(sigmaLp,alBp);
									j = symm.super().pack(alm1p,tmp1);
								}
								if (j<offset || j>=offset+total) continue;
								x[i] += y[j-offset]*l1.getValue(k1)*
										w.getValue(kw)*r1.getValue(k2)*fermionSign;
							} // k2 right
						} // kw Hamiltonian
					} // k1 left
				} // symmetry sector
			} // bl
		} // blm1
	}

	//! Used only for stored option
	void fullHamiltonian(SparseMatrixType& matrix) const
	{
		SizeType nsites = symmetry_.symmLocal().size();
		SizeType middle = nsites/2;
		const SymmetryFactorType& symm =  symmetry_.symmLocal()(siteForSymm_);
		SizeType offset = symm.super().partitionOffset(symmetrySector_);
		SizeType total = symm.super().partitionSize(symmetrySector_);

		SizeType leftIndex = currentSite_;
		SizeType rightIndex = (currentSite_ < middle) ? currentSite_ + 1: nsites - currentSite_ - 1;

		matrix.resize(total,total);
		const ContractedFactorType& cL = contractedPart_(leftIndex,
														 ProgramGlobals::PART_LEFT);
		const ContractedFactorType& cR = contractedPart_(rightIndex,
														 ProgramGlobals::PART_RIGHT);
		VectorType v(total,0);
		SizeType counter = 0;
		assert(hamiltonian_.n_row()>=cL.size());
		assert(hamiltonian_.n_col()>=cR.size());

		for (SizeType i=0;i<total;i++) {
			matrix.setRow(i,counter);
			for (SizeType blm1=0;blm1<cL.size();blm1++) {
				for (SizeType bl=0;bl<cR.size();bl++) {
					const OperatorType& wOp = hamiltonian_(blm1,bl);
					SparseMatrixType w = wOp.matrix();
					const SparseMatrixType& r1 = cR(bl);
					const SparseMatrixType& l1 = cL(blm1);
					PairType ab = symm.super().unpack(i+offset);
					SizeType alm1=0;
					SizeType sigmaL=0;
					SizeType alB=0;

					SizeType electronsLeft = symmetry_.electronsFromQn(symm.left().
																	   qn(ab.first));
					if (part_==PART_RIGHT) {
						PairType tmpPair1 = symm.left().unpack(ab.first);
						alm1=tmpPair1.first;
						sigmaL=tmpPair1.second;
						alB = ab.second;
						SizeType electronsBlock = symmetry_.electronsFromState(sigmaL);
						assert(electronsBlock<=electronsLeft);
						electronsLeft -= electronsBlock;
					} else {
						PairType tmpPair1 = symm.right().unpack(ab.second);
						alB=tmpPair1.second;
						sigmaL=tmpPair1.first;
						alm1=ab.first;
					}

					RealType fermionSign = (electronsLeft & 1) ? wOp.fermionSign() : 1.0;

					for (int k1=l1.getRowPtr(alm1);k1<l1.getRowPtr(alm1+1);k1++) {
						SizeType alm1p=l1.getCol(k1);
						for (int kw=w.getRowPtr(sigmaL);kw<w.getRowPtr(sigmaL+1);kw++) {
							SizeType sigmaLp=w.getCol(kw);
							for (int k2=r1.getRowPtr(alB);k2<r1.getRowPtr(alB+1);k2++) {
								SizeType alBp = r1.getCol(k2);
								SizeType j = 0;
								if (part_==PART_RIGHT) {
									SizeType tmp1 =  symm.left().pack(alm1p,sigmaLp);
									j = symm.super().pack(tmp1,alBp);
								} else {
									SizeType tmp1 = symm.right().pack(sigmaLp,alBp);
									j = symm.super().pack(alm1p,tmp1);
								}

								if (j<offset || j>=offset+total) continue;
								v[j-offset] += l1.getValue(k1)*w.getValue(kw)*
										r1.getValue(k2)*fermionSign;
							} // k2 right
						} // kw Hamiltonian
					} // k1 left
				} // bl
			} // blm1
			for (SizeType j=0;j<v.size();j++) {
				if (fabs(v[j])<1e-12) continue;
				matrix.pushCol(j);
				matrix.pushValue(v[j]);
				v[j]=0;
				counter++;
			}
		} // symmetry sector
		matrix.setRow(total,counter);
		matrix.checkValidity();
	}

private:

	const ContractedPartType& contractedPart_;
	SizeType symmetrySector_;
	SizeType currentSite_;
	SizeType part_;
	SizeType hilbertSize_;
	const MpoFactorType& hamiltonian_;
	const SymmetryHelperType& symmetry_;
	SizeType siteForSymm_;

}; // ModelHelper

} // namespace Mpspp

/*@}*/
#endif // MODEL_HELPER_H
