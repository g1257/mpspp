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

#ifndef CONTRACTED_FACTOR_H
#define CONTRACTED_FACTOR_H

#include "ProgramGlobals.h"

namespace Mpspp {

/* PSIDOC ContractedFactor
		ContractedFactor are type name of the L_ and R_ (eq. 192 and 193).
		This also contains the ``guts'' of moving "right" or "left" during
		the finite/infinite MPS sweep.
		*/

template<typename MatrixProductOperatorType>
class ContractedFactor {

	typedef typename MatrixProductOperatorType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::VectorType VectorType;
	typedef typename MpsLocalType::SparseMatrixType SparseMatrixType;
	typedef typename MpsLocalType::MpsFactorType  MpsFactorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MatrixProductOperatorType::SymmetryHelperType SymmetryHelperType;
	typedef typename MpsLocalType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef ContractedFactor<MatrixProductOperatorType> ThisType;
	typedef PsimagLite::Matrix<ComplexOrRealType> DenseMatrixType;
	typedef typename MpsLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename MpsFactorType::VectorIntegerType VectorIntegerType;
	typedef typename MpoFactorType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairForOperatorType;
	typedef typename SymmetryFactorType::SymmetryComponentType SymmetryComponentType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT,
		  TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT,
		  PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:

	typedef typename PsimagLite::Vector<SparseMatrixType>::Type DataType;

	ContractedFactor(SizeType leftOrRight)
	    : data_(1),leftOrRight_(leftOrRight)
	{
		data_[0].makeDiagonal(1,1);
	}

	//! NEEDED
	void build(const MpsFactorType& AorB,
	           const MpoFactorType& h,
	           const ThisType* prev,
	           const SymmetryHelperType& symmHelper,
	           SizeType siteForSymm)
	{
		const DataType* ptr = (prev) ? &prev->data_: 0;
		if (leftOrRight_==PART_RIGHT) {
			assert(AorB.type()==MpsFactorType::TYPE_B);

			data_.resize(h.n_row());

			moveRight(AorB,h,ptr,symmHelper,siteForSymm,true);
		} else {
			assert(AorB.type()==MpsFactorType::TYPE_A);

			data_.resize(h.n_col());
			moveLeft(AorB,h,ptr,symmHelper,siteForSymm);
		}
	}

	//! From As (or Bs) and Ws reconstruct *this
	//! NEEDED
	void move(const MpsFactorType& AorB,
	          const MpoFactorType& h,
	          const ThisType& dataPrev,
	          const SymmetryHelperType& symm,
	          SizeType currentSite)
	{
		if (leftOrRight_ == PART_RIGHT) {
			assert(AorB.type()==MpsFactorType::TYPE_B);
			moveRight(AorB,h,&dataPrev.data_,symm,currentSite,false);
		} else {
			assert(AorB.type()==MpsFactorType::TYPE_A);
			moveLeft(AorB,h,&dataPrev.data_,symm,currentSite);
		}
	}


	//!NEEDED
	template<typename SomeTruncationType>
	void truncate(SizeType cutoff,const SomeTruncationType& trunc)
	{
		for (SizeType i=0;i<data_.size();i++)
			trunc.matrixRowCol(data_[i],cutoff);
	}

	//!NEEDED
	const SparseMatrixType& operator()(SizeType b1) const
	{
		assert(b1<data_.size());
		return data_[b1];
	}

	SizeType size() const { return data_.size(); }

	SizeType row() const
	{
		assert(data_.size()>0);
		return data_[0].row();
	}

private:

	/* MOVE LEFT WILL BE NEEDED LATER ------------------------- */
	void moveLeft(SparseMatrixType& m,
	              const MpsFactorType& AorB,
	              const SparseMatrixType& Atranspose,
	              SizeType b,
	              const MpoFactorType& h,
	              const DataType* dataPrevPtr,
	              const SymmetryHelperType& symmHelper,
	              SizeType currentSite)
	{
		SizeType siteForSymm = currentSite;

		const SymmetryFactorType& symm = symmHelper.symmLocal()(siteForSymm);
		const SparseMatrixType& A = AorB();
		const DataType& dataPrev = *dataPrevPtr;

		m.resize(A.col(),A.col());
		SizeType counter=0;
		VectorIntegerType cols(m.row(),0);
		VectorType values(m.row(),0.0);
		assert(symm.left().split()==0 || symm.left().split()==dataPrev[0].row());
		assert(symm.left().size()==A.row());
		for (SizeType a2=0;a2<Atranspose.row();a2++) {
			m.setRow(a2,counter);
			for (int k3=Atranspose.getRowPtr(a2);k3<Atranspose.getRowPtr(a2+1);k3++) {
				PairType a1sigma2 = symm.left().unpack(Atranspose.getCol(k3));
				SizeType a1 = a1sigma2.first;
				SizeType sigma2 = a1sigma2.second;

				for (SizeType b1=0;b1<h.n_row();b1++) {
					const OperatorType& wOp = h(b1,b);
					const SparseMatrixType& w = wOp.matrix();
					if (w.row()==0) continue;
					const SparseMatrixType& l1 = dataPrev[b1];

					for (int k=l1.getRowPtr(a1);k<l1.getRowPtr(a1+1);k++) {
						SizeType a1p = l1.getCol(k);
						ComplexOrRealType tmp= l1.getValue(k)*
						        PsimagLite::conj(Atranspose.getValue(k3));
						for (int kp=w.getRowPtr(sigma2);kp<w.getRowPtr(sigma2+1);kp++) {
							SizeType sigma2p = w.getCol(kp);
							SizeType j = symm.left().pack(a1p,sigma2p);
							for (int k2=A.getRowPtr(j);k2<A.getRowPtr(j+1);k2++) {
								SizeType a2p = A.getCol(k2);
								values[a2p] += tmp*w.getValue(kp)*A.getValue(k2);
								cols[a2p]=1;
							} // k2
						} // kp
					} // k
				} //b1
			} // k3
			for (SizeType i=0;i<cols.size();i++) {
				if (cols[i]==0) continue;
				cols[i]=0;
				m.pushCol(i);
				m.pushValue(values[i]);
				values[i]=0;
				counter++;
			}
		} // a2

		m.setRow(m.row(),counter);
		m.checkValidity();
	}

	void moveLeft(const MpsFactorType& A,
	              const MpoFactorType& h,
	              const DataType* dataPrev,
	              const SymmetryHelperType& symm,
	              SizeType currentSite)
	{
		assert(leftOrRight_ == PART_LEFT);
		assert(A.type()==MpsFactorType::TYPE_A);
		SparseMatrixType Atranspose;
		transposeConjugate(Atranspose,A());
		if (data_.size()!=h.n_col())
			data_.resize(h.n_col());
		assert(data_.size()==h.n_col());
		for (SizeType b1=0;b1<data_.size();b1++)
			moveLeft(data_[b1],A,Atranspose,b1,h,dataPrev,symm,currentSite);
	}

	void moveRight(const MpsFactorType& B,
	               const MpoFactorType& h,
	               const DataType* dataPrev,
	               const SymmetryHelperType& symm,
	               SizeType siteForSymm,
	               bool isInitialGuess)
	{
		assert(leftOrRight_ == PART_RIGHT);
		assert(B.type()==MpsFactorType::TYPE_B);
		SparseMatrixType Btranspose;
		transposeConjugate(Btranspose,B());
		if (h.n_row()!=data_.size())
			data_.resize(h.n_row());
		assert(h.n_row()==data_.size());
		for (SizeType blm2=0;blm2<data_.size();blm2++)
			moveRight(data_[blm2],blm2,B,Btranspose,h,dataPrev,symm,siteForSymm,isInitialGuess);
	}

	void moveRight(SparseMatrixType& m,
	               SizeType blm2,
	               const MpsFactorType& B,
	               const SparseMatrixType& Btranspose,
	               const MpoFactorType& h,
	               const DataType* dataPrev,
	               const SymmetryHelperType& symm,
	               SizeType currentSite,
	               bool isInitialGuess)
	{
		SizeType nsites = symm.symmLocal().size();
		SizeType middle = static_cast<SizeType>(nsites/2);

		const DataType* dataPrevCopy = (dataPrev && dataPrev->size() == 1) ? 0 : dataPrev;

		if (currentSite >= middle && dataPrevCopy && isInitialGuess) {
			const SymmetryComponentType& symmC =
			        symm.symmLocal()(nsites-currentSite-1).right();
			SizeType someSize = B().row()/symmC.split();
			m.resize(someSize,someSize);

			SizeType counter = 0;
			VectorIntegerType cols(m.row(),0);
			VectorType values(m.row(),0.0);

			for (SizeType alm2=0;alm2<m.row();alm2++) {
				m.setRow(alm2,counter);

				for (SizeType sigma = 0; sigma < symmC.split(); ++sigma){

					SizeType sigmalm1alm1 = symmC.pack(sigma,alm2);

					moveRightRight(values,
					               cols,
					               sigmalm1alm1,
					               blm2,
					               B,
					               Btranspose,
					               h,
					               *dataPrevCopy,
					               symm,
					               currentSite);
				}

				for (SizeType i=0;i<cols.size();i++) {
					if (cols[i]==0) continue;
					cols[i]=0;
					m.pushCol(i);
					m.pushValue(values[i]);
					values[i]=0;
					counter++;
				}
			}
			m.setRow(m.row(),counter);
			m.checkValidity();
			return;
		}

		SizeType someSize = B().row();
		m.resize(someSize,someSize);

		SizeType counter = 0;
		VectorIntegerType cols(m.row(),0);
		VectorType values(m.row(),0.0);

		SizeType siteForSymm = (currentSite<middle) ? currentSite : nsites-currentSite-1;

		for (SizeType alm2=0;alm2<m.row();alm2++) {
			m.setRow(alm2,counter);
			moveRight(values,cols,alm2,blm2,B,Btranspose,h,dataPrevCopy,symm,siteForSymm);
			for (SizeType i=0;i<cols.size();i++) {
				if (cols[i]==0) continue;
				cols[i]=0;
				m.pushCol(i);
				m.pushValue(values[i]);
				values[i]=0;
				counter++;
			}
		}

		m.setRow(m.row(),counter);
		m.checkValidity();
	}


	void moveRight(VectorType& values,
	               VectorIntegerType& cols,
	               SizeType alm2,
	               SizeType blm2,
	               const MpsFactorType& B,
	               const SparseMatrixType& Btranspose,
	               const MpoFactorType& h,
	               const DataType* dataPrevPtr,
	               const SymmetryHelperType& symmHelper,
	               SizeType siteForSymm)
	{
		if (!dataPrevPtr)
			return moveRightFirst(values,cols,alm2,blm2,B,Btranspose,h,symmHelper,siteForSymm);

		const SymmetryFactorType& symm = symmHelper.symmLocal()(siteForSymm);
		const SymmetryComponentType& symmC = symm.right();
		const SparseMatrixType& Bmatrix = B();
		const DataType& dataPrev = *dataPrevPtr;

		assert(symmC.split()==0 || symmC.size()/symmC.split()==dataPrev[0].row());
		assert(Btranspose.row()==symmC.size());
		assert(dataPrev.size()<=h.n_col());

		for (int kb=Bmatrix.getRowPtr(alm2);kb<Bmatrix.getRowPtr(alm2+1);kb++) {
			PairType sigmalm1alm1 = symmC.unpack(Bmatrix.getCol(kb));
			SizeType sigmalm1 = sigmalm1alm1.first;
			SizeType alm1 = sigmalm1alm1.second;

			for (SizeType blm1=0;blm1<dataPrev.size();blm1++) {
				const OperatorType& wOp = h(blm2,blm1);
				const SparseMatrixType& w = wOp.matrix();
				if (w.row()==0) continue;
				RealType fermionSign = 1.0;
				for (int k=w.getRowPtr(sigmalm1);k<w.getRowPtr(sigmalm1+1);k++) {
					SizeType sigmaplm1 = w.getCol(k);
					ComplexOrRealType tmp = PsimagLite::conj(Bmatrix.getValue(kb))*w.getValue(k);
					for (int kprev=dataPrev[blm1].getRowPtr(alm1);
					     kprev<dataPrev[blm1].getRowPtr(alm1+1);
					     kprev++) {
						SizeType aplm1 = dataPrev[blm1].getCol(kprev);
						SizeType i = symmC.pack(sigmaplm1,aplm1);
						for (int kb2=Btranspose.getRowPtr(i);
						     kb2<Btranspose.getRowPtr(i+1);
						     kb2++) {
							SizeType aplm2 = Btranspose.getCol(kb2);
							values[aplm2] += tmp*Btranspose.getValue(kb2)*
							        dataPrev[blm1].getValue(kprev)*fermionSign;
							cols[aplm2]=1;
						} // kb2
					} // kprev
				} // k
			} // blm1
		} // kb
	}

	//! NEEDED
	void moveRightRight(VectorType& values,
	                    VectorIntegerType& cols,
	                    SizeType alm2,
	                    SizeType blm2,
	                    const MpsFactorType& B,
	                    const SparseMatrixType& Btranspose,
	                    const MpoFactorType& h,
	                    const DataType& dataPrev,
	                    const SymmetryHelperType& symmHelper,
	                    SizeType currentSite)
	{

		SizeType nsites = symmHelper.symmLocal().size();
		SizeType middle = static_cast<SizeType>(nsites/2);
		SizeType siteForSymm = (currentSite<middle) ? currentSite : nsites-currentSite-1;
		assert(currentSite >= middle);
		const SymmetryFactorType& symm = symmHelper.symmLocal()(siteForSymm);
		const SymmetryComponentType& symmC = symm.right();
		const SparseMatrixType& Bmatrix = B();

		assert(symmC.split()==0 ||
		       symmC.size()==dataPrev[0].row());
		assert(Btranspose.col()==symmC.size());
		assert(dataPrev.size()<=h.n_row());

		for (int kb=Bmatrix.getRowPtr(alm2);kb<Bmatrix.getRowPtr(alm2+1);kb++) {

			SizeType a = Bmatrix.getCol(kb);
			SizeType sigma = symmC.unpack(alm2).first;
			assert(sigma<symmC.split());
			for (SizeType blm1=0;blm1<dataPrev.size();blm1++) {
				const OperatorType& wOp = h(blm2,blm1);
				const SparseMatrixType& w = wOp.matrix();
				if (w.row()==0) continue;
				RealType fermionSign = 1.0;
				for (int k=w.getRowPtr(sigma);k<w.getRowPtr(sigma+1);k++) {
					SizeType sigmap = w.getCol(k);
					ComplexOrRealType tmp = PsimagLite::conj(Bmatrix.getValue(kb))*w.getValue(k);
					for (int kprev=dataPrev[blm1].getRowPtr(a);
					     kprev<dataPrev[blm1].getRowPtr(a+1);
					     kprev++) {
						SizeType ap = dataPrev[blm1].getCol(kprev);
						//SizeType i = symmC.pack(sigmap,ap);
						for (int kb2=Btranspose.getRowPtr(ap);
						     kb2<Btranspose.getRowPtr(ap+1);
						     kb2++) {
							SizeType col = Btranspose.getCol(kb2);
							PairType apbar = symmC.unpack(col);
							assert(apbar.first < symmC.split());
							if (apbar.first != sigmap) continue;
							values[apbar.second] += tmp*Btranspose.getValue(kb2)*
							        dataPrev[blm1].getValue(kprev)*fermionSign;
							cols[apbar.second]=1;
						} // kb2
					} // kprev
				} // k
			} // blm1
		} // kb
	}


	/* PSIDOC ContractedFactor - moveRightFirst
			First step of moving right. It is done if "previous" mps step is not defined.
	*/
	//! NEEDED
	void moveRightFirst(VectorType& values,
	                    VectorIntegerType& cols,
	                    SizeType alm2,
	                    SizeType blm2,
	                    const MpsFactorType& B,
	                    const SparseMatrixType& Btranspose,
	                    const MpoFactorType& h,
	                    const SymmetryHelperType& symmHelper,
	                    SizeType)
	{
		//SizeType nsites = symmHelper.symmLocal().size();
		//SizeType middle = static_cast<SizeType>(nsites/2);
		SizeType siteForSymm = 0;
		const SymmetryFactorType& symm = symmHelper.symmLocal()(siteForSymm);
		const SymmetryComponentType& symmC = symm.right();

		const SparseMatrixType& Bmatrix = B();

		assert(Btranspose.row()==symmC.size());

		for (int kb=Bmatrix.getRowPtr(alm2);kb<Bmatrix.getRowPtr(alm2+1);kb++) {
			PairType sigmalm1alm1 = symmC.unpack(Bmatrix.getCol(kb));
			SizeType sigmalm1 = sigmalm1alm1.first;
			assert(sigmalm1alm1.second == 0);
			const OperatorType& wOp = h(blm2,0);
			const SparseMatrixType& w = wOp.matrix();
			if (w.row()==0) continue;
			RealType fermionSign = 1.0;
			for (int k=w.getRowPtr(sigmalm1);k<w.getRowPtr(sigmalm1+1);k++) {
				SizeType sigmaplm1 = w.getCol(k);
				ComplexOrRealType tmp = PsimagLite::conj(Bmatrix.getValue(kb))*w.getValue(k);

				for (int kb2=Bmatrix.getRowPtr(sigmaplm1);
				     kb2<Bmatrix.getRowPtr(sigmaplm1+1);
				     kb2++) {
					SizeType aplm2 = Bmatrix.getCol(kb2);
					values[aplm2] += tmp*Bmatrix.getValue(kb2)*fermionSign;
					cols[aplm2]=1;
				} // kb2
			} // k
		} // kb
	}

	DataType data_;
	SizeType leftOrRight_;

}; // ContractedFactor

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_FACTOR_H

