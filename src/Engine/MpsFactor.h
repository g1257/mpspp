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

#ifndef MPS_FACTOR_TYPE_H
#define MPS_FACTOR_TYPE_H

#include "VectorWithOffset.h"
#include "ProgramGlobals.h"
#include "RandomForTests.h"
#include "Sort.h"

namespace Mpspp {

/* PSIDOC MpsFactor
		The main job of this class is to find the block of the
		Eigen-vector (in matrix format - M_{a_{l-1} \sigma_l \a_{l}})
		corresponding to the tarted quantum number of the superblock
		and perform singular value decomposition to new "A" or "B"
		operators, depending on the moveleft or moveright. If moving
		right, we perform SVD on M = USV^+ to obtain a new-optimised matrix
		A - which is the "U" of the SVD. Similarly, if moving left,
		same oparation to obtain a new-optimised result of "B" - which is
		the V^+ output of SVD. Note that, U and V^+ must be correctly
		normalized! MpsFactor is also the "typename" of matricies "A" and
		"B" at each site.
		*/

template<typename ComplexOrRealType,typename SymmetryFactorType>
class MpsFactor {

	typedef std::pair<SizeType,SizeType> PairType;

public:

	enum {TYPE_A,TYPE_B};

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef VectorWithOffset<ComplexOrRealType> VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename SymmetryFactorType::IoInputType IoInputType;
	typedef typename SymmetryFactorType::SymmetryComponentType SymmetryComponentType;
	typedef typename SymmetryComponentType::VectorIntegerType VectorIntegerType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::RandomForTests<RealType> RandomNumberGeneratorType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	MpsFactor(SizeType aOrB)
		: rng_(0),aOrB_(aOrB)
	{}

	void setRandom(SizeType site,SizeType n)
	{
		data_.resize(n,n);
		data_.makeDiagonal(n,1.0);
		//		assert(isNormalized(data_));
	}

	template<typename SomeTruncationType>
	void move(SomeTruncationType& truncation,
			  const VectorType& v,
			  SizeType symmetrySector,
			  const SymmetryFactorType& symm)
	{
		SizeType row = symm.left().size();
		SizeType col = symm.right().size();
		if (aOrB_==TYPE_B) std::swap(row,col);

		MatrixType m(row,col);

		SizeType offset = symm.super().partitionOffset(symmetrySector);
		SizeType total = symm.super().partitionSize(symmetrySector);
		SizeType qt = symm.super().qn(offset);
		for (SizeType i=0;i<total;i++) {
			PairType ab = symm.super().unpack(i+offset);
			SizeType qab = 0;
			if (symm.left().size() > 0) qab = symm.left().qn(ab.first);
			if (symm.right().size() > 0) qab += symm.right().qn(ab.second);
			assert(qab == qt);
			if (aOrB_==TYPE_A) {
				m(ab.first,ab.second) = v[i];
			} else {
				m(ab.second,ab.first) = v[i];
			}
		}

		moveFromVector(m,truncation,symm,qt);
	}

	template<typename SomeTruncationType>
	void truncate(SizeType cutoff,const SomeTruncationType& trunc)
	{
		trunc.matrixRow(data_,cutoff);
	}

	const SparseMatrixType& operator()() const { return data_; }

	const SizeType type() const { return aOrB_; }

private:

	template<typename SomeTruncationType>
	void moveFromVector(const MatrixType& m,
						SomeTruncationType& truncation,
						const SymmetryFactorType& symm,
						SizeType qt)
	{
		const SymmetryComponentType& summed = (aOrB_==TYPE_A) ? symm.left() : symm.right();
		const SymmetryComponentType& nonSummed = (aOrB_==TYPE_A) ? symm.right() : symm.left();

		MatrixType finalU(summed.size(),summed.size());
		assert(m.n_col() == nonSummed.size());
		assert(m.n_row() == summed.size());

#if 1 

		truncation.setSize(summed.size());
		for (SizeType i=0;i<summed.partitions()-1;i++) {
			SizeType istart = summed.partitionOffset(i);
			SizeType itotal = summed.partitionSize(i);
			SizeType qni = summed.qn(istart);
			for (SizeType j=0;j<nonSummed.partitions()-1;j++) {
				SizeType jstart = nonSummed.partitionOffset(j);
				SizeType jtotal = nonSummed.partitionSize(j);
				SizeType qnj = nonSummed.qn(jstart);

				if (qni + qnj != qt) continue;

				MatrixType u(itotal,jtotal);
				setThisSector(u,istart,itotal,jstart,jtotal,m);

				VectorRealType s;
				MatrixType vt;
				svd('A',u,s,vt);

				setFinalU(finalU,istart,itotal,jtotal,u);

				setFinalS(truncation,istart,itotal,s);
			}
		}
#else
		finalU = m;
		VectorRealType finalS(m.n_col());
		MatrixType finalVt(nonSummed.size(),nonSummed.size());
		svd('A',finalU,finalS,finalVt);
		resizeUAndNormalize(finalU,nonSummed.size());
		truncation.set(finalS);
#endif

		resizeU(finalU,std::min(summed.size(),nonSummed.size()));

		MatrixType mtranspose;
		if (aOrB_==TYPE_B)
			transposeConjugate(mtranspose,finalU);
		assert(isNormalized(finalU));
		//assert(respectsSymmetry(finalU,summed));
		//assert(isCorrectSvd(m,finalU,truncation,finalVt));
		fullMatrixToCrsMatrix(data_,(aOrB_==TYPE_A) ? finalU : mtranspose);
		assert(aOrB_ == TYPE_A);
	}

	void setThisSector(MatrixType& u,
					   SizeType istart,
					   SizeType itotal,
					   SizeType jstart,
					   SizeType jtotal,
					   const MatrixType& m) const
	{
		for (SizeType i=0;i<itotal;i++) {
			for (SizeType j=0;j<jtotal;j++) {
				u(i,j) = m(i+istart,j+jstart);
			}
		}
	}

	void setFinalU(MatrixType& finalU,
				   SizeType istart,
				   SizeType itotal,
				   SizeType jtotal,
				   const MatrixType& u) const
	{
		SizeType n = u.n_row();
		SizeType min = std::min(itotal,jtotal);
		assert(u.n_row()==u.n_col());
		for (SizeType i=0;i<n;i++) {
			for (SizeType j=0;j<min;j++) {
				finalU(i+istart,j+istart) = u(i,j);
			}
		}
	}

	template<typename SomeTruncationType>
	void setFinalS(SomeTruncationType& truncation,
				   SizeType jstart,
				   SizeType jtotal,
				   const VectorRealType& s) const
	{
		SizeType n = std::min(jtotal,static_cast<SizeType>(s.size()));
		for (SizeType j=0;j<n;j++) {
			//			assert(j+jstart<finalS.size());
			//			if (fabs(s[j])<1e-6) continue;
			if (j+jstart>=truncation.size()) continue;
			truncation(j+jstart) = s[j];
		}
	}


	bool isNormalized(const MatrixType& m) const
	{
		SizeType rows = m.n_row();
		SizeType cols = m.n_col();
		SizeType c = 0;
		for (SizeType i=0;i<cols;i++) {
			for (SizeType j=0;j<cols;j++) {
				ComplexOrRealType sum = 0;
				for (SizeType k=0;k<rows;k++) {
					sum += m(k,i) * PsimagLite::conj(m(k,j));
				}
				if (i==j && fabs(sum-1.0)>1e-5)
					c++;
				if (i!=j && fabs(sum)>1e-5)
					return false;
			}
		}

		return (c == 0);
	}

	void findNonZeroCols(VectorSizeType& vcols,
						 const MatrixType& m) const
	{
		SizeType rows = m.n_row();
		SizeType cols = m.n_col();
		SizeType c = 0;
		for (SizeType i=0;i<cols;++i) {
			ComplexOrRealType sum = 0;
			for (SizeType k=0;k<rows;++k)
				sum += m(k,i) * PsimagLite::conj(m(k,i));

			if (fabs(sum)<1e-6) continue;
			assert(c < vcols.size());
			vcols[c] = i;
			c++;
		}

		assert(c == vcols.size());
	}

	void resizeU(MatrixType& m,
				 SizeType smallSize)
	{
		if (m.n_col() == smallSize) return;
		SizeType rows = m.n_row();
		VectorSizeType vcols(smallSize,0);
		MatrixType tmp(rows,smallSize);
		findNonZeroCols(vcols,m);
		for (SizeType i=0;i<rows;++i) {
			for (SizeType j=0;j<smallSize;++j)
				tmp(i,j) = m(i,vcols[j]);
		}

		m = tmp;
	}

	RandomNumberGeneratorType rng_;
	SparseMatrixType data_;
	SizeType aOrB_;
}; // MpsFactor

} // namespace Mpspp

/*@}*/
#endif // MPS_FACTOR_TYPE_H

