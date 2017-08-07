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

#ifndef TRUNCATION_H
#define TRUNCATION_H

#include "AllocatorCpu.h"
#include "Sort.h"
#include "ProgramGlobals.h"

namespace Mpspp {

template<typename ContractedLocalType>
class Truncation {

	typedef typename ContractedLocalType::MatrixProductOperatorType MpoLocalType;
	typedef typename MpoLocalType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::MatrixType MatrixType;
	typedef typename MpsLocalType::SparseMatrixType SparseMatrixType;
	typedef typename MpsLocalType::VectorRealType VectorRealType;
	typedef typename VectorRealType::value_type RealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorIntegerType;
	typedef typename ContractedLocalType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;

	enum {PERMUTE_ROW=1, PERMUTE_COL=2};

public:

	Truncation(MpsLocalType& mps,
	           ContractedLocalType& contracted,
	           bool enabled,
	           SizeType numberOfSites)
	    : mps_(mps),contracted_(contracted),enabled_(enabled),numberOfSites_(numberOfSites)
	{}

	void setSize(SizeType size)
	{
		s_.resize(size);
		for (SizeType i=0;i<s_.size();i++)
			s_[i] = 0.0;
	}

	SizeType size() const { return s_.size(); }

	RealType& operator()(SizeType i)
	{
		assert(i<s_.size());
		return s_[i];
	}

	void operator()(SymmetryLocalType& symm,SizeType site,SizeType part,SizeType cutoff)
	{
		if (!enabled_) return;
		SizeType siteForSymm = (part==ProgramGlobals::PART_LEFT) ? site+1 : site;
		SizeType rightSize = symm(siteForSymm).right().size();
		SizeType leftSize  = symm(siteForSymm).left().size();

		if (part==ProgramGlobals::PART_LEFT) {
			if (leftSize<=cutoff) return;
		} else {
			if (rightSize<=cutoff) return;
		}

		order();
		SizeType nsites = numberOfSites_;
		mps_.truncate(site,part,cutoff,nsites,*this);
		contracted_.truncate(site,part,cutoff,nsites,*this);

		symm.truncate(siteForSymm,part,cutoff,*this);
	}

	void vector(VectorIntegerType& quantumNumbers,SizeType cutoff) const
	{
		cutoff = std::min(cutoff,static_cast<SizeType>(quantumNumbers.size()));
		SizeType toRemove =  quantumNumbers.size()-cutoff;
		VectorIntegerType q(cutoff);

		assert(quantumNumbers.size()==perm_.size());
		for (SizeType i=0;i<quantumNumbers.size();i++) {
			if (perm_[i]<toRemove) continue;
			q[perm_[i]-toRemove] = quantumNumbers[i];
		}
		quantumNumbers = q;
	}

	void matrixRow(SparseMatrixType& m,SizeType cutoff) const
	{
		cutoff = std::min(cutoff,m.cols());
		assert(m.cols()==perm_.size());
		SizeType toRemove = m.cols()-cutoff;
		SparseMatrixType newmatrix;
		permute(newmatrix,m, PERMUTE_COL);
		SparseMatrixType newdata(m.rows(),cutoff);
		SizeType counter = 0;
		for (SizeType i=0;i<newmatrix.rows();i++) {
			newdata.setRow(i,counter);
			for (int k=newmatrix.getRowPtr(i);k<newmatrix.getRowPtr(i+1);k++) {
				SizeType col = newmatrix.getCol(k);
				if (col<toRemove) continue;
				newdata.pushCol(col-toRemove);
				newdata.pushValue(m.getValue(k));
				counter++;
			}
		}
		newdata.setRow(newmatrix.rows(),counter);
		m=newdata;
		m.checkValidity();
	}

	void matrixRowCol(SparseMatrixType& m,SizeType cutoff) const
	{
		assert(m.rows()==m.cols());
		cutoff = std::min(cutoff,m.rows());
		assert(m.cols()==perm_.size());
		SparseMatrixType newmatrix;
		permute(newmatrix,m,PERMUTE_ROW | PERMUTE_COL);
		SizeType toRemove = m.cols()-cutoff;
		MatrixType dest(cutoff,cutoff);
		for (SizeType i=0;i<newmatrix.rows();i++) {
			if (i<toRemove) continue;
			for (int k=newmatrix.getRowPtr(i);k<newmatrix.getRowPtr(i+1);k++) {
				SizeType j = newmatrix.getCol(k);
				if (j<toRemove) continue;
				dest(i-toRemove,j-toRemove) = newmatrix.getValue(k);
			}
		}
		fullMatrixToCrsMatrix(m,dest);
		m.checkValidity();
	}

	void order()
	{
		perm_.resize(s_.size());
		PsimagLite::Sort<VectorRealType> sort;
		sort.sort(s_,perm_);
	}

private:

	void permute(SparseMatrixType& m,const SparseMatrixType& src,SizeType what) const
	{
		const VectorIntegerType& perm = perm_;
		SizeType row = src.rows();
		MatrixType dest(row,src.cols());
		for (SizeType i=0;i<row;i++) {
			SizeType ind = (what & PERMUTE_ROW) ? perm[i] : i;
			for (int k=src.getRowPtr(i);k<src.getRowPtr(i+1);k++) {
				SizeType j = src.getCol(k);
				SizeType jnd = (what & PERMUTE_COL) ? perm[j] : j;
				dest(ind,jnd) = src.getValue(k);
			}
		}
		fullMatrixToCrsMatrix(m,dest);
	}

	MpsLocalType& mps_;
	ContractedLocalType& contracted_;
	bool enabled_;
	SizeType numberOfSites_;
	VectorRealType s_;
	VectorIntegerType perm_;

}; // Truncation

} // namespace Mpspp

/*@}*/
#endif // TRUNCATION_H

