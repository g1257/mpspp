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

#ifndef SYMMETRY_LOCAL_H
#define SYMMETRY_LOCAL_H

#include "ProgramGlobals.h"
#include "SymmetryFactor.h"

namespace Mpspp {

class SymmetryLocal {

public:

	typedef SymmetryFactor SymmetryFactorType;
	typedef SymmetryFactorType::SymmetryComponentType SymmetryComponentType;
	typedef SymmetryFactorType::PairType PairType;
	typedef SymmetryFactorType::IoInputType IoInputType;
	typedef SymmetryFactorType::VectorIntegerType VectorIntegerType;

	SymmetryLocal()
	{}

	const SymmetryFactorType& operator()(SizeType site) const
	{
		assert(site<data_.size());
		return data_[site];
	}

	const SizeType size() const
	{
		return data_.size();
	}

	void initialGuess(SizeType site,
	                  const VectorIntegerType& quantumNumbers,
	                  SizeType nsites)
	{
		SizeType middle = nsites/2;
		SymmetryFactorType symmFactor;
		SymmetryFactorType* ptr = (data_.size() == 0) ? 0 : &data_[data_.size()-1];
		if (site < middle) {
			symmFactor.grow(site,quantumNumbers,ptr,nsites);
		} else {
			SizeType ref = nsites - site -1;
			if (ref+1 < data_.size()) {
				SymmetryComponentType tmp(SymmetryComponentType::COMPONENT_RIGHT);
				if (ref>0) tmp = data_[ref-1].right();
				symmFactor.set(data_[ref+1].left(),tmp);
			} else {
				assert(data_.size() > 0);
				SymmetryComponentType onesite(SymmetryComponentType::COMPONENT_RIGHT,
				                              0,
				                              site,
				                              quantumNumbers);
				SizeType max = data_.size() - 1;
				SymmetryComponentType l(SymmetryComponentType::COMPONENT_LEFT);
				l.combine(data_[max].left(),onesite);
				symmFactor.set(l,data_[ref-1].right());
			}
		}

		data_.push_back(symmFactor);
		std::cout<<symmFactor;
	}

	template<typename SomeTruncationType>
	void truncate(SizeType site,
	              SizeType part,
	              SizeType cutoff,
	              const SomeTruncationType& trunc)
	{
		assert(site<data_.size());
		data_[site].truncate(part,cutoff,trunc);
	}

	friend std::ostream& operator<<(std::ostream& os,const SymmetryLocal& symm);

private:

	PsimagLite::Vector<SymmetryFactorType>::Type data_;

}; // SymmetryLocal
} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_LOCAL_H

