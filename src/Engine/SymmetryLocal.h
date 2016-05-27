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

	//	SymmetryLocal(IoInputType& io)
	//	{
	//		SizeType n = 0;
	//		io.readline(n,"TotalNumberOfSites=");
	//		SizeType nk = 0;
	//		io.readline(nk,"HilbertOneSite=");
	//		assert(n>2);
	//		for (SizeType i=0;i<n-1;i++) {
	//			SymmetryFactorType f(io,nk);
	//			data_.push_back(f);
	////			if (i==0) data_.push_back(f); // left corner
	////			if (i==n-3) data_.push_back(f); // right corner
	//		}

	//		assert(data_.size()==n-1);
	////		data_[0].adjustCorner(SymmetryFactorType::CORNER_LEFT);
	////		data_[n-1].adjustCorner(SymmetryFactorType::CORNER_RIGHT);
	//	}

	const SymmetryFactorType& operator()(SizeType site) const
	{
		assert(site<data_.size());
		return data_[site];
	}

	void moveLeft(SizeType site,const VectorIntegerType& quantumNumbers)
	{
		if (site+1==data_.size()) return;
		assert(site+1<data_.size());
		SymmetryFactorType symmFactor = data_[site];
		SymmetryComponentType onesite(SymmetryComponentType::COMPONENT_LEFT,
		                              0,
		                              site,
		                              quantumNumbers);
		assert(site+1<data_.size());
		symmFactor.moveLeft(data_[site].left(),onesite, data_[site+1].right());
		data_[site] = symmFactor;
		std::cout<<symmFactor;
	}

	void moveRight(SizeType site,const VectorIntegerType& quantumNumbers)
	{
		if (site==0) return;
		assert(site>0);
		VectorIntegerType qn(1,0);

		if (site==data_.size()) {
			SizeType nsites = data_[site-1].super().block().size();
			assert(site<=nsites);
			assert(nsites-site<data_.size());
			data_.push_back(data_[nsites-site]);
			SymmetryComponentType onesite2(SymmetryComponentType::COMPONENT_RIGHT,0,site,qn);
			if (site==nsites)
				data_[data_.size()-1].set(SymmetryComponentType::COMPONENT_RIGHT,onesite2);
		}
		assert(site<data_.size());
		SymmetryFactorType symmFactor;
		SymmetryComponentType onesite(SymmetryComponentType::COMPONENT_RIGHT,
		                              0,
		                              site-1,
		                              quantumNumbers);
		symmFactor.moveRight(data_[site-1].left(),onesite,data_[site].right());

		data_[site] = symmFactor;
		std::cout<<symmFactor;
	}

    void initialGuess(SizeType site,const VectorIntegerType& quantumNumbers,SizeType nsites)
    {
            SymmetryFactorType symmFactor;
            SymmetryFactorType* ptr = (data_.size() == 0) ? 0 : &data_[data_.size()-1];
            symmFactor.grow(site,quantumNumbers,ptr,nsites);
            data_.push_back(symmFactor);
            std::cout<<symmFactor;
    }

	// left = prev.left + one site
	// right = prev.right + one site
	void grow(SizeType site,const VectorIntegerType& quantumNumbers,SizeType nsites)
	{

		if (data_.size()==0) {
			SymmetryFactorType symmFactor0;
			symmFactor0.growFirst(site,quantumNumbers,nsites);
			data_.push_back(symmFactor0);
		}
		SymmetryFactorType symmFactor;
        symmFactor.grow(site,quantumNumbers,&data_[data_.size()-1],nsites);
		data_.push_back(symmFactor);
		std::cout<<symmFactor;
	}

	template<typename SomeTruncationType>
	void truncate(SizeType site,SizeType part,SizeType cutoff,const SomeTruncationType& trunc)
	{
		assert(site<data_.size());
		data_[site].truncate(part,cutoff,trunc);
	}

	friend std::ostream& operator<<(std::ostream& os,const SymmetryLocal& symm);

private:

	PsimagLite::Vector<SymmetryFactorType>::Type data_;

}; // SymmetryLocal

std::ostream& operator<<(std::ostream& os,const SymmetryLocal& symm)
{
	os<<"symm.data.size= "<<symm.data_.size()<<"\n";
	for (SizeType i=0;i<symm.data_.size();i++)
		os<<symm.data_[i];
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_LOCAL_H

