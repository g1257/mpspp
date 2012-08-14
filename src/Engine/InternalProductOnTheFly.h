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

/*! \file InternalProductOnTheFly.h
 *
 *  A class to encapsulate the product x+=Hy, where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef	INTERNALPRODUCT_OTF_H
#define INTERNALPRODUCT_OTF_H

#include <vector>

namespace Mpspp {
	template<typename T,typename ModelType>
	class InternalProductOnTheFly {
	public:
		typedef T HamiltonianElementType;
		typedef T value_type;
		typedef typename ModelType::ModelHelperType ModelHelperType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;

		InternalProductOnTheFly(ModelType const *model,
					ModelHelperType const *modelHelper,
					ReflectionSymmetryType* rs=0)
		{
			model_ = model;
			modelHelper_=modelHelper;
			
		}

		size_t rank() const { return modelHelper_->size(); }

		template<typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType &x,SomeVectorType const &y) const
		{
			 model_->matrixVectorProduct(x,y,*modelHelper_);
		}

		size_t reflectionSector() const { return 0; }

		void reflectionSector(size_t p) {  }

	private:
		ModelType const *model_;
		ModelHelperType const *modelHelper_;
	}; // class InternalProductOnTheFly
} // namespace Mpspp

/*@}*/
#endif

