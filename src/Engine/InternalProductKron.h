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
///** \ingroup MPSPP */
///*@{*/

///*! \file InternalProductKron.h
// *
// *  A class to encapsulate the product x+=Hy, where x and y are vectors
// *  and H is the Hamiltonian matrix
// *
// */
//#ifndef	INTERNALPRODUCT_KRON_H
//#define INTERNALPRODUCT_KRON_H

//#include <vector>
//#include "InitKron.h"
//#include "KronMatrix.h"

//namespace Mpspp {
//template<typename T,typename ModelType>
//class InternalProductKron {
//public:
//	typedef T HamiltonianElementType;
//	typedef T value_type;
//	typedef typename ModelType::ModelHelperType ModelHelperType;
//	typedef typename ModelHelperType::RealType RealType;
//	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
//	typedef InitKron<ModelType,ModelHelperType> InitKronType;
//	typedef KronMatrix<InitKronType> KronMatrixType;

//	InternalProductKron(ModelType const *model,
//	                    ModelHelperType const *modelHelper,
//	                    ReflectionSymmetryType* rs=0)
//	    : initKron_(*model,*modelHelper),kronMatrix_(initKron_)
//	{
//	}

//	SizeType rank() const { return initKron_.size(); }

//	template<typename SomeVectorType>
//	void matrixVectorProduct(SomeVectorType &x,SomeVectorType const &y) const
//	{
//		kronMatrix_.matrixVectorProduct(x,y);
//	}

//	SizeType reflectionSector() const { return 0; }

//	void reflectionSector(SizeType p) {  }

//private:

//	InitKronType initKron_;
//	KronMatrixType kronMatrix_;
//}; // class InternalProductKron
//} // namespace Mpspp

///*@}*/
//#endif

