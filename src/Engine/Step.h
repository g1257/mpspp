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

#ifndef STEP_H
#define STEP_H

//#include "ContractedLocal.h"
#include "ProgressIndicator.h"
#include "ParametersForSolver.h"
#include "LanczosOrDavidsonBase.h"
#include "LanczosSolver.h"
#include "DavidsonSolver.h"
#include <vector>
#include "StatePredictor.h"
#include "Truncation.h"
#include "FermionSign.h"

namespace Mpspp {


/* PSIDOC Step
		As the name suggests, contains functions required to physically
		add a site to the system/environment then build the superblock
		(LWR) and diagonalize using Lanczos to out the eigen-vector and
		the groun-state energy. Ground state eigen-vector is used to
		build the "M" matrix on which the SVD will be performed in
		MpsSolver.
		*/

template<typename ModelType,
		 template<typename,typename> class InternalProductTemplate>
class Step {

	typedef typename ModelType::MpoLocalType MpoLocalType;
	typedef typename MpoLocalType::MpsLocalType MpsLocalType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef typename ModelType::ParametersSolverType ParametersSolverType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ModelType::VectorType VectorType;
	typedef typename ModelType::ContractedLocalType ContractedLocalType;
	typedef typename ContractedLocalType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::SymmetryComponentType SymmetryComponentType;
	typedef StatePredictor<RealType,VectorType> StatePredictorType;
	typedef typename MpsLocalType::VectorRealType VectorRealType;
	typedef Truncation<ContractedLocalType> TruncationType;
	typedef typename MpsLocalType::VectorIntegerType VectorIntegerType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::SymmetryHelperType SymmetryHelperType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename ParametersSolverType::InputValidatorType InputValidatorType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT,
		  TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	static int const MAX_ = 100;

public:

	Step(const ParametersSolverType& solverParams,
		 MpsLocalType& mps,
		 ContractedLocalType& contractedLocal,
		 const ModelType& model,
		 InputValidatorType& io)
		: progress_("Step"),
		  solverParams_(solverParams),
		  mps_(mps),
		  contractedLocal_(contractedLocal),
		  model_(model),
		  paramsForSolver_(io,"Lanczos"),
		  statePredictor_(),
		  truncation_(mps,
					  contractedLocal,
					  solverParams_.options.find("notruncation")==PsimagLite::String::npos)
	{}

	//! Moves the center of orthogonality by one to the left
	void moveLeft(SymmetryLocalType& symm,
				  SizeType currentSite,
				  const FiniteLoop& finiteLoop)
	{
		SizeType nsites = model_.geometry().numberOfSites();
		if (currentSite == nsites) return;
		VectorIntegerType quantumNumbers;
		model_.getOneSite(quantumNumbers,currentSite);

		//symm.moveLeft(currentSite,quantumNumbers);
		if (currentSite==0) return;

		FermionSign<ModelType> fermionSign(model_,currentSite);
		SymmetryHelperType symmetryHelper(fermionSign,symm);
		internalmove(TO_THE_LEFT,symmetryHelper,nsites - currentSite - 1);
		contractedLocal_.move(currentSite,TO_THE_LEFT,symmetryHelper);
		truncation_(symm,currentSite,ProgramGlobals::PART_RIGHT,finiteLoop.keptStates);
	}

	//! Moves the center of orthogonality by one to the right
	void moveRight(SymmetryLocalType& symm,
				   SizeType currentSite,
				   const FiniteLoop& finiteLoop)
	{
		FermionSign<ModelType> fermionSign(model_,currentSite);
		SymmetryHelperType symmetryHelper(fermionSign,symm);
		internalmove(TO_THE_RIGHT,symmetryHelper,currentSite);
		contractedLocal_.move(currentSite,TO_THE_RIGHT,symmetryHelper);
		truncation_(symm,
					currentSite,
					ProgramGlobals::PART_LEFT,
					finiteLoop.keptStates);
	}

	void initialGuess(SymmetryLocalType& symm,SizeType center)
	{

		SizeType nsites = model_.geometry().numberOfSites();
		VectorIntegerType quantumNumbers;
		model_.getOneSite(quantumNumbers,center);
		for (SizeType i = 0; i< nsites; ++i){
			symm.initialGuess(i,quantumNumbers,nsites);
		}
		for (SizeType i = 0; i< nsites; ++i){
			mps_.initialGuess(i,symm,nsites);
		}

		for (SizeType i = 0; i<nsites-1;++i){

			FermionSign<ModelType> fermionSign(model_,i);
			SymmetryHelperType symmetryHelper(fermionSign,symm);
			contractedLocal_.initialGuess(i,symmetryHelper,nsites);
		}

	}

private:

	void internalmove(SizeType direction,
					  const SymmetryHelperType& symmetryHelper,
					  SizeType siteForSymm)
	{
		SizeType part = (direction == TO_THE_RIGHT) ? ProgramGlobals::PART_RIGHT:ProgramGlobals::PART_LEFT;
		const SymmetryFactorType& symm = symmetryHelper.symmLocal()(siteForSymm);
		SizeType currentSite = symmetryHelper.currentSite();
		SizeType symmetrySector = getSymmetrySector(direction,symm.super());
		std::cerr<<"symmetrySector="<<symmetrySector<<"\n";
		SizeType total = symm.super().partitionSize(symmetrySector);
		VectorType v(total,0.0);
		RealType energy = internalmove(v,
									   currentSite,
									   part,
									   symmetryHelper,
									   symmetrySector,
									   siteForSymm);
		mps_.move(truncation_,currentSite,v,direction,symmetrySector,symm);
		statePredictor_.push(energy,v,symmetrySector);
		if (solverParams_.options.find("test")!=PsimagLite::String::npos)
			throw PsimagLite::LogicError("Exiting due to option test in the input\n");
	}

	RealType internalmove(VectorType& tmpVec,
						  SizeType currentSite,
						  SizeType part,
						  const SymmetryHelperType& symmetryHelper,
						  SizeType symmetrySector,
						  SizeType siteForSymm)
	{
		const ParametersSolverType& solverParams = model_.solverParams();

		typedef InternalProductTemplate<typename VectorType::value_type,ModelType>
				InternalProductType;
		typedef PsimagLite::LanczosOrDavidsonBase<ParametersForSolverType,
				InternalProductType,
				VectorType> LanczosOrDavidsonBaseType;

		ReflectionSymmetryType *rs = 0;
		ModelHelperType modelHelper(contractedLocal_,
									symmetrySector,
									currentSite,
									part,
									model_.hamiltonian()(currentSite),
									symmetryHelper,
									siteForSymm);
		InternalProductType lanczosHelper(&model_,&modelHelper,rs);

		LanczosOrDavidsonBaseType* lanczosOrDavidson = 0;

		bool useDavidson =
				(solverParams.options.find("useDavidson")!=PsimagLite::String::npos);
		if (useDavidson) {
			lanczosOrDavidson = new PsimagLite::DavidsonSolver<ParametersForSolverType,
					InternalProductType,
					VectorType>(lanczosHelper,paramsForSolver_);
		} else {
			lanczosOrDavidson = new PsimagLite::LanczosSolver<ParametersForSolverType,
					InternalProductType,
					VectorType>(lanczosHelper,paramsForSolver_);
		}

		RealType energyTmp = 0;
		assert(tmpVec.size()==lanczosHelper.rank());
		VectorType initialVector(lanczosHelper.rank());
		statePredictor_.createRandomVector(initialVector,0,initialVector.size());

		lanczosOrDavidson->computeGroundState(energyTmp,tmpVec,initialVector);
		if (lanczosOrDavidson) delete lanczosOrDavidson;
		return energyTmp;
	}

	SizeType getSymmetrySector(SizeType direction,
							   const SymmetryComponentType& super) const
	{
		SizeType sites = super.block().size();
		SizeType targetQuantumNumber = getQuantumSector(sites,direction);
		SizeType imin=0;
		SizeType minDiff=0;
		for (SizeType i=0;i<super.partitions()-1;i++) {
			SizeType state = super.partitionOffset(i);
			SizeType q = super.qn(state);
			SizeType diff = (q<targetQuantumNumber) ?
						targetQuantumNumber - q : q-targetQuantumNumber;
			if (diff<minDiff || i==0) {
				imin = i;
				minDiff = diff;
			}
		}

		if (minDiff!=0) {
			std::cerr<<__FILE__<<" "<<__LINE__<<"\n";
			std::cerr<<"getSymmetrySector ";
			std::cerr<<"WARNING: minDiff="<<minDiff<<" is non zero\n";
			throw std::runtime_error("error\n");
		}

		return imin;
	}

	SizeType getQuantumSector(SizeType sites,SizeType direction) const
	{
		return (solverParams_.targetQuantumNumbers.size()>0) ?
					getQuantumSectorT(sites,direction) :
					getQuantumSectorUd(sites,direction);
	}

	SizeType getQuantumSectorT(SizeType sites,SizeType direction) const
	{
		VectorIntegerType targetQuantumNumbers(solverParams_.targetQuantumNumbers.size());
		for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
			targetQuantumNumbers[ii] =
					SizeType(round(solverParams_.targetQuantumNumbers[ii]*sites));
		return getQuantumSector(targetQuantumNumbers,direction);
	}

	SizeType getQuantumSectorUd(SizeType sites,SizeType direction) const
	{
		VectorIntegerType targetQuantumNumbers(2);

		SizeType nsites = model_.geometry().numberOfSites();
		SizeType total = solverParams_.electronsUp + solverParams_.electronsDown;
		total *= sites;
		total /= nsites;
		targetQuantumNumbers[0] =
				SizeType(static_cast<RealType>(solverParams_.electronsUp*sites)/nsites);
		targetQuantumNumbers[1] =
				SizeType(static_cast<RealType>(solverParams_.electronsDown*sites)/nsites);

		if (targetQuantumNumbers[0] + targetQuantumNumbers[1] > total)
			targetQuantumNumbers[0]--;

		if (targetQuantumNumbers[0] + targetQuantumNumbers[1] < total)
			targetQuantumNumbers[0]++;

		return getQuantumSector(targetQuantumNumbers,direction);
	}

	SizeType getQuantumSector(const VectorIntegerType& targetQuantumNumbers,
							  SizeType direction) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Integer target quantum numbers are: ";
		for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
			msg<<targetQuantumNumbers[ii]<<" ";
		progress_.printline(msg,std::cout);
		return encodeQuantumNumber(targetQuantumNumbers);
	}

	static SizeType encodeQuantumNumber(const VectorIntegerType& v)
	{
		SizeType x= v[0] + v[1]*MAX_;
		if (v.size()==3) x += v[2]*MAX_*MAX_;
		return x;
	}

	PsimagLite::ProgressIndicator progress_;
	const ParametersSolverType& solverParams_;
	MpsLocalType& mps_;
	ContractedLocalType& contractedLocal_;
	const ModelType& model_;
	ParametersForSolverType paramsForSolver_;
	StatePredictorType statePredictor_;
	TruncationType truncation_;
}; // Step

} // namespace Mpspp

/*@}*/
#endif // STEP_H

