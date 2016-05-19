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

/*! \file ParametersMpsSolver.h
 *
 *  Contains the parameters for the DmrgSolver class and implements functionality to
 *  read them from a JSON file
 *
 */
#ifndef ParametersMpsSolver_HEADER_H
#define ParametersMpsSolver_HEADER_H

#include "TypeToString.h"
#include "Vector.h"
#include "Provenance.h"

namespace Mpspp {
/**
 \\subsubsection{Enabling finite loops}
 \\begin{itemize}
 \\item
 {\\bf options line:} \\verb=nofiniteloops=: Don't do finite loops, even if provided under {\\bf FiniteLoops} below.
 \\item {\\bf InfiniteLoopKeptStates} $m$ value for the infinite algorithm.
 \\item {\\bf FiniteLoops} A series of space-separated numbers. More than one space is allowed. The first
 number is the number of finite algorithm ``movements,'' followed by series of three numbers for
 each ``movement''. Of the three numbers, the first is the number of sites to go forward if positive
 or backward if negative. The second number is the $m$ for this movement and the last number
 is either 0 or 1, 0 will not save state data to disk and 1 will save all data to be able to calculate
 observables. The first movement starts from where the infinite loop left off, at the middle of the
 lattice.
 \\end{itemize}

 \\subsubsection{Example of a Finite loops line in the input file}
 \\begin{verbatim}
 FiniteLoops 4 7 200 0 -7 200 0 7 200 1 7 200 1
 \\end{verbatim}
 The number 4 implies 4 finite loops. The first fine loop is ``7 200 0'', meaning
 go forward 7 steps, use $m=200$ for this finite sweep, and 0: do not store transformation in disk.
 The next is ``-7 200 0'', which goes backwards 7 sites, etc.
 Remember that the finite loops start at the middle of the lattice, where the infinite loop left off.
 \\todo{ADD FIGURE SHOWING WHAT THIS DOES.}

 \\subsubsection{The third number in the triplet}
 The save option is a bitwise option where the
 first bit means  ``save or don't save,'' and the second bit
 ``compute the g.s. or WFT it.''
 So there are 4 combinations (as of today):

 \\begin{tabular}{ll}
 0 & Don't save, compute the g.~s.\\\\
 1 & Save, compute the g.~s.\\\\
 2 & Don't save, WFT the g.~s.\\\\
 3 & Save, WFT the g.~s.\\\\
 \\end{tabular}

 \\subsubsection{Caveats and Troubleshooting}
 \\begin{itemize}
 \\item If \\verb=nofiniteloops= is an option in the options line of the input file then
 the {\\bf FiniteLoops} line in the input file is ignored, and no finite loops are done.
 In this case, MPS++ stops when the infinite algorithm has finished.

 \\item Make sure the first number is the number of triplets that follow.
 \\item Make sure
 you don't fall off the lattice, by going forward or backwards too much.
 Remember that at least one site must remain for the ``system'' part of the lattice.
 So on a 16 site chain, when you start the finite loops you're at the middle, you
 can go forward at most 7 sites, and backwards at most 7 sites.
 \item There is some checking done to the finite loops input, see !PTEX\\_REF{139},
 but you might find that it's not comprehensive.
 \\end{itemize} */
struct FiniteLoop {
	int stepLength; // how much to go right (+) or left (-)
	size_t keptStates; // kept states
	int saveOption; // to save or not to save
	FiniteLoop(int sl,int ks,int so)
		: stepLength(sl),keptStates(ks),saveOption(so)
	{}
};

//!PTEX_LABEL{139}
inline void checkFiniteLoops(const typename ProgramGlobals::Vector<FiniteLoop>::Type& finiteLoop,size_t totalSites)
{
	std::string s = "checkFiniteLoops: I'm falling out of the lattice ";
	std::string loops = "";
	int x = totalSites/2-1; // must be signed
	if (finiteLoop[0].stepLength<0) x++;
	int prevDeltaSign = 1;
	size_t sopt = 0; // have we started saving yet?
	for (size_t i=0;i<finiteLoop.size();i++)  {
		size_t thisSaveOption = (finiteLoop[i].saveOption & 1);
		if (sopt == 1 && thisSaveOption ==0) {
			s = "Error for finite loop number " + ttos(i) + "\n";
			s += "Once you say 1 on a finite loop, then all";
			s += " finite loops that follow must have 1.";
			throw std::runtime_error(s.c_str());
		}
		if (sopt == 0 && thisSaveOption ==1) {
			sopt = 1;
			if (size_t(x) != 1 && size_t(x)!=totalSites-2) {
				s = __FILE__ + std::string(": FATAL: for finite loop number ")
						+ ttos(i) + "\n";
				s += "Saving finite loops must start at the left or";
				s += " right end of the lattice\n";
				throw std::runtime_error(s.c_str());
			}
		}
		// naive location:
		int delta = finiteLoop[i].stepLength;
		x += delta;
		loops = loops + ttos(delta) + " ";

		// take care of bounces:
		if (i>0 && delta*prevDeltaSign < 0) x += prevDeltaSign;
		prevDeltaSign = 1;
		if (delta<0) prevDeltaSign = -1;

		// check that we don't fall out
		bool flag = false;
		if (x<=0) {
			s = s + "on the left end\n";
			flag = true;
		}
		if (size_t(x)>=totalSites-1) {
			s = s + "on the right end\n";
			flag = true;
		}
		if (flag) {
			// complain and die if we fell out:
			s = s + "Loops so far: " + loops + "\n";
			s =s + "x=" + ttos(x) + " last delta=" +
					ttos(delta);
			s =s + " sites=" + ttos(totalSites);
			throw std::runtime_error(s.c_str());
		}
	}

}

std::istream &operator>>(std::istream& is,FiniteLoop& fl)
{
	is>>fl.stepLength;
	is>>fl.keptStates;
	is>>fl.saveOption;
	return is;
}

std::ostream &operator<<(std::ostream& os,const FiniteLoop& fl)
{
	os<<fl.stepLength<<" ";
	os<<fl.keptStates<<" ";
	os<<fl.saveOption;
	return os;
}

struct DmrgCheckPoint {
	bool enabled;
	std::string filename;
};

std::istream &operator>>(std::istream& is,DmrgCheckPoint& c)
{
	is>>c.enabled;
	is>>c.filename;
	return is;
}

//! Structure that contains the Dmrg parameters
/**
 \\inputItem{Model}
 A string indicating the model, be it HubbardOneBand HeisenbergSpinOneHalf, etc.

 \\inputItem{Options}
 A comma-separated list of strings. At least one of the following strings must be provided:
 \\inputSubItem{none}  Use this when no options are given, since the list of strings must be non-null.
 Note that ``none'' does not disable other options.

 \\inputSubItem{hasQuantumNumbers} If this option is given, the program will read the line ``QNS''
 described below and act accordingly. It is recommended that you set this option.  \\\\\n
 \\inputSubItem{wft}  Use the Wave Function Transformation speed-up, which is disabled by default.

 \\inputSubItem{useSu2Symmetry} Use the SU(2) symmetry for the model, and interpret quantum numbers in
 the line ``QNS'' appropriately. The option ``hasQuantumNumbers'' must be set for this to work.

 \\inputSubItem{nofiniteloops}  Don't do finite loops, even if provided under ``FiniteLoops'' below.

 \\inputItem{version}  A mandatory string that is read and ignored. Usually contains the result
 of doing ``git rev-parse HEAD''.

 \\inputItem{outputfile}  The output file. This file will be created if non-existent, and if it
 exits it will be truncated.

 \\inputItem{InfiniteLoopKeptStates}  ``m'' value for the infinite algorithm.

 \\inputItem{FiniteLoops} A series of space-separated numbers. More than one space is allowed.
 The first number is the number of finite algorithm ``movements'', followed by series
 of three numbers for each ``movement''. Of the three numbers, the first
 is the number of sites to go forward if positive or backward if negative.
 The second number is the ``m'' for this ``movement' and the last number is either 0 or 1,
 0 will not save state data to disk and 1 will save all data to be able to calculate observables.
 The first ``movement'' starts from where the infinite loop left off, at the middle of the lattice.
 \\inputItem{QNS}  A space-separated list of numbers. More than one space is allowed.
 The first number is the number of numbers to follow, these numbers being the density of quantum
 numbers for each conserved quantum number to be used.
 In a simpler way, usually this is 3 followed by $n_\\uparrow n_\\downarrow 0$  if not using
 SU(2) symmetry, where  $n_\\uparrow$, and $n_\\downarrow$ are the densities of up and down
 electrons respectively. If there is SU(2) symmetry then this is 3 followed by $n_\\uparrow n_\\downarrow j$,
 where $n_\\uparrow$, and $n_\\downarrow$ are the densities of up and down
 electrons respectively, and $j$ is twice the angular momentum divided by the number of sites.
 */
template<typename RealType_,typename ComplexOrRealType_,typename InputValidatorType>
struct ParametersMpsSolver {

	typedef RealType_ RealType;
	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename ProgramGlobals::Vector<FiniteLoop>::Type FiniteLoopsType;
	typedef typename ProgramGlobals::Vector<RealType>::Type VectorRealType;

	std::string filename;
	size_t keptStatesInfinite;
	FiniteLoopsType finiteLoops;
	std::string version;
	std::string options;
	std::string model;
	VectorRealType targetQuantumNumbers;
	size_t electronsUp,electronsDown;
//	std::string initialMps;
//	RealType tolerance;
//	DmrgCheckPoint checkpoint;
	size_t nthreads;
//	int useReflectionSymmetry;
//	std::string fileForDensityMatrixEigs;

	//! Read Dmrg parameters from inp file
	ParametersMpsSolver(InputValidatorType& io)
	{
		io.readline(model,"Model=");
		io.readline(options,"SolverOptions=");
		io.readline(version,"Version=");
		io.readline(filename,"OutputFile=");
		io.readline(keptStatesInfinite,"InfiniteLoopKeptStates=");
		VectorRealType tmpVec;
		io.read(tmpVec,"FiniteLoops");
		for (size_t i=0;i<tmpVec.size();i+=3) {
			typename ProgramGlobals::Vector<int>::Type xTmp(3);
			for (size_t j=0;j<xTmp.size();j++) xTmp[j]=int(tmpVec[i+j]);
			FiniteLoop fl(xTmp[0],xTmp[1],xTmp[2]);
			finiteLoops.push_back(fl);
		}

		size_t repeat = 0;

		try {
			io.read(repeat,"RepeatFiniteLoopsTimes=");
		}  catch (std::exception& e) {}

		size_t fromFl = 0;
		try {
			io.read(fromFl,"RepeatFiniteLoopsFrom=");
		}  catch (std::exception& e) {}

		size_t upToFl = finiteLoops.size()-1;
		try {
			io.read(upToFl,"RepeatFiniteLoopsTo=");
		}  catch (std::exception& e) {}

		if (upToFl>=finiteLoops.size()) {
			std::string s (__FILE__);
			s += "\nFATAL: RepeatFiniteLoopsTo=" + ttos(upToFl) + " is larger than current finite loops\n";
			s += "\nMaximum is " + ttos(finiteLoops.size())+ "\n";
			throw std::runtime_error(s.c_str());
		}
		if (fromFl>upToFl) {
			std::string s (__FILE__);
			s += "\nFATAL: RepeatFiniteLoopsFrom=" + ttos(fromFl) + " is larger than RepeatFiniteLoopsTo\n";
			s += "\nMaximum is " + ttos(upToFl)+ "\n";
			throw std::runtime_error(s.c_str());
		}
		upToFl++;

		for (size_t i=0;i<repeat;i++) {
			for (size_t j=fromFl;j<upToFl;j++) {
				FiniteLoop fl = finiteLoops[j];
				finiteLoops.push_back(fl);
			}
		}

		if (options.find("hasQuantumNumbers")!=std::string::npos) {
			std::string s = "*** WARNING: hasQuantumNumbers ";
			s += "option is obsolete in input file\n";
			std::cerr<<s;
		}
		try {
			io.read(targetQuantumNumbers,"TargetQuantumNumbers");
		} catch (std::exception& e) {}

		bool hasElectrons = false;
		try {
			io.readline(electronsUp,"TargetElectronsUp");
			io.readline(electronsDown,"TargetElectronsDown");
			hasElectrons = true;
		} catch (std::exception& e) {}

		if (hasElectrons && targetQuantumNumbers.size()>0) {
			std::string s (__FILE__);
			s += "\nFATAL: Specifying both TargetElectronsUp/Down and TargetQuantumNumbers is an error.";
			s += "\nSpecify one or the other only.\n";
			throw std::runtime_error(s.c_str());
		}

		if (!hasElectrons && targetQuantumNumbers.size()==0) {
			std::string s (__FILE__);
			s += "\nFATAL: Either TargetElectronsUp/Down or TargetQuantumNumbers must be specified.\n";
			throw std::runtime_error(s.c_str());
		}

		if (options.find("useSu2Symmetry")!=std::string::npos && hasElectrons) {
			std::string s (__FILE__);
			s += "\nFATAL: TargetElectronsUp/Down cannot be specified while using SU(2) symmetry\n";
			s += "\nTargetQuantumNumbers must be specified instead.\n";
			throw std::runtime_error(s.c_str());
		}

//		io.readline(initialMps,"InitialMps=");
//		tolerance = -1.0;
//		try {
//			io.readline(tolerance,"TruncationTolerance=");
//		} catch (std::exception& e) {}

//		if (options.find("checkpoint")!=std::string::npos)
//			io.readline(checkpoint.filename,"CheckpointFilename=");
//		else if (options.find("restart")!=std::string::npos)
//			io.readline(checkpoint.filename,"RestartFilename=");

		nthreads=1; // provide a default value
		try {
			io.readline(nthreads,"Threads=");
		} catch (std::exception& e) {}

		if (nthreads==0) {
			std::string s (__FILE__);
			s += "\nFATAL: nthreads cannot be zero\n";
			throw std::runtime_error(s.c_str());
		}
	}

};

//! print dmrg parameters
template<typename RealType,typename ComplexOrRealType,typename InputValidatorType>
std::ostream &operator<<(std::ostream &os,
			 ParametersMpsSolver<RealType,ComplexOrRealType,InputValidatorType> const &parameters)
{
	os<<"#This is MPS++\n";
	Provenance provenance;
	os<<provenance;
	os<<"parameters.version="<<parameters.version<<"\n";
	os<<"parameters.model="<<parameters.model<<"\n";
	os<<"parameters.filename="<<parameters.filename<<"\n";
	os<<"parameters.options="<<parameters.options<<"\n";
	os<<"parameters.keptStatesInfinite="<<parameters.keptStatesInfinite<<"\n";
	os<<"finiteLoop\n";
	os<<parameters.finiteLoop;

	if (parameters.targetQuantumNumbers.size()>0) {
		os<<"parameters.targetQuantumNumbers=";
		for (size_t i=0;i<parameters.targetQuantumNumbers.size();i++)
			os<<parameters.targetQuantumNumbers[i]<<" ";
		os<<"\n";
	} else {
		os<<"parameters.electronsUp="<<parameters.electronsUp<<"\n";
		os<<"parameters.electronsDown="<<parameters.electronsDown<<"\n";
	}
	if (parameters.tolerance>0)
		os<<"parameters.tolerance="<<parameters.tolerance<<"\n";
	os<<"parameters.nthreads="<<parameters.nthreads<<"\n";
	os<<"parameters.useReflectionSymmetry="<<parameters.useReflectionSymmetry<<"\n";
	if (parameters.checkpoint.filename!="")
		os<<"parameters.restartFilename="<<parameters.checkpoint.filename<<"\n";
	if (parameters.fileForDensityMatrixEigs!="")
		os<<"parameters.fileForDensityMatrixEigs="<<parameters.fileForDensityMatrixEigs<<"\n";
	return os;
}
} // namespace Dmrg
/*@}*/

#endif
