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

/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "TypeToString.h"
#include "Options.h"

namespace Mpspp {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	InputCheck() : optsReadable_(0) {}

	~InputCheck()
	{
		if (optsReadable_!=0) delete optsReadable_;
	}

	PsimagLite::String import() const
	{
		return "";
	}

	bool checkSimpleLabel(const PsimagLite::String& label,
	                      SizeType line) const
	{
		return true;
	}


	bool check(const PsimagLite::String& label,const VectorStringType& vec,size_t line) const
	{
		if (label=="JMVALUES") {
			if (vec.size()!=2) return error1("JMVALUES",line);
			return true;
		} else if (label=="RAW_MATRIX") {
			size_t row = atoi(vec[0].c_str());
			size_t col = atoi(vec[1].c_str());
			size_t n = row*col;
			if (vec.size()!=n+2) return error1("RAW_MATRIX",line);
			return true;
		} else if (label=="Connectors") {
			return true;
		} else if (label=="MagneticField") {
			return true;
		} else if (label=="FiniteLoops") {
			size_t n = atoi(vec[0].c_str());
			if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
			return true;
		}
		return false;
	}

	void check(const PsimagLite::String& label,const PsimagLite::String& val,size_t line)
	{
		if (label!="SolverOptions") return;
		VectorStringType registerOpts;

		//			registerOpts.push_back("restart");
		registerOpts.push_back("nolocalsymm");
		registerOpts.push_back("debugmatrix");
		registerOpts.push_back("test");
		registerOpts.push_back("notruncation");
		registerOpts.push_back("InternalProductStored");

		PsimagLite::Options::Writeable optWriteable(registerOpts,PsimagLite::Options::Writeable::PERMISSIVE);
		optsReadable_ = new  OptionsReadableType(optWriteable,val);
	}

	bool isSet(const PsimagLite::String& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	//		void checkForThreads(size_t nthreads) const
	//		{
	//			if (nthreads==1) return;

	//			PsimagLite::String message1(__FILE__);
	//			message1 += " FATAL: You are requesting nthreads>0 but you did not compile with USE_PTHREADS enabled\n";
	//			message1 += " Either set Threads=1 in the input file (you won't have threads though) or\n";
	//			message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile and recompile\n";
	//			throw std::runtime_error(message1.c_str());
	//		}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

private:

	bool error1(const PsimagLite::String& message,size_t line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw std::runtime_error(s.c_str());

	}

	OptionsReadableType* optsReadable_;

}; // class InputCheck
} // namespace Mpspp

/*@}*/
#endif

