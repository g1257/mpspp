/*
Copyright (c) 2012, UT-Battelle, LLC
All rights reserved

[Mpspp++, Version 0.1]
[by K. Al-Hassanieh, Oak Ridge National Laboratory]
[by J. Rincon, Oak Ridge National Laboratory]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."
 
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
/** \ingroup Mpspp */
/*@{*/

/*! \file VectorWithOffset.h
 *
 *  A class to represent a vector like this 000000 XXXXXXXX 0000000000000
 *  offset_ is where the first X (non-zero element) is.
 *  data_ contains the nonzero part.
 *  sizE_ is the size of the vector
 */
#ifndef VECTOR_WITH_OFFSET_H
#define VECTOR_WITH_OFFSET_H
#include "Vector.h"

namespace Mpspp {
	template<typename FieldType>
	struct VectorWithOffset {
	public:
		typedef FieldType value_type;
		typedef std::pair<size_t,size_t> PairType;
		typedef std::vector<FieldType> VectorType;

		static const FieldType zero_;

		VectorWithOffset()  : size_(0),offset_(0),m_(0) { }

		template<typename SomeBasisType>
		VectorWithOffset(const std::vector<size_t>& weights,
		                 const SomeBasisType& someBasis)
		: size_(someBasis.size())
		{
			bool found = false;
			for (size_t i=0;i<weights.size();i++) {
				if (weights[i]>0) {
					if (found)
						throw std::runtime_error("VectorWithOffset::"
						  " more than one non-zero sector found. "
						  " Maybe you should be using VectorWithOffsets instead?\n");
					data_.resize(weights[i]);
					offset_ = someBasis.partition(i);
					m_ = i;
					found = true;
				}
			}
		}

		void resize(size_t x)
		{
			size_ = x;
			data_.clear();
			offset_=0;
			m_=0;
		}

		template<typename SomeBasisType>
		void set(const std::vector<VectorType>& v,
		         const SomeBasisType& someBasis)
		{
			bool found = false;
			size_ = someBasis.size();
			for (size_t i=0;i<v.size();i++) {
				if (v[i].size()>0) {
					if (found) throw std::runtime_error("VectorWithOffset::"
						" more than one non-zero sector found. "
						" Maybe you should be using VectorWithOffsets instead?\n");
					data_ = v[i];
					offset_ = someBasis.partition(i);
					m_ = i;
					found = true;
				}
			}
			if (!found) throw std::runtime_error("Set failed\n");
		}
		
		template<typename SomeBasisType>
		void fromFull(const VectorType& v,const SomeBasisType& someBasis)
		{
			size_ = someBasis.size();
			try {
				m_ = findPartition(v,someBasis);
				offset_ = someBasis.partition(m_);
				size_t total = someBasis.partition(m_+1) - offset_;
				data_.resize(total);
				for (size_t i=0;i<total;i++) data_[i] = v[i+offset_];
			} catch (std::exception& e) {
				std::cerr<<e.what();
				m_=0;
				offset_=0;
				data_.resize(0);
			}
		}

		size_t sectors() const { return 1; }

		size_t sector(size_t dummy) const { return m_; }

		size_t offset(size_t dummy) const { return offset_; }

		size_t effectiveSize(size_t dummy) const { return data_.size(); }

		void setDataInSector(const VectorType& v,size_t dummy)
		{
			data_=v;
		}

		void extract(VectorType& v, size_t dummy = 0) const
		{
			v=data_;
		}

		template<typename SparseVectorType>
		void toSparse(SparseVectorType& sv) const
		{
			sv.resize(size_);
			for (size_t i=0;i<data_.size();i++)
				sv[i+offset_] = data_[i];
		}

		template<typename IoOutputter>
		void save(IoOutputter& io,const std::string& label) const
		{
			io.print(label);
			io.print("#size=",size_);
			io.print("#offset=",offset_);
			io.print("#m=",m_);
			io.printVector(data_,"#data");
		}

		template<typename IoInputter>
		void load(IoInputter& io,const std::string& label,size_t counter=0)
		{
			io.advance(label,counter);
			int x = 0;
			io.readline(x,"#size=");
			if (x<0) throw std::runtime_error("VectorWithOffset::load(...): size<0\n");
			size_ = x;
			io.readline(x,"#offset=");
			if (x<0) throw std::runtime_error("VectorWithOffset::load(...): offset<0\n");
			offset_ = x;
			io.readline(x,"#m=");
			if (x<0) throw std::runtime_error("VectorWithOffset::load(...): m<0\n");
			m_ = x;
			io.read(data_,"#data");
		}

		size_t size() const { return size_; }

		size_t effectiveSize() const { return data_.size(); }

		size_t offset() const { return offset_; }

		 const FieldType& operator[](size_t i) const
		{
			if (i<offset_ || i>= (offset_+data_.size())) return zero_;
			//assert(i>=offset_ && i<offset_+data_.size());
			return data_[i-offset_];
		}

		FieldType& operator[](size_t i)
		{
// 			if (i<offset_ || i>= (offset_+data_.size()))
// 				throw std::runtime_error("VectorWithOffset\n");
			return data_[i-offset_];
		}
		
		const FieldType& fastAccess(size_t i,size_t j) const 
		{
			return data_[j];
		}

		template<typename FieldType2>
		friend FieldType2 std::norm(const Mpspp::VectorWithOffset<FieldType2>& v);

		template<typename FieldType2>
		friend FieldType2 std::norm(const Mpspp::VectorWithOffset<std::complex<FieldType2> >& v);

		template<typename FieldType2>
		friend FieldType2 operator*(const Mpspp::VectorWithOffset<FieldType2>& v1,
									const Mpspp::VectorWithOffset<FieldType2>& v2);

		template<typename FieldType3,typename FieldType2>
		friend VectorWithOffset<FieldType2> operator*(const FieldType3& value,const VectorWithOffset<FieldType2>& v);
		
		template<typename FieldType2>
		friend FieldType2 multiply(const VectorWithOffset<FieldType2>& v1,
		                                         const VectorWithOffset<FieldType2>& v2);

	private:
		template<typename SomeBasisType>
		size_t findPartition(const VectorType& v,const SomeBasisType& someBasis)
		{
			bool found = false;
			size_t p = 0;
			for (size_t i=0;i<someBasis.partition()-1;i++) {
				if (nonZeroPartition(v,someBasis,i)) {
					if (found) throw std::runtime_error("VectorWithOFfset::"
						" More than one partition found\n");
					found = true;
					p = i;
				}
			}
			if (!found)
				throw std::runtime_error("VectorWithOFfset::"
				                         " No partition found\n");
			return p;
		}

		template<typename SomeBasisType>
		bool nonZeroPartition(const VectorType& v,const SomeBasisType& someBasis,size_t i)
		{
			typename VectorType::value_type zero = 0;
			for (size_t j=someBasis.partition(i);j<someBasis.partition(i+1);j++) {
				if (v[j]!=zero) return true;
			}
			return false;
		}

		size_t size_;
		VectorType data_;
		size_t offset_;
		size_t m_; // partition 
	}; // class VectorWithOffset

	template<typename FieldType>
	const FieldType VectorWithOffset<FieldType>::zero_=0;

//	template<typename FieldType>
//	std::ostream& operator<<(std::ostream& os,const VectorWithOffset<FieldType>& s)
//	{
//		s.print(os,"VectorWithOffset");
//		return os;
//	}
	
	template<typename FieldType>
	inline FieldType operator*(const Mpspp::VectorWithOffset<FieldType>& v1,
							   const Mpspp::VectorWithOffset<FieldType>& v2)
	{
		return (v1.data_ * v2.data_);
	}

	template<typename FieldType,typename FieldType2>
	inline VectorWithOffset<FieldType2> operator*(const FieldType& value,const VectorWithOffset<FieldType2>& v)
	{
		VectorWithOffset<FieldType2> w = v;
		w.data_ *= value;
		return w;
	}
	
	template<typename FieldType>
	inline FieldType multiply(const VectorWithOffset<FieldType>& v1,
	                                        const VectorWithOffset<FieldType>& v2)
	{
		return v1.data_*v2.data_; // call to * will conj()
	}

} // namespace Mpspp

namespace std {
	template<typename FieldType>
	inline FieldType norm(const Mpspp::VectorWithOffset<FieldType>& v)
	{
		return PsimagLite::norm(v.data_);
	}

	template<typename FieldType>
	inline FieldType norm(const Mpspp::VectorWithOffset<std::complex<FieldType> >& v)
	{
		return PsimagLite::norm(v.data_);
	}

	template<typename FieldType>
	inline std::complex<FieldType> operator*(const Mpspp::VectorWithOffset<std::complex<FieldType> >& v1,
											 const Mpspp::VectorWithOffset<std::complex<FieldType> >& v2)
	{
		throw std::runtime_error("Unimplemented\n");
	}

}
/*@}*/
#endif
