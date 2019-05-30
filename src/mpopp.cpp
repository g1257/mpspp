#include "InputCheck.h"
#include "Vector.h"
#include "InputNg.h"
#include "ParametersMpsSolver.h"
#include "ModelSelector.h"
#include "Geometry/Geometry.h"
#include "SymmetryLocal.h"
#include "Matrix.h"

typedef double RealType;
typedef std::complex<RealType> ComplexOrRealType;
//typedef RealType ComplexOrRealType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::InputNg<Mpspp::InputCheck> InputNgType;
typedef Mpspp::SymmetryLocal SymmetryLocalType;
typedef PsimagLite::Geometry<RealType,InputNgType::Readable,Mpspp::ProgramGlobals> GeometryType;
typedef Mpspp::ParametersMpsSolver<RealType,ComplexOrRealType,InputNgType::Readable>
ParametersSolverType;
typedef PsimagLite::InputNg<Mpspp::InputCheck>::Readable InputValidatorType;
typedef Mpspp::MpoLocal<ComplexOrRealType, SymmetryLocalType> MpoLocalType;
typedef MpoLocalType::MpoFactorType MpoFactorType;
typedef MpoFactorType::SparseMatrixType SparseMatrixType;
typedef Mpspp::ModelBase<ParametersSolverType,
InputValidatorType,
SymmetryLocalType,
GeometryType> ModelBaseType;
typedef Mpspp::ModelSelector<ModelBaseType> ModelSelectorType;
typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
typedef PsimagLite::Matrix<MatrixType> MatrixMatrixType;
typedef PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
typedef std::pair<SizeType, SizeType> PairSizeType;
typedef PsimagLite::Vector<RealType>::Type VectorRealType;

SizeType pack(SizeType i, SizeType j, SizeType sizeOfI)
{
	assert(i < sizeOfI);
	return i + sizeOfI*j;
}

PairSizeType unpack(SizeType ij, SizeType sizeOfI)
{
	div_t d = div(ij, sizeOfI);
	return PairSizeType(d.rem, d.quot);
}

void toMatrixMatrix(MatrixMatrixType& mm, const MpoFactorType& mpof)
{
	const SizeType rows = mpof.rows();
	const SizeType cols = mpof.cols();
	mm.resize(rows, cols);
	for (SizeType i = 0; i < rows; ++i)
		for (SizeType j = 0; j < rows; ++j)
			mm(i, j) = mpof(i, j).matrix().toDense();
}

SizeType findKeptSize(const VectorRealType& s)
{
	const SizeType n = s.size();
	for (SizeType i = 0; i < n; ++i)
		if (fabs(s[i]) < 1e-10) return i;

	return n;
}

void algoMpoVerticalOnlySvd(const ModelBaseType& model, SizeType maxPower)
{
	const SizeType site = 1;
	const MpoFactorType& wsite = model.hamiltonian()(site);
	MatrixMatrixType w;
	toMatrixMatrix(w, wsite);
	MatrixMatrixType powers = w;
	const SizeType hilbert = w(0, 0).rows();
	const SizeType hilbertSq = hilbert*hilbert;
	PsimagLite::Svd<ComplexOrRealType> svd;
	MatrixType zero(hilbert, hilbert);
	zero.setTo(0.0);
	for (SizeType n = 1; n < maxPower; ++n) {
		std::cout<<n<<" "<<powers.rows()<<"\n";

		// new(wfat, vfat)_{sigma, tau} = \sum_t powers(w1, v1)_{sigma, t} * w(w2, v2)_{t, tau}
		// where wfat = pack(w1, w2)
		// and   vfat = pack(v1, v2)
		// and set matrix1(row1, vfat) = new(wfat, vfat)_{sigma, tau}
		// where row1 = pack(wfat, sigma, tau)
		MatrixType tmp(hilbert, hilbert);
		const SizeType np = powers.rows();
		assert(powers.cols() == np);
		const SizeType nw = w.rows();
		assert(w.cols() == nw);
		MatrixType matrix(hilbert*hilbert*nw*np, nw*np);
		VectorRealType s;
		MatrixType vt;
		for (SizeType w1 = 0; w1 < np; ++w1) {
			for (SizeType v1 = 0; v1 < np; ++v1) {
				const MatrixType& p = powers(w1, v1);
				for (SizeType w2 = 0; w2 < nw; ++w2) {
					const SizeType wfat = pack(w2, w1, nw);
					for (SizeType v2 = 0; v2 < nw; ++v2) {
						const SizeType vfat = pack(v2, v1, nw);
						const MatrixType& ww = w(w2, v2);
						tmp = p*ww;

						for (SizeType sigma = 0; sigma < hilbert; ++sigma) {
							for (SizeType tau = 0; tau < hilbert; ++tau) {
								const SizeType row1 = pack(sigma + tau*hilbert, wfat, hilbertSq);
								matrix(row1, vfat) = tmp(sigma, tau);
							}
						}
					}
				}
			}
		}

		// and then do the SVD matrix1 = U1 S1 V1^\dagger
		vt.clear();
		svd('A', matrix, s, vt);
		const SizeType kept = findKeptSize(s);
		matrix.resize(matrix.rows(), kept);

		// and end up with U1(pack(wfat, sigma, tau), vNoSoFat)
		// where vNoSoFat is truncated to non-zeros of S1
		// now we set matrix2(row2, wfat) = U1(pack(wfat, sigma, tau), vNoSoFat)
		// were row2 = pack(sigma, tau, vNoSoFat)
		const SizeType rows = matrix.rows();
		const SizeType cols = matrix.cols();
		MatrixType matrix2(cols*hilbertSq, rows);
		for (SizeType rowOld = 0; rowOld < rows; ++rowOld) {
			const PairSizeType sigmatauWfat = unpack(rowOld, hilbertSq);
			for (SizeType vNoSoFat = 0; vNoSoFat < cols; ++vNoSoFat) {
				const SizeType row2 = pack(sigmatauWfat.first, vNoSoFat, hilbertSq);
				matrix2(row2, sigmatauWfat.second) += matrix(rowOld, vNoSoFat);
			}
		}

		// and then we SVD  matrix2(row, wfat) = U2 S2 V2^\dagger
		vt.clear();
		svd('A', matrix2, s, vt);
		const SizeType kept2 = findKeptSize(s);
		matrix2.resize(matrix2.rows(), kept2);

		// and end up with U2(pack(sigma, tau, vNoSoFat), wNoSoFat)
		// where wNoSoFat is truncated to non-zeros of S2
		// finally we set powers(wNoSoFat,  vNoSoFat)_{sigma, tau} = U2(row2, wNoSoFat)
		// were row2 = pack(sigma, tau, vNoSoFat)
		const SizeType rows2 = matrix2.rows();
		const SizeType cols2 = matrix2.cols();
		const SizeType colsMin = std::min(cols, cols2);

		powers.clear();
		powers.resize(colsMin, colsMin);
		powers.setTo(zero);
		for (SizeType row2 = 0; row2 < rows2; ++row2) {
			const PairSizeType sigmatauVnoSoFat = unpack(row2, hilbertSq);
			div_t sigmatau = div(sigmatauVnoSoFat.first, hilbert);
			if (sigmatauVnoSoFat.second >= colsMin) continue;
			for (SizeType wNoSoFat = 0; wNoSoFat < colsMin; ++wNoSoFat) {
				powers(wNoSoFat, sigmatauVnoSoFat.second)(sigmatau.rem, sigmatau.quot)
				        += matrix2(row2, wNoSoFat);
			}
		}
	}
}

int main(int argc, char** argv)
{
	Mpspp::InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename";
	int precision = 6;
	bool versionOnly = false;
	SizeType maxPower = 2;
	while ((opt = getopt(argc, argv,"f:p:P:V")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'P':
			maxPower = atoi(optarg);
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename=="" && !versionOnly) {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	//! setup distributed parallelization
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);

	// print license
	if (ConcurrencyType::root()) {
		std::cout<<Mpspp::ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryType geometry(io);

	ParametersSolverType mpsSolverParams(io);

	ModelSelectorType modelSelector(mpsSolverParams.model);

	const ModelBaseType& model = modelSelector(mpsSolverParams,io,geometry);

	algoMpoVerticalOnlySvd(model, maxPower);
}
