// [[Rcpp::depends("RcppEigen")]]
#include "RcppEigen.h"
#include <time.h>
#include <R_ext/BLAS.h>

#define BILLION 1000000000L

// #include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

class MatrixReplacement;
using Eigen::MatrixXd;

namespace Eigen {
namespace internal {
	// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
	template<>
	struct traits<MatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double> >
	{};
}
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
	// Required typedefs, constants, and method:
	typedef double Scalar;
	typedef double RealScalar;
	typedef int StorageIndex;
	enum {
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor = false
	};

	Index rows() const { return m_M->rows() * m_M->rows(); }
	Index cols() const { return m_M->cols() * m_M->rows(); }
	template<typename Rhs>
	Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
		return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
	}

	// Custom API:
	MatrixReplacement() : m_M(0) {}
	void attachMatrix(const Eigen::MatrixXd &M) { m_M = &M; }
	const Eigen::MatrixXd& getM() const { return *m_M; }

private:
	const Eigen::MatrixXd *m_M;
};

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
	template<typename Rhs>
	struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct>
	: generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> >
	{
		typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;
		template<typename Dest>
		static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
		{
			// This method should implement "dst += alpha * lhs * rhs" inplace,
			// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
			assert(alpha == Scalar(1) && "scaling is not implemented");

			time_t now;
			// time(&now);
			// Rprintf("%s - Computing matrix-free Ax\n", Rcpp::Datetime(now).format().c_str());

			struct timespec start, end;
			clock_gettime(CLOCK_MONOTONIC, &start);

			const Eigen::MatrixXd& M = lhs.getM();
			size_t m = M.rows();
			Eigen::VectorXd x_j(m);
			Eigen::VectorXd v(m);

			double* v_ptr = v.data();
			double* x_j_ptr = x_j.data();
			double* dst_ptr = dst.data();
			const double* rhs_ptr = rhs.data();
			const double* M_ptr = M.data();

			// For call to BLAS
			double zero = 0.0;
			double one = 1.0;
    		int ione = 1;
    		int int_m = int(m);

			dst += rhs;
			for (size_t j = 0; j < m; ++j) {
				for (size_t l = 0; l < m; ++l) {
					x_j_ptr[l] = rhs_ptr[l + j*m];
				}

				// Directly calling BLAS seems to be significantly faster than using
				// straightfoward multiplication operation in Eigen.
				// v = M * x_j;
    			F77_CALL(dgemv)("n", &int_m, &int_m, &one, M_ptr, &int_m, x_j_ptr, &ione, &zero, v_ptr, &ione);

				const Eigen::VectorXd& w = M.col(j);
				const double* w_ptr = w.data();

				for (size_t i = 0; i < m; ++i) {
					double M_ij = w_ptr[i];
					for (size_t l = 0; l < m; ++l) {
						// Using C arrays directly seems to perform significantly faster than
						// indexing into Eigen arrays via element accessors.
						// dst(l + i*m) -= v(l) * M_ij;
						dst_ptr[l + i*m] -= v_ptr[l] * M_ij;
					}
				}
			}

			time(&now);
			clock_gettime(CLOCK_MONOTONIC, &end);
			double diff = double(end.tv_sec - start.tv_sec) + double(end.tv_nsec - start.tv_nsec) / BILLION;
			Rprintf("%s - Ax rep took %f sec\n", Rcpp::Datetime(now).format().c_str(), diff);
		}
	};
}
}

// [[Rcpp::export]]
Rcpp::List solve_Gamma0(const Eigen::MatrixXd& M, const Eigen::MatrixXd& Sigma)
{
	size_t m = Sigma.rows();

	MatrixReplacement A;
	A.attachMatrix(M);

	Eigen::VectorXd vecSigma(m*m);
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++) {
			vecSigma(i + m*j) = Sigma(i,j);
		}
	}

	Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
	bicg.compute(A);
	Eigen::VectorXd x = bicg.solve(vecSigma);

	return Rcpp::List::create(
		Rcpp::Named("x") = x,
		Rcpp::Named("iterations") = bicg.iterations(),
		Rcpp::Named("error") = bicg.error()
	);
}
