#include "sor.h"

#include "vector.h"
#include "sparse_matrix.h"

#include "macros.h"

#include <cmath>


using namespace pdes::Math;

LinearSolver::SOR::
SOR(const SparseMatrix& A,
    const double tolerance,
    const size_t max_iterations,
    const double omega) :
  A(A), tolerance(tolerance),
  max_iterations(max_iterations), omega(omega)
{
  Assert(omega >= 0 && omega <= 1, "Invalid relaxation parameter.");
  Assert(A.n_rows() == A.n_cols(), "Square matrix required.");
  Assert(tolerance > 0.0, "Illegal negative tolerance specified.");
}


void LinearSolver::SOR::
solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  Assert(b.size() == n, "Dimension mismatch error.");
  Assert(x.size() == n, "Dimension mismatrch error.");

  double diff;
  size_t nit;
  bool converged = false;

  //======================================== Iteration loop
  for (nit = 0; nit < max_iterations; ++nit)
  {
    diff = 0.0;
    for (size_t i = 0; i < A.n_rows(); ++i)
    {
      //==================== Compute element-wise update
      double a_ii = *A.diagonal(i);
      double value = (1.0 - omega)*x[i]*a_ii + omega*b[i];
      for (const auto el : A.const_row_iterator(i))
        if (el.column != i)
          value -= el.value * x[el.column];
      value /= a_ii;

      //==================== Increment difference
      diff += std::fabs(value - x[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    //==================== Check convergence
    if (diff < tolerance)
    { converged = true; break; }
  }

  std::stringstream ss;
  ss << "SOR Solver Status: "
     << (converged? "CONVERGED,  " : "NOT CONVERGED,  ")
     << (converged? "# Iterations: " : "Difference: ")
     << (converged? nit : diff) << std::endl;
  std::cout << ss.str();
}
