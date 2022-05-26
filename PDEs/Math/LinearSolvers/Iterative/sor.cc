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
    const double omega,
    const bool verbose) :
  A(A), tolerance(tolerance),
  max_iterations(max_iterations), omega(omega),
  LinearSolverBase(verbose)
{
  Assert(omega > 0 && omega < 2, "Invalid relaxation parameter.");
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
      double factor = omega / *A.diagonal(i);
      double value = (1.0 - omega) * x[i] + b[i] * factor;
      for (const auto el : A.const_row_iterator(i))
        if (el.column != i)
          value -= el.value * x[el.column] * factor;

      //==================== Increment difference
      diff += std::fabs(value - x[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    //==================== Check convergence
    if (diff < tolerance)
    { converged = true; break; }
  }

  if (verbose)
  {
    std::stringstream ss;
    ss << "SOR Solver Status:\n"
       << (converged ? "  CONVERGED\n" : "  NOT CONVERGED\n")
       << (converged ? "  # Iterations: " : "  Difference: ")
       << (converged ? nit : diff) << std::endl;
    std::cout << ss.str();
  }
}
