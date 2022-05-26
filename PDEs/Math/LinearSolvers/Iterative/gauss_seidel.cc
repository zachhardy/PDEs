#include "gauss_seidel.h"

#include "vector.h"
#include "sparse_matrix.h"

#include "macros.h"

#include <cmath>


using namespace pdes::Math;

LinearSolver::GaussSeidel::
GaussSeidel(const SparseMatrix& A,
            const double tolerance,
            const size_t max_iterations,
            const bool verbose) :
  SOR(A, tolerance, max_iterations, 1.0, verbose)
{
  Assert(A.n_rows() == A.n_cols(), "Square matrix required.");
  Assert(tolerance > 0.0, "Illegal negative tolerance specified.");
}
