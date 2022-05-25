#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "sparse_matrix.h"
#include "vector.h"

namespace pdes::Math
{

class GaussSeidelSolver
{
private:
  double tol = 1.0e-8;
  size_t max_iter = 1000;

public:
  /** Default constructor. */
  GaussSeidelSolver() = default;

  /** Constructor with specified iteration controls. */
  GaussSeidelSolver(const double tolerance,
                    const size_t max_iterations);

  /** Solve the system using the Gauss Seidel iterative method. */
  void
  solve(const SparseMatrix& A, const Vector& b, Vector& x);
};

}


#endif //GAUSS_SEIDEL_H
