#include "gauss_seidel.h"


using namespace pdes::Math;

LinearSolver::GaussSeidel::
GaussSeidel(const SparseMatrix& A,
            const double tolerance,
            const size_t max_iterations,
            const bool verbose) :
  SOR(A, 1.0, tolerance, max_iterations, verbose)
{}
