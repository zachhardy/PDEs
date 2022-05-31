#include "gauss_seidel.h"


using namespace Math;

LinearSolver::GaussSeidel::
GaussSeidel(const SparseMatrix& A, const Options& opts) :
  SOR(A, opts, "Gauss-Seidel")
{ omega = 1.0; }
