#ifndef SOR_H
#define SOR_H

#include "linear_solver.h"
#include "sparse_matrix.h"
#include "vector.h"


namespace pdes::Math
{

/**
 * Implementation of the SOR iterative method.
 */
class SORSolver
{
public:
  using value_type = SparseMatrix::value_type;

private:
  const SparseMatrix& A;
  value_type  omega;

  value_type tol;
  size_t maxiter;

public:

};


}


#endif //SOR_H
