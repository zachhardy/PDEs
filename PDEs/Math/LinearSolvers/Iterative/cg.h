#ifndef CG_H
#define CG_H

#include "sparse_matrix.h"
#include "vector.h"
#include "linear_solver.h"


namespace pdes::Math
{


/**
 * A class which implements a conjugate gradient solver.
 */
class ConjugateGradientSolver
{
public:
  using value_type = SparseMatrix::value_type;

private:
  const SparseMatrix& A;
  val

};


}

#endif //CG_H
