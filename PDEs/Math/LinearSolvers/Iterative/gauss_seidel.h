#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "LinearSolvers/linear_solver.h"
#include "LinearSolvers/Iterative/sor.h"

#include <cstddef>

namespace pdes::Math::LinearSolver
{

/**
 * Implementation of the Gauss Seidel iterative method.
 *
 * \note This method corresponds to the SOR method where the relaxation
 *       parameter \f$ \omega = 1.0 \f$. For this reason, implementations
 *       are borrowed from the \ref SOR class.
 */
class GaussSeidel : public SOR
{
public:

  /**
   * Default constructor.
   */
  GaussSeidel(const SparseMatrix& A,
              const double tolerance = 1.0e-8,
              const size_t max_iterations = 1000,
              const bool verbose = false);
};

}


#endif //GAUSS_SEIDEL_H
