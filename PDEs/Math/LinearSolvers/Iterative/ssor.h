#ifndef SSOR_H
#define SSOR_H

#include "LinearSolvers/Iterative/sor.h"


namespace pdes::Math::LinearSolver
{


/**
 * Implementation of the symmetric successive over relaxation (SSOR)
 * iterative method.
 */
class SSOR : public SOR
{
public:
  /**
   * Default constructor.
   */
  SSOR(const SparseMatrix& A,
       const double omega = 1.5,
       const double tolerance = 1.0e-8,
       const size_t max_iteration = 1000,
       const bool verbose = false);

  /**
   * Solve the system using the SSOR iterative method.
   */
  void
  solve(Vector& x, const Vector& b) const override;
};

}

#endif //SSOR_H
