#ifndef SSOR_H
#define SSOR_H

#include "LinearSolvers/Iterative/sor.h"


namespace Math::LinearSolver
{


  /**
   * Implementation of the symmetric successive over relaxation (SSOR)
   * iterative method.
   */
  class SSOR : public SOR
  {
  public:
    /** Default constructor. */
    SSOR(const double omega = 1.5,
         const Options& opts = Options());

    /** Solve the system using the SSOR iterative method. */
    void solve(Vector& x, const Vector& b) const override;
  };

}

#endif //SSOR_H
