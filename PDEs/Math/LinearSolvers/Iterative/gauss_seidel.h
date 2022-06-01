#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "LinearSolvers/Iterative/sor.h"

#include <cstddef>


namespace Math::LinearSolver
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
                const Options& opts = Options());
  };

}


#endif //GAUSS_SEIDEL_H
