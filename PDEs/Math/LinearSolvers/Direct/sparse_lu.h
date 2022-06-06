#ifndef SPARSE_LU_H
#define SPARSE_LU_H

#include "../linear_solver.h"
#include "Sparse/sparse_matrix.h"

#include <vector>
#include <cstddef>


namespace Math::LinearSolver
{

  /**
   * Implementation of a sparse LU solver.
   */
  class SparseLU : public DirectSolverBase<SparseMatrix>
  {
  private:
    bool pivot_flag = true;

    /**
     * The pivot mapping vector.
     * The index corresponds to the initial row number and the value to the
     * pivoted row number. This is used to map the right-hand side vector to the
     * correct row when solving.
     */
    std::vector<size_t> row_pivots;

  public:
    /**
     * Default constructor.
     */
    SparseLU(const bool pivot = true);

    /**
     * Attach a matrix to the solver.
     */
    void
    set_matrix(const SparseMatrix& matrix) override;

    /**
     * Factor the matrix \f$ \boldsymbol{A} \f$ in-place.
     *
     * \see LU::factorize
     */
    void
    factorize() override;

    /**
     * Solve the LU factored linear system.
     *
     * \see LU::solve
     */
    void
    solve(Vector& x, const Vector& b) const override;

    /**
     * Set the pivot option.
     */
    void
    pivot(const bool flag);

    /**
     * Return the pivot option.
     */
    bool
    pivot() const;


    using LinearSolverBase::solve;
  };

}

#endif //SPARSE_LU_H
