#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H


namespace PDEs
{
  namespace Math
  {
    // forward declarations
    class Vector;


    namespace LinearSolvers
    {

      /**
       * Solve a system using Gaussian elimination.
       *
       * This routine factors the matrix and right-hand side into row-echelon
       * form, then uses back substitution to directly obtain the solution.
       * Row-echelon form creates an upper triangular system. This is done by
       * traversing through each column and performing row-wise operations to
       * eliminate all sub-diagonal coefficients. For numerical stability,
       * pivoting can be used to ensure that for each subsequent column,
       * division by the largest element is carried out.
       *
       * \param A An \f$ n \times n \f$ matrix.
       * \param b A vector of length \f$ n \f$.
       * \param pivot A flag for whether pivoting is performed.
       * \return The solution \f$ \vec{x} \f$ of
       *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
       */
      template<class MatrixType>
      Vector
      gaussian_elimination(MatrixType& A, Vector& b, const bool pivot = true);
    }
  }
}
#endif //GAUSSIAN_ELIMINATION_H
