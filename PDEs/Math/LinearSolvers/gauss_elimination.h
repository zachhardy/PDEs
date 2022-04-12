#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H

#include "linear_solver.h"


//######################################################################
class GaussElimination : public LinearSolver
{
public:
  bool with_pivoting = false;
  bool with_normalization = false;

public:
  /// Default constructor with a matrix and right-hand side vector.
  GaussElimination(Matrix& matrix, Vector& rhs)
    : LinearSolver(matrix, rhs)
  {}

  /// Default constructor plus options.
  GaussElimination(Matrix& matrix, Vector& rhs, bool pivoting, bool normalize)
    : LinearSolver(matrix, rhs),
      with_pivoting(pivoting),
      with_normalization(normalize)
  {}


  void setup() override
  {
    if (not initailized)
    {
      std::stringstream err;
      err << "GaussElimination::" << __FUNCTION__ << ": "
          << "No matrix available to compute the row-echelon form.";
      throw std::runtime_error(err.str());
    }
    row_echelon();
  }


  Vector solve() override { return backward_substitution_solve(); }

private:
  void row_echelon();
  Vector backward_substitution_solve();




};


#endif //GAUSS_ELIMINATION_H
