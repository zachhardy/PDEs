#ifndef LU_H
#define LU_H

#include "LinearSolvers/linear_solver.h"

namespace math
{

/** A class for an LU decomposition solver. */
template<typename number>
class LU : public LinearSolver<number>
{
private:
 bool pivot = true;

 /**
  * The pivot mapping vector.
  * The index corresponds to the initial row number and the value to the
  * pivoted row number. This is used to map the right-hand side vector to the
  * correct row when solving.
  */
 std::vector<size_t> row_pivots;

public:
 LU(Matrix<number>& matrix, const bool pivot_option = true)
   : LinearSolver<number>(matrix), pivot(pivot_option)
 {}

public:
 void set_pivot_option(const bool pivot_option) { pivot = pivot_option; }
 bool get_pivot_option() const { return pivot; }

  void setup() override;
  Vector<number> solve(const Vector<number>& b) override;
};

}
#endif //LU_H
