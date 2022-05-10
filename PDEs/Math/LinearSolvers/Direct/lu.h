#ifndef LU_H
#define LU_H

#include "LinearSolvers/linear_solver.h"

namespace math
{

/** A class for an LU decomposition solver. */
template<typename value_type>
class LU : public LinearSolver<value_type>
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
 LU(Matrix<value_type>& matrix, const bool pivot_option = true)
   : LinearSolver<value_type>(matrix), pivot(pivot_option)
 {}

public:
 void set_pivot_option(const bool pivot_option) { pivot = pivot_option; }
 bool get_pivot_option() const { return pivot; }

  void setup() override;
  Vector<value_type> solve(const Vector<value_type>& b) override;
};

}
#endif //LU_H
