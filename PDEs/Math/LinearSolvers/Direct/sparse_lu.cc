#include "sparse_lu.h"


//template<typename number>
//void math::SparseLU<number>::setup()
//{
//  auto& A = this->A;
//  uint64_t n = A.n_rows();
//
//  // Initialize the pivot mappings such that each row maps to itself
//  for (uint64_t i = 0; i < n; ++i)
//    row_pivots.emplace_back(i);
//
//  //======================================== Apply Doolittle algorithm
//  for (uint64_t i = 0; i < n; ++i)
//  {
//    number* a_ii = A.locate(i, i);
//    Assert(a_ii != nullptr && *a_ii != 0.0,
//           "Diagonal terms must be nonzero.");
//
//    if (pivot)
//    {
//      /* Find the row containing the largest magnitude entry for column i.
//       * This is only done for the sub-diagonal elements. */
//      uint64_t argmax = i;
//      number max = std::fabs(*a_ii);
//      for (uint64_t k = i + 1; k < n; ++k)
//      {
//        number* a_ki = A.locate(k, i);
//        if (a_ki != nullptr && *a_ki > max)
//        {
//          argmax = k;
//          max = std::fabs(*a_ki);
//        }
//      }
//
//      // If the sub-diagonal is uniformly zero, throw error
//      Assert(max != 0.0, "Singular matrix error.");
//
//      /* Swap the current row and the row containing the largest magnitude
//       * entry corresponding for the current column. This is done to improve
//       * the numerical stability of the algorithm. */
//      if (argmax != i)
//      {
//        std::swap(row_pivots[i], row_pivots[argmax]);
//        A.swap_row(i, argmax);
//      }
//    }//if pivot
//
//    // Compute the elements of the LU decomposition
//    for (uint64_t j = i + 1; j < n; ++j)
//    {
//      number* a_ji = A.locate(j, i);
//      if (a_ji != nullptr && *a_ji != 0.0)
//      {
//        /* Lower triangular components. This represents the row operations
//         * performed to attain the upper-triangular, row-echelon matrix. */
//        A.set(j, i, (*a_ji) / (*a_ii));
//
//        /* Upper triangular components. This represents the row-echelon form of
//         * the original matrix. */
//        for (auto entry = A.begin(i); entry != A.end(i); ++entry)
//        {
//          uint64_t k = entry.column();
//          if (k > i)
//            A.add(j, k, -(*a_ji) * entry.value());
//        }
//      }
//    }
//  }
//  this->initialized = true;
//}
//
//
//template<typename number>
//math::Vector<number>
//math::SparseLU<number>::solve(const Vector<number>& b)
//{
//  auto& A = this->A;
//  Assert(b.size() == A.n_rows(), "Dimension mismatch error.");
//  if (not this->initialized)
//    this->setup();
//
//  uint64_t n = A.n_rows();
//  number value = 0.0;
//
//  // Forward solve
//  Vector x(n, 0.0);
//  for (uint64_t i = 0; i < n; ++i)
//  {
//    value = b[row_pivots[i]];
//    for (auto it = A.begin(i); it != A.end(i); ++it)
//    {
//      const uint64_t j = it.column();
//      value -= (j < i)? it.value() * x[j] : 0.0;
//    }
//    x[i] = value;
//  }
//
//  // Backward solve
//  for (int i = n - 1; i >= 0; --i)
//  {
//    value = x[i];
//    for (auto it = A.begin(i); it != A.end(i); ++it)
//    {
//      const uint64_t j = it.column();
//      value -= (j > i)? it.value() * x[j] : 0.0;
//    }
//    x[i] = value / A(i, i);
//  }
//  return x;
//}
//
//
//template class math::SparseLU<double>;
