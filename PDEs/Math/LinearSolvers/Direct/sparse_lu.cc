#include "sparse_lu.h"


template<typename value_type>
void math::SparseLU<value_type>::setup()
{
  auto& A = this->A;
  uint64_t n = A.n_rows();

  // Initialize the pivot mappings such that each row maps to itself
  for (uint64_t i = 0; i < n; ++i)
    row_pivots.emplace_back(i);

  //======================================== Apply Doolittle algorithm
  for (uint64_t i = 0; i < n; ++i)
  {
    Assert(A.exists(i, i) && A(i, i) != 0.0,
           "Diagonal terms must be nonzero.");

    if (pivot)
    {
      /* Find the row containing the largest magnitude entry for column i.
       * This is only done for the sub-diagonal elements. */
      uint64_t argmax = i;
      value_type max = std::fabs(A(i, i));
      for (uint64_t k = i + 1; k < n; ++k)
      {
        if (A.exists(k, i) && A(k, i) > max)
        {
          argmax = k;
          max = std::fabs(A(k, i));
        }
      }

      // If the sub-diagonal is uniformly zero, throw error
      Assert(max != 0.0, "Singular matrix error.");

      /* Swap the current row and the row containing the largest magnitude
       * entry corresponding for the current column. This is done to improve
       * the numerical stability of the algorithm. */
      if (argmax != i)
      {
        std::swap(row_pivots[i], row_pivots[argmax]);
        A.swap_row(i, argmax);
      }
    }//if pivot

    // Compute the elements of the LU decomposition
    for (uint64_t j = i + 1; j < n; ++j)
    {
      if (A.exists(j, i) && A(j, i) != 0.0)
      {
        /* Lower triangular components. This represents the row operations
         * performed to attain the upper-triangular, row-echelon matrix. */
        A.insert(j, i, A(j, i) / A(i, i), false);

        /* Upper triangular components. This represents the row-echelon form of
         * the original matrix. */
        for (uint64_t k = i + 1; k < n; ++k)
          if (A.exists(i, k))
            A.insert(j, k, -A(j, 1) * A(i, k));
      }
    }
  }
  this->initialized = true;
}


template<typename value_type>
math::Vector<value_type>
math::SparseLU<value_type>::solve(const Vector<value_type>& b)
{
  auto& A = this->A;
  Assert(b.size() == A.n_rows(), "Dimension mismatch error.");
  if (not this->initialized)
    this->setup();

  uint64_t n = A.n_rows();
  value_type value = 0.0;

  // Forward solve
  Vector x(n, 0.0);
  for (uint64_t i = 0; i < n; ++i)
  {
    value = b[row_pivots[i]];
    for (auto it = A.begin(i); it != A.end(i); ++it)
    {
      const uint64_t j = it.column();
      value -= (j < i)? it.value() * x[j] : 0.0;
    }
    x[i] = value;
  }

  // Backward solve
  for (int i = n - 1; i >= 0; --i)
  {
    value = x[i];
    for (auto it = A.begin(i); it != A.end(i); ++it)
    {
      const uint64_t j = it.column();
      value -= (j > i)? it.value() * x[j] : 0.0;
    }
    x[i] = value / A(i, i);
  }
  return x;
}


template class math::SparseLU<double>;
