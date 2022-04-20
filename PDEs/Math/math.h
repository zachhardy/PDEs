#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <vector>

namespace math
{

// Forward declarations
class Vector;
class Matrix;

// Factorizations
void row_echelon_form(Matrix& A, Vector& b, const bool pivot = true);
std::vector<size_t> lu_factorization(Matrix& A, const bool pivot = true);

// Solve routines
Vector back_substitution(const Matrix& A, const Vector& b);
Vector forward_substitution(const Matrix& A, const Vector& b);
Vector lu_solve(const Matrix& A, const Vector& b, const std::vector<size_t> P);

}

#endif //MATH_H
