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
void cholesky_factorization(Matrix& A);

// Solve routines
Vector back_substitution(const Matrix& A, const Vector& b);
Vector forward_substitution(const Matrix& A, const Vector& b);

Vector gaussian_elimination(Matrix& A, Vector& b, const bool pivot = true);
Vector lu_solve(const Matrix& A, const Vector& b, const std::vector<size_t> P);
Vector cholesky_solve(const Matrix& A, const Vector& b);

}
#endif //MATH_H
