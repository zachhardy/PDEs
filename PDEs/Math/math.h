#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <vector>

namespace math
{

// Forward declarations
template<typename value_type>
class Vector;

template<typename value_type>
class Matrix;

// Factorizations
void row_echelon_form(Matrix<double>& A,
                      Vector<double>& b,
                      const bool pivot = true);

std::vector<size_t> lu_factorization(Matrix<double>& A,
                                     const bool pivot = true);

void cholesky_factorization(Matrix<double>& A);

// Solve routines
Vector<double> back_substitution(const Matrix<double>& A,
                                 const Vector<double>& b);
Vector<double> forward_substitution(const Matrix<double>& A,
                                    const Vector<double>& b);

Vector<double> gaussian_elimination(Matrix<double>& A,
                                    Vector<double>& b,
                                    const bool pivot = true);

Vector<double> lu_solve(const Matrix<double>& A,
                        const Vector<double>& b,
                        const std::vector<size_t> P);

Vector<double> cholesky_solve(const Matrix<double>& A,
                              const Vector<double>& b);

}
#endif //MATH_H
