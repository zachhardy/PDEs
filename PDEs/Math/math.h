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

template<typename value_type>
Vector<value_type> gaussian_elimination(Matrix<value_type>& A,
                                        Vector<value_type>& b,
                                        const bool pivot = true);
}
#endif //MATH_H
