#ifndef JACOBI_H
#define JACOBI_H

#include "linear_solver.h"

namespace math
{

template<typename number>
class JacobiSolver : public SparseLinearSolver<number>
{

};
}

#endif //JACOBI_H
