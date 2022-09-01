#include "LinearSolvers/linear_solver.h"

#include "vector.h"
#include "matrix.h"
#include "Math/sparse_matrix.h"

#include <iomanip>
#include <cassert>


using namespace PDEs;
using namespace Math;
using namespace LinearSolvers;


template<class MatrixType>
Vector
LinearSolverBase<MatrixType>::solve(const Vector& b) const
{
  Vector x(b.size(), 0.0);
  solve(x, b);
  return x;
}


template<class MatrixType>
void
LinearSolverBase<MatrixType>::set_matrix(const MatrixType& matrix)
{
  assert(matrix.n_rows() == matrix.n_cols());
}


template class LinearSolvers::LinearSolverBase<Matrix>;
template class LinearSolvers::LinearSolverBase<SparseMatrix>;


//################################################## DirectSolverBase


template<class MatrixType>
DirectSolverBase<MatrixType>::
DirectSolverBase() :
    LinearSolverBase<MatrixType>()
{}


template<class MatrixType>
void
DirectSolverBase<MatrixType>::set_matrix(const MatrixType& matrix)
{
  LinearSolverBase<MatrixType>::set_matrix(matrix);
  A = matrix;
  factorize();
}

template class LinearSolvers::DirectSolverBase<Matrix>;
template class LinearSolvers::DirectSolverBase<SparseMatrix>;


Options::Options(const double tolerance,
                 const unsigned int max_iterations,
                 const unsigned int verbosity) :
    tolerance(tolerance),
    max_iterations(max_iterations),
    verbosity(verbosity)
{}


IterativeSolverBase::
IterativeSolverBase(const Options& opts, const std::string name) :
    tolerance(opts.tolerance),
    max_iterations(opts.max_iterations),
    verbosity(opts.verbosity),
    solver_name(name)
{}


void
IterativeSolverBase::
set_matrix(const SparseMatrix& matrix)
{
  LinearSolverBase<SparseMatrix>::set_matrix(matrix);
  A = &matrix;
}


bool
IterativeSolverBase::
check(const unsigned int iteration, const double value) const
{
  bool converged = value <= tolerance;

  if (verbosity > 1)
    std::cout << solver_name << "::"
              << "Iteration   " << std::setw(4) << iteration << "    "
              << "Value   " << value
              << (converged ? "  CONVERGED\n" : "\n");

  if (converged && verbosity == 1)
    std::cout << solver_name << "::  CONVERGED   "
              << "Iteration   " << iteration << "    "
              << "Value   " << value << std::endl;

  if (iteration == max_iterations && !converged)
    throw_convergence_error(iteration, value);

  return converged;
}


void
IterativeSolverBase::
throw_convergence_error(const unsigned int iteration,
                        const double value) const
{
  std::stringstream err;
  err << "!!*!! " << solver_name << " FAILURE !!*!!\n"
      << "# of Iterations:   " << iteration << std::endl
      << "Final Difference:  " << value << std::endl;
  throw std::runtime_error(err.str());
}
