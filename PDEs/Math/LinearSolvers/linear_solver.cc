#include "LinearSolvers/linear_solver.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include "macros.h"

#include <iomanip>


using namespace pdes::Math;


Vector
LinearSolver::LinearSolverBase::
solve(const Vector& b) const
{
  Vector x(b.size(), 0.0);
  solve(x, b);
  return x;
}


LinearSolver::DirectSolverBase::
DirectSolverBase(SparseMatrix& A) : A(A)
{}


LinearSolver::IterativeSolverBase::
IterativeSolverBase(const SparseMatrix& A,
                    const Options& opts,
                    const std::string solver_name) :
  A(A), tolerance(opts.tolerance),
  max_iterations(opts.max_iterations),
  verbose_history(opts.verbose_history),
  verbose_result(opts.verbose_result),
  solver_name(solver_name)
{
  Assert(A.n_rows() == A.n_cols(), "Square matrix required.");
  Assert(tolerance > 0, "Invalid tolerance specified.");
}


bool
LinearSolver::IterativeSolverBase::
check(const size_t iteration, const double value) const
{
  bool converged = value <= tolerance;

  if (verbose_history)
    std::cout << solver_name << "::"
              << "Iteration:  " << std::setw(4) << iteration << "    "
              << "Value:  " << value
              << (converged? "  CONVERGED\n" : "\n");

  if (converged && verbose_result && !verbose_history)
    std::cout << solver_name << "::CONVERGED:  "
              << "Iteration:  " << iteration << "    "
              << "Value:  " << value << std::endl;

  if (iteration == max_iterations && !converged)
    throw_convergence_error(iteration, value);

  return converged;
}


void
LinearSolver::IterativeSolverBase::
throw_convergence_error(const size_t iteration,
                        const double value) const
{
  std::stringstream err;
  err << solver_name << "!!*!! FAILURE !!*!!\n"
      << "# of Iterations:   " << iteration << std::endl
      << "Final Difference:  " << value << std::endl;
  throw std::runtime_error(err.str());
}
