#include "LinearSolvers/linear_solver.h"

#include "vector.h"
#include "sparse_matrix.h"

#include "macros.h"


using namespace pdes::Math;


Vector
LinearSolver::LinearSolverBase::
solve(const Vector& b) const
{
  Vector x(b.size(), 0.0);
  solve(x, b);
  return x;
}


LinearSolver::IterativeSolverBase::
IterativeSolverBase(const SparseMatrix& A,
                    const Options& opts,
                    const std::string solver_name) :
  A(A), tolerance(opts.tolerance),
  max_iterations(opts.max_iterations),
  verbose(opts.verbose), solver_name(solver_name)
{
  Assert(A.n_rows() == A.n_cols(), "Square matrix required.");
  Assert(tolerance > 0, "Invalid tolerance specified.");
}


void
LinearSolver::IterativeSolverBase::
throw_convergence_error(const size_t iteration,
                        const double difference) const
{
  std::stringstream err;
  err << solver_name << " Solver did not converge!\n"
      << "# of Iterations:   " << iteration << std::endl
      << "Final Difference:  " << difference << std::endl;
  throw std::runtime_error(err.str());
}
