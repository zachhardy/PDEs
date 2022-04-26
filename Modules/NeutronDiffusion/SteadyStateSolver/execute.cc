#include "steadystate_solver.h"

#include <fstream>

/// Run the steady state multigroup diffusion simulation.
void neutron_diffusion::SteadyStateSolver::execute()
{
  std::cout << "\nExecuting solver...\n";

  switch (solution_method)
  {
    case SolutionMethod::DIRECT:
    {
      assemble_matrix();
      set_source();

      linear_solver->setup();
      phi = linear_solver->solve(system_rhs);
      break;
    }

    case SolutionMethod::ITERATIVE:
    {
      assemble_matrix();
      linear_solver->setup();

      size_t nit = 0;
      double diff = 1.0;
      math::Vector phi_ell = phi;
      while (diff > 1.0e-8 and nit < 100)
      {
        ++nit;

        set_source();
        phi = linear_solver->solve(system_rhs);

        diff = math::l2_norm(phi - phi_ell);
        phi_ell = phi;

        std::cout << "Iteration Number:  " << std::setw(4) << nit << " "
                  << "Differrence:       " << diff << "\n";
      }
      break;
    }
  }

  write_solution("solution", phi);

  std::cout << "\nDone executing solver.\n";
}

//######################################################################

/** \brief Assemble the matrix with the routine associated with the
 *         discretization type. */
void neutron_diffusion::SteadyStateSolver::assemble_matrix()
{
  switch (discretization->type)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    { assemble_fv_matrix(); break; }
    default: return;
  }
}

/** \brief Assemble the right-hand side vector with the routine associated with
 *         the discretization type. */
void neutron_diffusion::SteadyStateSolver::set_source()
{
  switch (discretization->type)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    { set_fv_source(); break; }
    default: return;
  }
}
