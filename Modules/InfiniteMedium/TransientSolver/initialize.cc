#include "transient_solver.h"


using namespace InfiniteMedium;


void TransientSolver::initialize()
{
  KEigenvalueSolver::initialize();

  psi_old = psi;
  phi_old = phi;

  // Set the initial condition
  for (const auto& ic : initial_conditions)
    for (unsigned int n = 0; n < n_angles; ++n)
      psi_old[n*n_groups + ic.first] = ic.second;
}
