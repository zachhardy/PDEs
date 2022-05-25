#include "steadystate_solver_fv.h"
#include "Discretization/FiniteVolume/fv.h"

void
NeutronDiffusion::SteadyStateSolver_FV::
initialize_discretization()
{
  fv = std::make_shared<FiniteVolume>(mesh);
  discretization = std::static_pointer_cast<Discretization>(fv);
}
