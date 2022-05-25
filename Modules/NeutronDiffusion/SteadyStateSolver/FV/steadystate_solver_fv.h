#ifndef STEADYSTATE_SOLVER_FV_H
#define STEADYSTATE_SOLVER_FV_H

#include "NeutronDiffusion/SteadyStateSolver/Base/steadystate_solver.h"
#include "Math/Discretization/FiniteVolume/fv.h"

using namespace pdes::Math;


namespace NeutronDiffusion
{

class SteadyStateSolver_FV : public SteadyStateSolver
{
private:
  std::shared_ptr<FiniteVolume> fv;

protected:

  void
  initialize_discretization() override;

protected:
  void
  assemble_matrix(Groupset& groupset) override;

  void
  set_source(Groupset& groupset, Vector& b,
             SourceFlags source_flags) override;
};

}


#endif //STEADYSTATE_SOLVER_FV_H
