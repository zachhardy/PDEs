#ifndef STEADYSTATE_SOLVER_FV_H
#define STEADYSTATE_SOLVER_FV_H

#include "NeutronDiffusion/SteadyStateSolver/steadystate_solver.h"
#include "Math/Discretization/FiniteVolume/fv.h"

using namespace pdes::Math;


namespace NeutronDiffusion
{

class SteadyStateSolver_FV : public SteadyStateSolver
{
private:
  std::shared_ptr<FiniteVolume> fv;

protected:

  void initialize_discretization() override;

protected:
  void assemble_matrix(Groupset& groupset) override;

  void set_source(Groupset& groupset, SourceFlags source_flags) override;

  void compute_precursors() override;
};



}


#endif //STEADYSTATE_SOLVER_FV_H
