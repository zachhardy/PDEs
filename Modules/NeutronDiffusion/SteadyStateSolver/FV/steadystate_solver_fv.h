#ifndef STEADYSTATE_SOLVER_FV_H
#define STEADYSTATE_SOLVER_FV_H

#include "NeutronDiffusion/SteadyStateSolver/Base/steadystate_solver.h"
#include "Math/Discretization/finite_volume.h"

namespace neutron_diffusion
{

class SteadyStateSolver_FV : public SteadyStateSolver
{
private:
  std::shared_ptr<math::FiniteVolume> fv;

protected:
  void initialize_discretization() override
  {
    fv = std::make_shared<math::FiniteVolume>(mesh);
    discretization = std::static_pointer_cast<Discretization>(fv);
  }

protected:
  void assemble_matrix(Groupset& groupset) override;
  void set_source(Groupset& groupset, math::Vector& b,
                  SourceFlags source_flags) override;

protected:
  void scoped_copy(const Groupset& groupset,
                   const math::Vector& x,
                   math::Vector& y) override;
};

}


#endif //STEADYSTATE_SOLVER_FV_H
