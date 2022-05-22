#ifndef STEADYSTATE_SOLVER_FV_H
#define STEADYSTATE_SOLVER_FV_H

#include "NeutronDiffusion/SteadyStateSolver/Base/steadystate_solver.h"
#include "Math/Discretization/FiniteVolume/fv.h"

namespace NeutronDiffusion
{

class SteadyStateSolver_FV : public SteadyStateSolver
{
private:
  std::shared_ptr<Math::FiniteVolume> fv;

protected:
  void initialize_discretization() override
  {
    fv = std::make_shared<Math::FiniteVolume>(mesh);
    discretization = std::static_pointer_cast<Discretization>(fv);
  }

protected:
  void assemble_matrix(Groupset& groupset) override;

  void set_source(Groupset& groupset, Math::Vector& b,
                  SourceFlags source_flags) override;

protected:
  void scoped_transfer(const Groupset& groupset,
                       const Math::Vector& x,
                       Math::Vector& destination) override;

  void scoped_copy(const Groupset& groupset,
                   const Math::Vector& x,
                   Math::Vector& destination) override;

  double compute_change(const Groupset& groupset) override;
};

}


#endif //STEADYSTATE_SOLVER_FV_H
