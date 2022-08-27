#ifndef TRANSIENT_SOLVER_H
#define TRANSIENT_SOLVER_H

#include "../KEigenvalueSolver/keigenvalue_solver.h"

#include <map>


namespace InfiniteMedium
{

  /**
   * Implementation of a transient solver.
   */
  class TransientSolver : public KEigenvalueSolver
  {
  public:
    /*-------------------- Initial Conditions --------------------*/

    /**
     * A map containing group-wise initial conditions. The initial
     * conditions should be generated external to the solver.
     */
    std::map<unsigned int, double> initial_conditions;

    /*-------------------- Time Stepping --------------------*/

    /**
     * The simulation end time.
     */
    double t_end = 1.0;

    /**
     * The time step size.
     */
    double dt = 0.01;

    /*-------------------- Initialization Routines --------------------*/

    /**
     * Initialize the transient multi-group infinite medium solver.
     */
    void
    initialize() override;

    /**
     * Execute the transient multi-group infinite medium solver.
     */
    void
    execute() override;

  protected:
    /*-------------------- Solve Time Step --------------------*/

    /**
     * Execute a time step.
     */
    void
    execute_time_step();

    /**
     * Perform a transient sweep. See #sweep.
     */
    void
    transient_sweep();

    /*-------------------- Book-Keeping Quantities --------------------*/

    /**
     * The current simulation time.
     */
    double time = 0.0;

    /**
     * The multi-group angular flux from last time step. See #psi.
     */
    Vector psi_old;

    /**
     * The multi-group scalar flux from last time step. See #phi.
     */
    Vector phi_old;
  };
}

#endif //TRANSIENT_SOLVER_H
