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

    bool write_outputs = false;
    std::string output_directory;

    bool lag_density = false;

    /*-------------------- Initial Conditions --------------------*/

    /**
     * A map containing group-wise initial conditions. The initial
     * conditions should be generated external to the solver.
     */
    std::map<unsigned int, double> initial_conditions;

    /*-------------------- Time Stepping --------------------*/

    double t_end = 1.0;
    double dt = 0.01;

  protected:
    double time = 0.0;

    Vector psi_old;
    Vector phi_old;

  public:
    /*-------------------- Public Routines --------------------*/

    void initialize() override;
    void execute() override;

  protected:
    /*-------------------- Solve Time Step --------------------*/

    void execute_time_step();
    void transient_sweep();
    double balance_check();

    /*-------------------- Write Operations --------------------*/

    void write_snapshot(const unsigned int index) const;

  };
}

#endif //TRANSIENT_SOLVER_H
