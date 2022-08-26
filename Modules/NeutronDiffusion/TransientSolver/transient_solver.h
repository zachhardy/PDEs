#ifndef TRANSIENT_SOLVER_H
#define TRANSIENT_SOLVER_H

#include "../KEigenvalueSolver/keigenvalue_solver.h"

#include <map>
#include <functional>


namespace PDEs
{
  namespace Math
  {
    /**
     * Time-stepping method options.
     */
    enum class TimeSteppingMethod
    {
      BACKWARD_EULER = 0,
      CRANK_NICHOLSON = 1,
      TBDF2 = 2
    };
  }
}


using namespace PDEs;
using namespace Math;


namespace NeutronDiffusion
{

  /**
   * Normalization method for the initial condition.
   */
  enum class NormalizationMethod
  {
    NONE = 0, ///< Use the raw initial condition data.
    TOTAL_POWER = 1, ///< Normalize to the total reactor power.
    AVERAGE_POWER = 2 ///< Normalize to the avererage power density.
  };


  /**
   * Implementation of a transient neutron diffusion solver.
   */
  class TransientSolver : public KEigenvalueSolver
  {
  public:
    /**
     * Shorthand for the time stepping method.
     */
    using TSMethod = TimeSteppingMethod;

    /**
     * Shorthand for the initial condition normalization method.
     */
    using NormMethod = NormalizationMethod;

    /**
     * A convenient typedef for group-wise initial condition functions. Each
     * function should take a point object as input and returns a the
     * evaluated initial condition.
     */
    using IC = std::function<double(const Grid::CartesianVector p)>;

    /*-------------------- Outputting --------------------*/

    /**
     * A flag for writing outputs. When true, outputs are written for each
     * time step.
     */
    bool write_outputs = false;

    /**
     * The frequency that outputs should be written. If -1.0, the solution
     * is written at each time step, otherwise, solutions are written at
     * the specified interval. If adaptive time stepping is on, time steps
     * are modified to ensure coincidence with the output times. If adaptive
     * time stepping is off, and the output frequency is not a multiple of
     * the time step size, an error is thrown. If the specified initial time
     * step size is smaller than the output frequency, the time step is
     * changed to the output frequency.
     */
    double output_frequency = -1.0;

    /**
     * The output directory where the time step snapshots should be written.
     */
    std::string output_directory;

    /*-------------------- Precursor Options --------------------*/

    /**
     * A flag for lagging precursors to the previous time step.
     */
    bool lag_precursors = false;

    /*-------------------- Initial Conditions --------------------*/

    /**
     * A map holding the initial condition functions used to initialize
     * the transient solver. If the map is empty, a k-eigenvalue solver is
     * used for the initial condition. Otherwise, the provided initial
     * conditions are evaluated. Groups which do not have a specified initial
     * condition are assumed to be uniformly zero.
     */
    std::map<size_t, IC> initial_conditions;

    /**
     * The normalization method for the initial condition. This is used to map
     * the magnitude of the flux profile to a specified reactor power or
     * average power density.
     */
    NormMethod normalization_method = NormMethod::TOTAL_POWER;

    /**
     * A flag for normalizing fission cross-sections to the \f$k\f$-eigenvalue.
     */
    bool normalize_fission_xs = false;

    /**
     * The initial power level in the reactor.
     */
    double initial_power = 1.0;

    /**
     * The initial temperature in the reactor.
     */
    double initial_temperature = 300.0;

    /*-------------------- Time Stepping --------------------*/

    /**
     * The simulation start time.
     */
    double t_start = 0.0;

    /**
     * The simulation end time.
     */
    double t_end =1.0;

    /**
     * The initial time step size to use. If adaptivity is off, this remains the
     * time step size through the simulation.
     */
    double dt = 0.1;

    /**
     * The minimum time step to allow during the simulation. This is only
     * relevant when adaptive time stepping is turned on.
     */
    double dt_min = 1.0e-6;

    /**
     * The time stepping method to use.
     */
    TSMethod time_stepping_method = TSMethod::CRANK_NICHOLSON;

    /**
     * A flag for whether to use adaptive time stepping or not. Currently,
     * adaptive time stepping is based only on the relative change in power
     * over a time step. When certain bounds are exceeded, the time step size
     * is either halved or doubled.
     */
    bool adaptive_time_stepping = false;

    /**
     * The relative power change that triggers refinement.
     */
    double refine_threshold = 0.01;

    /**
     * The relative power change that triggers coarsening.
     */
    double coarsen_threshold = 0.05;

    /**
     * Initialize the transient multi-group diffusion solver.
     */
    void
    initialize() override;

    /**
     * Execute the transient multi-group diffusion solver.
     */
    void
    execute() override;

    /**
     * Write the current simulation state to an output file.
     *
     * \param output_index The output number. This defines the file name that
     *      the system state is saved to.
     */
    void
    write(const unsigned int output_index) const;

  protected:
    /**
     * Evaluate the specified initial conditions.
     */
    void
    compute_initial_values();

    /*-------------------- Time Step Routines --------------------*/

    /**
     * Execute a time step by computing the end of time step solutions.
     *
     * \param reconstruct_matrices A flag for whether the matrices need to be
     *      reconstructed this time step or not.
     */
    void
    execute_time_step(bool reconstruct_matrices = false);

    /**
     * Solve the time step system by iterating on the specified SourceFlags.
     *
     * \param source_flags Bitwise flags defining the source terms to iterate
     *      on and converge.
     */
    void
    iterative_time_step_solve(SourceFlags source_flags);

    /**
     * Refine the time step when the relative change in power is greater than
     * \p refine_threshold. When this routine is called, the previous time step
     * results are discarded and the time step is rerun until the new time step
     * is accepted.
     */
    void
    refine_time_step();

    /**
     * Coarsen the time step when the relative change in power is less than the
     * \p coarsen_threshold. When this routine is called, the previous  time
     * step is accepted and the increase time step size goes into effect in the
     * following time step.
     */
    void
    coarsen_time_step();

    /**
     * Set the last time step quantities to the current values.
     */
    void
    step_solutions();

    /*-------------------- Assembly Routines --------------------*/

    /**
     * Assemble the multi-group matrix according to the specified
     * AssemblerFlags. By default, the within-group terms (time derivative,
     * total interaction, buckling, diffusion, and boundary) are included in
     * the matrix. When specified, the cross-group scattering and fission terms
     * may be included.
     *
     * \param assembler_flags Bitwise flags used to specify which cross-group
     *      terms to include in the matrix.
     */
    void
    assemble_transient_matrix(AssemblerFlags assembler_flags);


    /**
     * Set the right-hand side source vector. This is an additive routine which
     * will only add the specified sources to the source vector. Source options
     * include the inhomogeneous source, scattering source, and fission source.
     * By default, the previous time step contributions are added.
     *
     * \param source_flags Bitwise flags used to specify which sources are
     *      added to the source vector.
     */
    void
    set_transient_source(SourceFlags source_flags);

    /**
     * Rebuild the system matrix. This is a shorthand notation to avoid
     * repetitive condition checking throughout the code.
     */
    void
    rebuild_matrix();

    /*-------------------- Auxiliary Quantities --------------------*/

    /**
     * Compute the fission rate with the most recent multi-group scalar flux.
     */
    void
    update_fission_rate();

    /**
     * Take a time step for the delayed neutron precursors concentrations using
     * the most recent multi-group scalar flux.
     */
    void
    update_precursors();

    /**
     * Take a time step for the adiabatic model using the most recent
     * multi-group scalar flux.
     */
    void
    update_temperature();

    /** Update the averaged and peak powers and temperatures. */
    void
    compute_bulk_properties();

    /**
     * Return the effective time step size based on the time stepping method.
     * For example, when using Crank-Nicholson, the effective time step is
     * half the true time step.
     */
    double
    effective_time_step();

    /**
     * A flag for whether or not the problem has dynamic cross-sections
     * or not. This is used to decide whether to call update functions, where
     * necessary.
     */
    bool has_dynamic_xs = false;

    /*-------------------- Constants --------------------*/

    /**
     * Energy release per fission (J/fission).
     */
    const double energy_per_fission = 3.204e-11;

    /**
     * A conversion factor to convert fission energy release to a change in
     * material temperature (K (\f$ cm^3 \f$).
     */
    const double conversion_factor = 3.83e-11;

    /*-------------------- Macro Quantities --------------------*/

    /**
     * The current simulation time.
     */
    double time = 0.0;

    /**
     * The current total reactor power.
     */
    double power = 1.0;

    /**
     * The total reactor power last time step.
     */
    double power_old = 1.0;

    /**
     * The current average power density in the reactor.
     */
    double average_power_density;

    /**
     * The current peak power density in the reactor.
     */
    double peak_power_density;

    /**
     * The current average fuel temperature in the reactor.
     */
    double average_fuel_temperature;

    /**
     * The current peak fuel temperature in the reactor.
     */
    double peak_fuel_temperature;

    /*-------------------- Extra System Vectors --------------------*/

    /**
     * The multi-group scalar flux last time step. See #phi
     */
    Vector phi_old;

    /**
     * The delayed neutron precursors last time step. See #precursor.
     */
    Vector precursors_old;

    /**
     * The fission rate defined at cell centers.
     */
    Vector fission_rate;

    /**
     * The temperature defined at cell centers.
     */
    Vector temperature;

    /**
     * The temperature last time step.
     */
    Vector temperature_old;
  };
}

#endif //TRANSIENT_SOLVER_H
