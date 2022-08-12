#ifndef TRANSIENT_SOLVER_H
#define TRANSIENT_SOLVER_H

#include "../KEigenvalueSolver/keigenvalue_solver.h"

#include <map>
#include <functional>


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

  //############################################################

  /**
   * Implementation of a transient neutron diffusion solver.
   */
  class TransientSolver : public KEigenvalueSolver
  {
  protected:
    typedef TimeSteppingMethod TSMethod;
    typedef NormalizationMethod NormMethod;

    /*-------------------- Constants --------------------*/
  protected:
    /**
     * Energy release per fission (J/fission).
     */
    const double energy_per_fission = 3.204e-11;

    /**
     * A conversion factor to convert fission energy release to a change in
     * material temperature (K (\f$ cm^3 \f$).
     */
    const double conversion_factor = 3.83e-11;

    /*-------------------- Options --------------------*/
  public:
    /**
     * A flag for normalizing fission cross-sections to a precomputed
     * k-eigenvalue.
     */
    bool normalize_fission_xs = false;

    /**
     * A flag for lagging precursors to the previous time step.
     */
    bool lag_precursors = false;

    /**
     * A flag for whether or not the problem has dynamic cross-sections
     * or not. This is used to decide whether to call update functions, where
     * necessary.
     */
    bool has_dynamic_xs = false;

    /**
     * The simulation start time.
     */
    double t_start = 0.0;

    /**
     * The simulation end time.
     */
    double t_end =1.0;

    /**
     * The nominal time step size.
     */
    double dt = 0.1;

    /**
     * The minimum time step size when using adaptive time stepping.
     */
    double dt_min = 1.0e-6;

    /**
     * The time discretization methd.
     */
    TSMethod time_stepping_method = TSMethod::CRANK_NICHOLSON;

    /**
     * A flag for adaptive time stepping.
     */
    bool adaptive_time_stepping = false;

    /**
     * The relative power change that triggers refinement when exceeded.
     */
    double refine_threshold = 0.01;

    /**
     * The relative power change that triggers coarsening when exceeded.
     */
    double coarsen_threshold = 0.05;

    /**
     * A convenient typedef for group-wise initial condition functions. Each
     * function should take a point object as input and returns a the
     * evaluated initial condition.
     */
    typedef std::function<double(const Grid::Point p)> InitialCondition;

    /**
     * A map holding the initial condition functions used to initialize
     * the transient solver. If the map is empty, a k-eigenvalue solver is
     * used for the initial condition. Otherwise, the provided initial
     * conditions are evaluated. Groups which do not have a specified initial
     * condition are assumed to be uniformly zero.
     */
    std::map<size_t, InitialCondition> initial_conditions;

    /**
     * The normalization method for the initial condition. This is used to map
     * the magnitude of the flux profile to a specified reactor power or
     * average power density.
     */
    NormMethod normalization_method = NormMethod::TOTAL_POWER;

    /**
     * A flag for whether to write the time step results to output files.
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
     * The path to the output directory to write outputs to.
     */
    std::string output_directory;


    /*-------------------- Macro System Data --------------------*/

    /**
     * The current simulation time.
     */
    double time = 0.0;

    /**
     * The current reactor power.
     */
    double power = 1.0;

    /**
     * The reactor power from last time step.
     */
    double power_old = 1.0;

    /**
     * The average reactor power density.
     */
    double average_power_density;

    /**
     * The peak power density in the reactor.
     */
    double peak_power_density;

    /**
     * The initial fuel temperature in K.
     */
    double initial_temperature = 300.0;

    /**
     * The average temperature in the fuel.
     */
    double average_fuel_temperature;

    /**
     * The peak temperature in the fuel.
     */
    double peak_fuel_temperature;

    /*-------------------- System Vectors --------------------*/

    /**
     * The fission rate defined at cell centers.
     */
    Vector fission_rate;

    /**
     * The temperature defined at cell centers.
     */
    Vector temperature;

  protected:
    /**
     * The multi-group scalar flux from last time step.
     */
    Vector phi_old;

    /**
     * The precursor concentrations from last time step.
     */
    Vector precursor_old;

    /**
     * The temperature from last time step.
     */
    Vector temperature_old;

    /*-------------------- Interface Routines --------------------*/
  public:
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
    write(const size_t output_index) const;

    /*-------------------- Initialization --------------------*/
  protected:
    /**
     * Compute the initial conditions for the transient.
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
     * the \p refine_threshold attribute. When this routine is called, the
     * previous time step results are discarded and the time step is rerun until
     * the new time step is accepted.
     */
    void
    refine_time_step();

    /**
     * Coarsen the time step when teh relative change in power is less than the
     * \p coarsen_threshold attribute. When this routine is called, the previous
     * time step is accepted and the increase time step size goes into effect
     * in the following time step.
     */
    void
    coarsen_time_step();

    /**
     * Set the current time step solutions to the old time step solutions to
     * prepare for the next time step.
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
     * Rebuild the transient matrix.
     */
    void
    rebuild_matrix();

    /*-------------------- Auxiliary Quantities --------------------*/

    /**
     * Compute the fission rate using the most up-to-date multi-group
     * scalar flux solution.
     */
    void
    update_fission_rate();

    /**
     * Compute the precursor concentrations using the most up-to-date
     * multi-group scalar flux solution.
     */
    void
    update_precursors();

    /**
     * Compute the temperature profile using the most up-to-date multi-group
     * scalar flux solution.
     */
    void
    update_temperature();

    /**
     * Compute the bulk properties using the most up-to-date multi-group
     * scalar flux solution.
     */
    void
    compute_bulk_properties();

    /**
     * Return the effective time step size based on the time stepping method.
     *
     * \return The effective time step size for the linear system.
     */
    double
    effective_time_step();
  };
}

#endif //TRANSIENT_SOLVER_H
