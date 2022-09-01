#ifndef TRANSIENT_SOLVER_H
#define TRANSIENT_SOLVER_H

#include "../KEigenvalueSolver/keigenvalue_solver.h"

#include <map>
#include <functional>


namespace PDEs
{
  namespace Math
  {
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
  protected:
    using TSMethod = TimeSteppingMethod;
    using NormMethod = NormalizationMethod;

  public:
    /**
     * A convenient typedef for group-wise initial condition functions. Each
     * function should take a point object as input and returns a the
     * evaluated initial condition.
     */
    using IC = std::function<double(const Grid::CartesianVector p)>;

  public:

    /*-------------------- General Options --------------------*/

    bool write_outputs = false;
    std::string output_directory;

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

    bool lag_precursors = false;

    /**
     * The normalization method for the initial condition. This is used to map
     * the magnitude of the flux profile to a specified reactor power or
     * average power density.
     */
    NormMethod normalization_method = NormMethod::TOTAL_POWER;

    /** A flag for normalizing \f$ \sigma_f \f$ to the \f$ k \f$-eigenvalue. */
    bool normalize_fission_xs = false;

    /*-------------------- Initial Conditions --------------------*/

    double initial_power = 1.0;
    double initial_temperature = 300.0;

    /**
     * A map holding the initial condition functions used to initialize
     * the transient solver. If the map is empty, a k-eigenvalue solver is
     * used for the initial condition. Otherwise, the provided initial
     * conditions are evaluated. Groups which do not have a specified initial
     * condition are assumed to be uniformly zero.
     */
    std::map<unsigned int, IC> initial_conditions;

    /*-------------------- Time Stepping --------------------*/

    double t_start = 0.0;
    double t_end = 1.0;
    double dt = 0.1;

    TSMethod time_stepping_method = TSMethod::CRANK_NICHOLSON;

    /**
     * A flag for whether to use adaptive time stepping or not. Currently,
     * adaptive time stepping is based only on the relative change in power
     * over a time step. When the change in power is greater than the \p
     * refine_threshold, the time step is halved. When it is less than the \p
     * coarsen_threshold, it is doubled. If the time-step becomes smaller than
     * \p dt_min, it is defaulted to \p dt_min and not refined further.
     */
    bool adaptive_time_stepping = false;
    double refine_threshold = 0.05;
    double coarsen_threshold = 0.01;
    double dt_min = 1.0e-6;

  protected:
    /*-------------------- Problem Information --------------------*/

    /**
     * A flag for whether or not the problem has dynamic cross-sections
     * or not. This is used to decide whether to call update functions, where
     * necessary.
     */
    bool has_dynamic_xs = false;

    /**
     * Energy release per fission (J/fission).
     */
    const double energy_per_fission = 3.204e-11;

    /**
     * A conversion factor to convert fission energy release to a change in
     * material temperature (K (\f$ cm^3 \f$).
     */
    const double conversion_factor = 3.83e-11;

    double time = 0.0;

    double power = 1.0;
    double power_old = 1.0;


    double average_power_density;
    double peak_power_density;

    double average_fuel_temperature;
    double peak_fuel_temperature;


    Vector phi_old; ///< The multi-group scalar flux last time step.
    Vector precursors_old; ///< The precursor concentrations last time step.

    Vector fission_rate; ///< The cell-wise fission rate.

    Vector temperature; ///< The cell-wise temperature.
    Vector temperature_old; ///< The temperature last time step.

  public:
    /*-------------------- Public Routines --------------------*/

    void initialize() override;
    void execute() override;

    /**
     * Write the current simulation state to an output file at the specified
     * \p output_frequency.
     *
     * This routine writes the simulation state to an output file within
     * the specified \p output_directory. The \p output_index corresponds to
     * the number of outputs that have already been printed. The output files
     * are saved as <tt><output_index>.data</tt> where the integer is padded
     * with to ensure that all output file names have the same number of
     * characters.
     */
    void write(const unsigned int output_index) const;

  protected:
    /*-------------------- Initialization Routines --------------------*/

    /**
     * Evaluate the specified initial conditions.
     *
     * If \p initial_conditions is empty, the initial condition will default
     * to the result of a \f$ k \f$-eigenvalue simulation, otherwise, the
     * group-wise initial conditions are evaluated. This routine also normalizes
     * the multi-group scalar flux to the specified power level and the fission
     * cross-sections to the \f$ k \f$-eigenvalue, if specified.
     */
    void compute_initial_values();

    /*-------------------- Time Step Routines --------------------*/

    /**
     * Execute a time step. When the \p reconstruct_matrices flag is \p true,
     * this routine will rebuild the multi-group operator using the current
     * simulation state.
     */
    void execute_time_step(bool reconstruct_matrices = false);

    /**
     * Lag the specified sources within \p source_flags and iteratively solve
     * the multi-group system over a time step.
     */
    void iterative_time_step_solve(SourceFlags source_flags);

    /**
     * Refine the time step when the relative change in power is greater than
     * \p refine_threshold.
     *
     * When this routine is called, the previous time step
     * results are discarded and the time step is rerun until the new time step
     * is accepted. If the time step sizes becomes smaller than \p dt_min, the
     * refinement is stopped.
     */
    void refine_time_step();

    /**
     * Coarsen the time step when the relative change in power is less than the
     * \p coarsen_threshold. When this routine is called, the previous  time
     * step is accepted and the increase time step size goes into effect in the
     * following time step.
     */
    void coarsen_time_step();

    /** Set the last time step quantities to the current values. */
    void step_solutions();

    /*-------------------- Assembly Routines --------------------*/

    /**
     * Assemble the multi-group matrix according to the specified
     * \p assembler_flags.
     *
     * By default, the within-group terms (time derivative, total interaction,
     * buckling, diffusion, and boundary) are included in the matrix. When
     * specified, the cross-group scattering and fission terms may be included.
     */
    void assemble_transient_matrix(
        AssemblerFlags assembler_flags = NO_ASSEMBLER_FLAGS);

    /**
     * Accumulate sources into the right-hand side according to the specified
     * \p source_flags.
     *
     * This routine is additive, implying that if the right-hand side needs to
     * be cleared, this must be done external to this routine. Available
     * sources are the inhomogeneous, scattering, fission, and boundary source.
     */
    void set_transient_source(SourceFlags source_flags);

    /**
     * Rebuild the multi-group matrix.
     *
     * This is shorthand for the conditional construction based on the
     * \p algorithm option.
     */
    void rebuild_matrix();

    /*-------------------- Auxiliary Quantities --------------------*/

    void update_fission_rate();
    void update_precursors();
    void update_temperature();

    /** Update the averaged and peak powers and temperatures. */
    void compute_bulk_properties();

    /**
     * Return the effective time step size based on the time stepping method.
     *
     * For example, when using Crank-Nicholson, the effective time step is
     * half the true time step.
     */
    double effective_time_step();

    /*-------------------- Write Routines --------------------*/

    void write_snapshot(const unsigned int index) const;
    void write_temperature(const std::string directory = ".",
                           const std::string file_prefix = "temperature");

  };
}

#endif //TRANSIENT_SOLVER_H
