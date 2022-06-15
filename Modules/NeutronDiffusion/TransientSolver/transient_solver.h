#ifndef TRANSIENT_SOLVER_H
#define TRANSIENT_SOLVER_H

#include "../KEigenvalueSolver/keigenvalue_solver.h"

#include <map>
#include <functional>


namespace Math
{
  enum class TimeSteppingMethod
  {
    BACKWARD_EULER = 0,
    CRANK_NICHOLSON = 1,
    TBDF2 = 2
  };
}


//######################################################################


using namespace Math;


namespace NeutronDiffusion
{

  /** Normalization method for the initial condition. */
  enum class NormalizationMethod
  {
    NONE = 0,
    TOTAL_POWER = 1,
    AVERAGE_POWER = 2
  };


  /** Implementation of a transient neutron diffusion solver. */
  class TransientSolver : public KEigenvalueSolver
  {
  public:

    /*-------------------- Constants --------------------*/

    /** Energy release per fission (J/fission). */
    const double energy_per_fission = 3.2e-11;

    /**
     * A conversion factor to convert fission energy release to a change in
     * material temperature (K (\f$ cm^3 \f$).
     */
    const double conversion_factor = 3.83e-11;

    /*-------------------- Options --------------------*/

    /**
     * A flag for normalizing fission cross-sections to a precomputed
     * k-eigenvalue.
     */
    bool normalize_fission_xs = false;

    /**
     * A flag for lagging precursors or using the same temporal discretization
     * method as the scalar flux.
     */
    bool lag_precursors = false;

    /*-------------------- Time Stepping --------------------*/

    double time = 0.0;
    double dt = 0.1;

    double t_start = 0.0;
    double t_end = 1.0;

    typedef TimeSteppingMethod SteppingMethod;
    SteppingMethod time_stepping_method = SteppingMethod::BACKWARD_EULER;

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
    NormalizationMethod normalization_method = NormalizationMethod::TOTAL_POWER;

    /*-------------------- Outputting Options --------------------*/

    /** A flag for whether to write time step results. */
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

    /** The path to the directory to write output files. */
    std::string output_directory;

    /*-------------------- Solutions --------------------*/

    Vector phi_old;
    Vector precursor_old;

    /*-------------------- Power Quantities --------------------*/

    Vector fission_rate;

    double power = 1.0;
    double power_old = 1.0;

    double average_power_density;
    double peak_power_density;

    /*-------------------- Temperature Quantities --------------------*/

    Vector temperature;
    Vector temperature_old;

    double initial_temperature = 300.0;

    double average_fuel_temperature;
    double peak_fuel_temperature;

    /*-------------------- Public Facing Routines --------------------*/

    void initialize() override;
    void execute() override;

  private:

    void compute_initial_values();
    void evaluate_initial_conditions();

    /*-------------------- Time Step Routines --------------------*/

    void solve_time_step();
    void solve_groupset_time_step(Groupset& groupset,
                                  SourceFlags source_flags);
    void solve_full_system_time_step(SourceFlags source_flags);

    void step_solutions();

    /*-------------------- Assembly Routines --------------------*/

    void assemble_transient_matrix(Groupset& groupset,
                                   AssemblerFlags assembler_flags);

    void set_transient_source(Groupset& groupset,
                              SourceFlags source_flags);

    void assemble_matrices();

    /*-------------------- Auxiliary Quantities --------------------*/

    void compute_fission_rate();
    void compute_power();

    void update_precursors();

    double effective_time_step();
  };
}

#endif //TRANSIENT_SOLVER_H
