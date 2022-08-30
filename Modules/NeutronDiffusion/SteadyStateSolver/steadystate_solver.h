#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "../boundaries.h"

#include "mesh.h"
#include "Discretization/discretization.h"

#include "vector.h"
#include "Math/sparse_matrix.h"
#include "LinearSolvers/linear_solver.h"

#include "material.h"
#include "CrossSections/cross_sections.h"
#include "lightweight_xs.h"

#include <string>


using namespace PDEs;
using namespace Grid;
using namespace Math;
using namespace Physics;


namespace NeutronDiffusion
{
  /**
   * Algorithms available to solve the multi-group system.
   */
  enum class Algorithm
  {
    DIRECT = 0,   ///< Solve the full multi-group system.
    ITERATIVE = 1  ///< Iterate on the cross-group terms.
  };


  /**
   * Bitwise source flags for right-hand side vector construction.
   */
  enum SourceFlags : int
  {
    NO_SOURCE_FLAGS = 0,
    APPLY_MATERIAL_SOURCE = (1 << 0),
    APPLY_SCATTER_SOURCE = (1 << 1),
    APPLY_FISSION_SOURCE = (1 << 2),
    APPLY_BOUNDARY_SOURCE = (1 << 3)
  };

  /**
   * Bitwise source flag operator.
   */
  inline SourceFlags
  operator|(const SourceFlags f1, const SourceFlags f2)
  {
    return static_cast<SourceFlags>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
  }


  /**
   * Bitwise assembler flags for matrix construction.
   */
  enum AssemblerFlags : int
  {
    NO_ASSEMBLER_FLAGS = 0,
    ASSEMBLE_SCATTER = (1 << 0),
    ASSEMBLE_FISSION = (1 << 1)
  };


  /**
   * Bitwise assembler flag operator.
   */
  inline AssemblerFlags
  operator|(const AssemblerFlags f1, const AssemblerFlags f2)
  {
    return static_cast<AssemblerFlags>(static_cast<int>(f1) |
                                       static_cast<int>(f2));
  }


  /**
   * A steady-state multi-group diffusion solver.
   */
  class SteadyStateSolver
  {
  protected:
    using SDMethod = SpatialDiscretizationMethod;
    using IsotropicMGSource = IsotropicMultiGroupSource;
    using BndryPtr = std::shared_ptr<Boundary>;
    using LinearSolver = LinearSolvers::LinearSolverBase<SparseMatrix>;

  public:
    /*-------------------- General Options --------------------*/

    Algorithm algorithm = Algorithm::DIRECT;
    SDMethod discretization_method = SDMethod::FINITE_VOLUME;

    /**
     * A flag for whether or not delayed neutron precursors should be used
     * or not. If delayed neutron data is present and this flag is false, the
     * total fission data derived from the prompt and delayed quantities are
     * used.
     */
    bool use_precursors = false;


    double inner_tolerance = 1.0e-6;
    unsigned int max_inner_iterations = 100;

    unsigned int verbosity = 0;

    /*-------------------- Spatial Domain --------------------*/

    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<Discretization> discretization;

    /*-------------------- Physics Information --------------------*/

     std::vector<std::shared_ptr<Material>> materials;

     /**
      * A list of the group IDs to be used in the simulation. While this may
      * seem useless, it allows for a subset of available groups within
      * a cross-section library to be considered.
      */
     std::vector<unsigned int> groups;

     /*-------------------- Boundary Information --------------------*/

    /**
     * A list specifying the different boundary conditions in the problem.
     * Each boundary condition is given by a pair. The first element contains
     * the boundary type. The second contains an index which points to the
     * position within the auxiliary \p boundary_values vector that the
     * multi-group boundary values are located.
     */
    std::vector<std::pair<BoundaryType, unsigned int>> boundary_info;

    /**
     * The multi-group boundary values. The outer index is used by the
     * \p boundary_info attribute to access a boundary value, the middle points
     * to particular energy groups, and the inner to the boundary values. For
     * non-Robin boundaries, this always has one entry. For Robin boundaries,
     * three entries in the order of <tt>(a, b, f)</tt> are used.
     */
    std::vector<std::vector<std::vector<double>>> boundary_values;

    /*-------------------- Solver Information --------------------*/

    std::shared_ptr<LinearSolver> linear_solver;

  protected:
    /*-------------------- Problem Information --------------------*/

    unsigned int n_groups;

    /**
     * The total number of delayed neutron precursors across all materials.
     * This is defined by the total number of unique decay constants.
     */
    unsigned int n_precursors = 0;

    /**
     * The maximum number of precursors on a material. This is used to size
     * the precursor vector so that only a limited number of precursors are
     * stored per cell.
     */
    unsigned int max_precursors = 0;

    std::vector<std::shared_ptr<CrossSections>> material_xs;
    std::vector<std::shared_ptr<IsotropicMGSource>> material_src;
    std::vector<LightWeightCrossSections> cellwise_xs;

    /**
     * Map a material ID to particular cross-sections. This mapping alleviates
     * the need to store multiple copies of cross-sections when used on
     * several materials.
     */
    std::vector<unsigned int> matid_to_xs_map;

    /**
     * Map a material ID to a particular isotropic multi-group source object.
     * This mapping alleviates the need to store multiple copies of an
     * isotropic multi-group source object when used on several materials.
     */
    std::vector<unsigned int> matid_to_src_map;

    /**
     * The multi-group boundary conditions. This is a vector of vectors of
     * pointers to boundary conditions. The outer indexing corresponds to the
     * boundary index and the inner index to the group. These are created at
     * solver initialization.
     */
    std::vector<std::vector<BndryPtr>> boundaries;

    /**
     * The multi-group scalar flux vector.
     *
     * This vector stores group values contiguously. This means that for each
     * node, all groups for that node are contiguous.
     */
    Vector phi;
    Vector phi_ell;  ///< The multi-group scalar flux last iteration.

    /**
     * The cell-wise delayed neutron precursor concentration.
     *
     * This vector stores concentrations contiguously. This means that for
     * each cell, all concentrations for that node are contiguous. Further,
     * to eliminate wasteful storage in problems with many materials and
     * precursors, only enough data to hold the maximum number of precursors
     * on a material is allocated per cell. To obtain particular precursors
     * from particular materials, knowledge of the cell material IDs must be
     * incorporated.
     */
    Vector precursors;

    SparseMatrix A;  ///< The multi-group matrix.
    Vector b; ///< The right-hand side vector.

  public:
    /*-------------------- Public Routines --------------------*/

    /** Initialize the multi-group diffusion solver. */
    virtual void initialize();

    /** Execute the multi-group diffusion solver. */
    virtual void execute();

    /**
     * Write the result of the simulation to an output file with the
     * specified prefix. This creates a file named <tt><file_prefix>.data</tt>
     * in the \p output_directory
     */
    virtual void
    write(const std::string output_directory,
          const std::string file_prefix) const;

  protected:
    /*-------------------- Initialization Routines --------------------*/

    void initialize_materials();
    void initialize_boundaries();

    /*-------------------- Solve Routines --------------------*/

    /**
     * Lag the specified sources within \p source_flags and iteratively solve
     * the multi-group system.
     *
     * Each iteration, this routine adds the sources specified by \p
     * source_flags to the right-hand side and solves the multi-group system.
     * In many cases, both scattering and fission are lagged which results
     * in an symmetric positive definite linear system.
     *
     * The number of iterations and the final convergence check value are
     * returned as a pair.
     */
    std::pair<unsigned int, double>
    iterative_solve(SourceFlags source_flags);

    /**
     * Compute the steady-state precursor concentration profile.
     *
     * In the steady-state formulation, the delayed neutron precursor
     * concentrations can be written as a function of only the multi-group
     * scalar flux and substituted. This eliminates the precursor concentrations
     * from the system, resulting in this routine being a post-processing step.
     */
    void compute_precursors();

    /*-------------------- Assembly Routines --------------------*/

    /**
     * Assemble the multi-group matrix according to the specified
     * \p assembler_flags.
     *
     * By default, the within-group terms (total  interaction, buckling,
     * diffusion, and boundary) are included in the matrix. When specified,
     * the cross-group scattering and fission terms may be included.
     */
    void assemble_matrix(AssemblerFlags assembler_flags = NO_ASSEMBLER_FLAGS);

    /**
     * Accumulate sources into the right-hand side according to the specified
     * \p source_flags.
     *
     * This routine is additive, implying that if the right-hand side needs to
     * be cleared, this must be done external to this routine. The purpose of
     * the additive nature is to accommodate for different iteration levels,
     * such as in a \f$ k \f$-eigenvalue problem where the fission source is
     * held constant while the scattering source is iterated on. Available
     * sources are the inhomogeneous, scattering, fission, and boundary source.
     */
    void set_source(SourceFlags source_flags = NO_SOURCE_FLAGS);
  };

}

#endif //STEADYSTATE_SOLVER_H
