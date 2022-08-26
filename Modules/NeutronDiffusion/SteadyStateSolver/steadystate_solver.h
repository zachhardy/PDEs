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
    /**
     * Shorthand for the spatial discretization method.
     */
    using SDMethod = SpatialDiscretizationMethod;

    /**
     * Shorthand for an isotropic multi-group source.
     */
    using IsotropicMGSource = IsotropicMultiGroupSource;

    /**
     * Shorthand for Robin boundary values. Recall that generalized Robin
     * boundaries have three components.
     */
    using RobinBndryVals = std::vector<double>;

    /**
     * Shorthand for a pointer to a boundary condition.
     */
    using BndryPtr = std::shared_ptr<Boundary>;

    /**
     * Shorthand for a linear solver.
     */
    using LinearSolver = LinearSolvers::LinearSolverBase<SparseMatrix>;

    /*-------------------- General Information --------------------*/
  public:
    /**
     * The algorithm used to solve the multi-group problem.
     */
    Algorithm algorithm = Algorithm::DIRECT;

    /**
     * The spatial discretization method used.
     */
    SDMethod discretization_method = SDMethod::FINITE_VOLUME;

    /**
     * A flag for whether or not delayed neutron precursors should be used
     * or not. If delayed neutron data is present and this flag is false, the
     * total fission data derived from the prompt and delayed quantities are
     * used.
     */
    bool use_precursors = false;

    /**
     * The inner iteration convergence tolerance. For fixed source computations
     * this is somewhat of a misnomer because there are no outer iterations.
     * Outer iterations are introduced in eigenvalue computations.
     */
    double inner_tolerance = 1.0e-6;

    /**
     * The maximum number of inner iterations allowed.
     */
    unsigned int max_inner_iterations = 100;

    /**
     * The level of screen output from the solver.
     */
    unsigned int verbosity = 0;


    /*-------------------- Spatial Domain --------------------*/
    /**
     * A pointer to the mesh the problem is defined on.
     */
    std::shared_ptr<Mesh> mesh;

    /**
     * A pointer to the discretization of the mesh.
     */
    std::shared_ptr<Discretization> discretization;

    /*-------------------- Physics Information --------------------*/
    /**
     * A list of pointers to materials, each of which must contain
     * cross-sections and, optionally, an isotropic multi-group sources.
     */
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

    /**
     * A pointer to a linear solver to be used to solve the multi-group system.
     */
    std::shared_ptr<LinearSolver> linear_solver;

    /**
     * Initialize the multi-group diffusion solver.
     */
    virtual void
    initialize();

    /**
     * Execute the multi-group diffusion solver.
     */
    virtual void
    execute();

    /**
     * Write the result of the simulation to an output file.
     *
     * \param output_directory The directory where the output should be placed.
     * \param file_prefix The name of the file without a suffix. By default,
     *      the suffix \p .data will be added to this input.
     */
    virtual void
    write(const std::string& output_directory,
          const std::string& file_prefix) const;

  protected:
    /*-------------------- Problem Information --------------------*/
    /**
     * The number of energy groups in the simulation.
     */
    unsigned int n_groups;

    /**
     * The total number of delayed neutron precursors across all materials.
     */
    unsigned int n_precursors = 0;

    /**
     * The maximum number of precursors on a material. This is used to size
     * the precursor vector so that only a limited number of precursors are
     * stored per cell.
     */
    unsigned int max_precursors = 0;

    /**
     * A list of the cross-sections used in the problem. These are obtained from
     * the materials list during initialization so that the properties need not
     * be searched during execution.
     */
    std::vector<std::shared_ptr<CrossSections>> material_xs;

    /**
     * A list of the isotropic multi-group sources used in the problem. These
     * are obtained from the materials list during initialization so that the
     * properties need not be searched during execution.
     */
    std::vector<std::shared_ptr<IsotropicMGSource>> material_src;

    /**
     * A list of cell-wise light-weight cross-sections used for functional
     * cross-sections.
     */
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

    /*-------------------- Linear System Storage --------------------*/

    /**
     * The multi-group flux solution vector. This vector is group-contiguous,
     * or \f$ \phi = [\phi_{i_1, g_1}, \ldots, \phi_{i_1, g_G}, \ldots,
     * \phi_{i_N, g_1}, \ldots, \phi_{i_N, g_G}] \f$, where \f$ i_j \f$ is node
     * \f$ j \f$ and \f$ g_i \f$ is group \f$ i \f$.
     */
    Vector phi;

    /**
     * The last iterate of the multi-group flux solution vector.
     */
    Vector phi_ell;

    /**
     * The delayed neutron precursor solution vector. This vector is
     * precursor-contiguous with up to the maximum number of precursors on
     * a material. See #phi and #max_precursors.
     */
    Vector precursors;

    /**
     * The sparse matrix for the multi-group system.
     */
    SparseMatrix A;

    /**
     * The right-hand side source vector for the multi-group system.
     */
    Vector b;

    /*-------------------- Initialization Routines --------------------*/

    /**
     * Parse the cross-sections and isotropic multi-group sources objects from
     * the materials and set the appropriate internal data.
     */
    void
    initialize_materials();

    /**
     * Parse the boundary specification and define the boundary conditions
     * used in the simulation. This initializes \p n_groups Boundary objects
     * for each of the spatial domain boundaries.
     */
    void
    initialize_boundaries();

    /*-------------------- Solve Routines --------------------*/

    /**
     * Solve the system iteratively by iterating on the specified SourceFlags.
     * The primary utility of this for k-eigenvalue solvers where the fission
     * source is held constant (outer iterations) while the scattering source
     * is converged (inner iterations).
     *
     * \param source_flags Bitwise flags defining the source terms to iterate
     *      on and converge.
     */
    unsigned int
    iterative_solve(SourceFlags source_flags);

    /**
     * Compute the steady-state precursor concentration profile.
     */
    void
    compute_precursors();

    /*-------------------- Assembly Routines --------------------*/

    /**
     * Assemble the multi-group matrix according to the specified
     * AssemblerFlags. By default, the within-group terms (total interaction,
     * buckling, diffusion, and boundary) are included in the matrix. When
     * specified, the cross-group scattering and fission terms may be included.
     *
     * \param assembler_flags Bitwise flags used to specify which cross-group
     *      terms to include in the matrix.
     */
    void
    assemble_matrix(AssemblerFlags assembler_flags = NO_ASSEMBLER_FLAGS);

    /**
     * Set the right-hand side source vector. This is an additive routine which
     * will only add the specified sources to the source vector. Source options
     * include the inhomogeneous source, scattering source, fission source,
     * and boundary source
     *
     * \param source_flags Bitwise flags used to specify which sources are
     *      added to the source vector.
     */
    void
    set_source(SourceFlags source_flags = NO_SOURCE_FLAGS);
  };

}

#endif //STEADYSTATE_SOLVER_H
