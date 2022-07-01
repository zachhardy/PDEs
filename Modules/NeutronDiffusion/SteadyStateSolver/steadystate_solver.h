#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "../boundaries.h"
#include "../groupset.h"

#include "mesh.h"
#include "Discretization/discretization.h"

#include "vector.h"
#include "LinearSolvers/linear_solver.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include <string>


using namespace Math;


namespace NeutronDiffusion
{

  /** Algorithms to solve the multigroup diffusion problem. */
  enum class SolutionTechnique
  {
    FULL_SYSTEM = 0,   ///< Solve the full multigroup system.
    GROUPSET_WISE = 1  ///< Iteratively solve by groupset.
  };

  //######################################################################

  /** Bitwise source flags. */
  enum SourceFlags : int
  {
    NO_SOURCE_FLAGS = 0,
    APPLY_MATERIAL_SOURCE = (1 << 0),
    APPLY_WGS_SCATTER_SOURCE = (1 << 1),
    APPLY_AGS_SCATTER_SOURCE = (1 << 2),
    APPLY_WGS_FISSION_SOURCE = (1 << 3),
    APPLY_AGS_FISSION_SOURCE = (1 << 4)
  };


  inline SourceFlags operator|(const SourceFlags f1,
                               const SourceFlags f2)
  {
    return static_cast<SourceFlags>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
  }


  /** Bitwise assembler flags */
  enum AssemblerFlags : int
  {
    NO_ASSEMBLER_FLAGS = 0,
    ASSEMBLE_SCATTER = (1 << 0),
    ASSEMBLE_FISSION = (1 << 1)
  };


  inline AssemblerFlags operator|(const AssemblerFlags f1,
                                  const AssemblerFlags f2)
  {
    return static_cast<AssemblerFlags>(static_cast<int>(f1) |
                                       static_cast<int>(f2));
  }

  //######################################################################

  /** A steady state solver for multigroup neutron diffusion applications. */
  class SteadyStateSolver
  {
  protected:
    typedef Grid::Mesh Mesh;

    typedef Physics::Material Material;
    typedef Physics::MaterialPropertyType MaterialPropertyType;
    typedef Physics::CrossSections CrossSections;
    typedef Physics::IsotropicMultiGroupSource IsotropicMGSource;

    typedef std::vector<double> RobinBndryVals;
    typedef std::shared_ptr<Boundary> BndryPtr;

    typedef LinearSolver::LinearSolverBase<SparseMatrix> LinearSolverBase;

  public:

    /*-------------------- General Information --------------------*/

    size_t n_groups = 0;
    size_t n_precursors = 0;

    bool use_precursors = false;

    size_t verbosity = 0;

    /*-------------------- Solver Information --------------------*/

    SolutionTechnique solution_technique = SolutionTechnique::GROUPSET_WISE;

    std::shared_ptr<LinearSolverBase> linear_solver;

    /*-------------------- Groupsets and Groups --------------------*/

    std::vector<size_t> groups;
    std::vector<Groupset> groupsets;

    /*-------------------- Spatial Grid Information --------------------*/

    std::shared_ptr<Mesh> mesh;

    DiscretizationMethod discretization_method = DiscretizationMethod::FINITE_VOLUME;
    std::shared_ptr<Discretization> discretization;

    /*-------------------- Physics Information --------------------*/

    std::vector<std::shared_ptr<Material>> materials;
    std::vector<std::shared_ptr<CrossSections>> material_xs;
    std::vector<std::shared_ptr<IsotropicMGSource>> material_src;

  protected:
    /**
     * The maximum number of precursors that live on a material.
     * This is used to promote sparsity in the precursor vector for problems with
     * many materials and precursor sets, such as in burn-up applications.
     */
    size_t max_precursors = 0;

    /**
     * Map a material ID to a particular CrossSection object.
     * This mapping alleviates the need to store multiple copies of the
     * CrossSections objects when one property appears more than once.
     */
    std::vector<int> matid_to_xs_map;

    /**
     * Map a material ID to a particular IsotropicMultiGroupSource object.
     * This mapping alleviates the need to store multiple copies of the
     * IsotropicMultiGroupSource objects when one property appears more than
     * once.
     */
    std::vector<int> matid_to_src_map;

  public:

    /*-------------------- Boundary Information --------------------*/

    /**
     * A list containing a pair with the boundary type and index corresponding
     * to the location of the boundary values within the boundary values vector.
     * This is similar to the matid_to_xs_map attribute.
     */
    std::vector<std::pair<BoundaryType, size_t>> boundary_info;

    /**
     * The multigroup boundary values. The outer index corresponds to the
     * boundary index, the middle to the group, and the last to the boundary
     * value index. For non-Robin boundaries, this always has one entry at the
     * innermost lever. For Robin boundaries, three entries in the order of
     * <tt>(a, b, f)</tt> are used.
     */
    std::vector<std::vector<std::vector<double>>> boundary_values;

  protected:

    /**
     * The multigroup boundary conditions. This is a vector of vectors of
     * pointers to Boundary objects. The outer indexing corresponds to the
     * boundary index and the inner index to the group. These are created at
     * solver initialization.
     */
    std::vector<std::vector<BndryPtr>> boundaries;

  public:

    /*-------------------- Solutions --------------------*/

    Vector phi;
    Vector phi_ell;
    Vector precursors;

  public:

    /*-------------------- Public Facing Routines --------------------*/

    /**
     * Initialize the solver. This routine ensures that the specified setup is
     * valid. For example, the mesh, groups, groupsets, materials, and
     * boundaries are all checked and relevant properties are initialized.
     */
    virtual void initialize();

    /** Run the steady state multigroup diffusion simulation. */
    virtual void execute();

  protected:
    /*-------------------- Solve Routines --------------------*/

    /** \name Solve Routines */
    // @{

    /**
     * Solve the system for a single groupset. This routine iteratively converges
     * the scattering and fission term. In this particular case, matrices will
     * uniformly be SPD.
     */
    void solve_groupset(Groupset& groupset, SourceFlags source_flags);

    /**
     * Solve the full multigroup system. This routine uses a direct solver to
     * solve the <tt>often</tt> asymmetric system.
     */
    void solve_full_system(SourceFlags source_flags);

    // @}

    /*-------------------- Assembly Routines --------------------*/

    /** \name Assembly Routines */
    // @{

    /**
     * Assemble the matrix for the specified \p groupset.
     *
     * If solving the full system, this routine assembles the full multigroup
     * operator for the groupset, including off-diagonal scattering and fission
     * coupling terms. Otherwise, this routine assembles the within-group system
     * for the groups within the groupset.
     *
     * \param groupset The groupset to construct the matrix for.
     */
    void
    assemble_matrix(Groupset& groupset,
                    AssemblerFlags assembler_flags = NO_ASSEMBLER_FLAGS);

    /**
     * Set the right-hand side source vector for the specified groupset.
     *
     * The available source flags include the material source, across groupset
     * scattering and fission terms, and within groupset scattering and fission
     * terms.
     *
     * \param groupset The groupset to construct the source for.
     * \param source_flags The terms to assemble into the destination vector.
     */
    void set_source(Groupset& groupset, SourceFlags source_flags);

    /** Compute the steady-state delayed neutron precursor concentrations. */
    void compute_precursors();

    //@}

    /*-------------------- Initialization Routines --------------------*/

    /** \name Initialization Routines */
    // @{

    /**
     * Validate the general setup of the simulation.
     */
    void input_checks();

    /**
     * Parse the appropriate properties from the materials list, validate
     * compatibility with other properties and the simulation, and set relevant
     * attributes of the simulation obtained from materials.
     */
    void initialize_materials();

    /**
     * Create a boundary condition for each boundary and each group.
     */
    void initialize_boundaries();

    // @}

    /*-------------------- Vector Operations --------------------*/

    /** \name Vector Operations */
    // @{

    /**
     * Transfer a groupset vector to a full multigroup vector.
     *
     * \param groupset The groupset the vector belongs to.
     * \param x The groupset vector to be transferred.
     * \param dst The destination multigroup vector.
     */
    void scoped_transfer(const Groupset& groupset,
                         const Vector& x, Vector& dst);

    /**
     * Copy the elements corresponding to the specified groupset from one full
     * multigroup vector to another.
     *
     * \param groupset The groupset to copy data from.
     * \param x The multigroup vector to be copied.
     * \param dst The destination multigroup vector.
     */
    void scoped_copy(const Groupset& groupset,
                     const Vector& x, Vector& dst);

    /**
     * Return the \f$\ell_2\f$-norm between the last two iterates of the
     * multigroup flux vectors for elements belonging to the specified groupset.
     *
     * \param groupset The groupset to compute the change within.
     */
    double compute_change(const Groupset& groupset);

    // @}

    /*-------------------- File I/O --------------------*/

    /** Write the solution data to a binary file. */
    void write_snapshot(const std::string& file_base) const;
  };

}

#endif //STEADYSTATE_SOLVER_H
