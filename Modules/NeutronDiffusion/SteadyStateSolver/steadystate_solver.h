#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "../boundaries.h"

#include "mesh.h"
#include "Discretization/discretization.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"
#include "LinearSolvers/linear_solver.h"

#include "material.h"
#include "CrossSections/cross_sections.h"
#include "lightweight_xs.h"

#include <string>


using namespace Math;


namespace NeutronDiffusion
{
  enum class Algorithm
  {
    DIRECT = 0,   ///< Solve the full multi-group system.
    ITERATIVE = 1  ///< Iterate on the cross-group terms.
  };

  //######################################################################

  /** Bitwise source flags for right-hand side vector construction. */
  enum SourceFlags : int
  {
    NO_SOURCE_FLAGS = 0,
    APPLY_MATERIAL_SOURCE = (1 << 0),
    APPLY_SCATTER_SOURCE = (1 << 1),
    APPLY_FISSION_SOURCE = (1 << 2),
    APPLY_BOUNDARY_SOURCE = (1 << 3)
  };


  inline SourceFlags operator|(const SourceFlags f1,
                               const SourceFlags f2)
  {
    return static_cast<SourceFlags>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
  }


  /** Bitwise assembler flags for matrix construction. */
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

  /** A steady-state multi-group diffusion solver. */
  class SteadyStateSolver
  {
  protected:
    typedef Grid::Mesh Mesh;
    typedef SpatialDiscretizationMethod SDMethod;

    typedef Physics::Material Material;
    typedef Physics::MaterialPropertyType MaterialPropertyType;
    typedef Physics::CrossSections CrossSections;
    typedef Physics::LightWeightCrossSections LightWeightCrossSections;
    typedef Physics::IsotropicMultiGroupSource IsotropicMGSource;

    typedef std::vector<double> RobinBndryVals;
    typedef std::shared_ptr<Boundary> BndryPtr;

    typedef LinearSolver::LinearSolverBase<SparseMatrix> LinearSolverBase;

    /*-------------------- General Information --------------------*/
  public:
    unsigned int verbosity = 0;

    Algorithm algorithm = Algorithm::DIRECT;
    SDMethod discretization_method = SDMethod::FINITE_VOLUME;

    bool use_precursors = false;

    double inner_tolerance = 1.0e-6;
    unsigned int max_inner_iterations = 100;

  protected:
    unsigned int n_groups;

    /** The total number of delayed neutron precursors across all materials. */
    unsigned int n_precursors = 0;

    /**
     * The maximum number of precursors on a material. This is used to size
     * the precursor vector so that only a limited number of precursors are
     * stored per cell.
     */
    unsigned int max_precursors = 0;

    /*-------------------- Spatial Domain --------------------*/
  public:
    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<Discretization> discretization;

    /*-------------------- Physics Informations --------------------*/
  public:
    /**
     * A list of materials, each of which must contain CrossSections and,
     * optionally, an IsotropicMultiGroupSource.
     */
     std::vector<std::shared_ptr<Material>> materials;

     /**
      * A list of the group IDs to be used in the simulation. While this may
      * seem useless, it allows for a subset of available groups within
      * a cross-section library to be considered.
      */
     std::vector<unsigned int> groups;

  protected:
    std::vector<std::shared_ptr<CrossSections>> material_xs;
    std::vector<std::shared_ptr<IsotropicMGSource>> material_src;
    std::vector<LightWeightCrossSections> cellwise_xs;

    /**
     * Map a material ID to a particular CrossSections. This mapping alleviates
     * the need to store multiple copies of CrossSections when used on
     * several materials.
     */
    std::vector<int> matid_to_xs_map;

    /**
     * Map a material ID to a particular IsotropicMultiGroupSource object.
     * This mapping alleviates the need to store multiple copies of an
     * IsotropicMultiGroupSource object when used on several materials.
     */
    std::vector<int> matid_to_src_map;

     /*-------------------- Boundary Information --------------------*/
  public:
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

  protected:
    /**
     * The multi-group boundary conditions. This is a vector of vectors of
     * pointers to Boundary objects. The outer indexing corresponds to the
     * boundary index and the inner index to the group. These are created at
     * solver initialization.
     */
    std::vector<std::vector<BndryPtr>> boundaries;

    /*-------------------- Linear Solver --------------------*/
  public:
    std::shared_ptr<LinearSolverBase> linear_solver;

    /*-------------------- System Storage --------------------*/
  public:

    Vector phi;

    /** The delayed neutron precursors defined at the cell centers. */
    Vector precursors;

  protected:

    SparseMatrix A;
    Vector b;

    /*-------------------- Public Facing Routines --------------------*/
  public:

    virtual void initialize();
    virtual void execute();

    /**
     * Write the result of the simulation to an output file.
     *
     * \param output_directory The directory where the output should be placed.
     * \param file_prefix The name of the file without a suffix. By default,
     *      the suffix \p .data will be added to this input.
     */
    virtual void write(const std::string& output_directory,
                       const std::string& file_prefix) const;

    /*-------------------- Initialization Routines --------------------*/
  protected:
    /**
     * Parse the CrossSections and IsotropicMultiGroupSource objects from the
     * Material list and set the appropriate internal data.
     */
    void initialize_materials();

    /**
     * Parse the boundary specification and define the boundary conditions
     * used in the simulation. This initializes \p n_groups Boundary objects
     * for each of the spatial domain boundaries.
     */
    void initialize_boundaries();


    /*-------------------- Solve Routines --------------------*/
  protected:
    /**
     * Solve the system iteratively by iterating on the specified SourceFlags.
     * The primary utility of this for k-eigenvalue solvers where the fission
     * source is held constant (outer iterations) while the scattering source
     * is converged (inner iterations).
     *
     * \param source_flags Bitwise flags defining the source terms to iterate
     *      on and converge.
     */
    void iterative_solve(SourceFlags source_flags);

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
    void assemble_matrix(AssemblerFlags assembler_flags = NO_ASSEMBLER_FLAGS);

    /**
     * Set the right-hand side source vector. This is an additive routine which
     * will only add the specified sources to the source vector. Source options
     * include the inhomogeneous source, scattering source, fission source,
     * and boundary source
     *
     * \param source_flags Bitwise flags used to specify which sources are
     *      added to the source vector.
     */
    void set_source(SourceFlags source_flags = NO_SOURCE_FLAGS);

    /** Compute the steady-state precursor concentration profile. */
    void compute_precursors();
  };

}

#endif //STEADYSTATE_SOLVER_H
