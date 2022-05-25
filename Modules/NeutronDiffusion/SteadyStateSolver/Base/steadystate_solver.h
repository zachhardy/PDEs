#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "NeutronDiffusion/boundaries.h"
#include "NeutronDiffusion/Groupset/groupset.h"

#include "Grid/mesh.h"
#include "Discretization/discretization.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "vector.h"
#include "matrix.h"

using namespace pdes;
using namespace Math;
using namespace Grid;


namespace NeutronDiffusion
{

/**
 * Algorithms to solve the multigroup diffusion problem.
 */
enum class SolutionTechnique
{
  FULL_SYSTEM = 0,   ///< Solve the full multigroup system.
  GROUPSET_WISE = 1  ///< Iteratively solve by groupset.
};

//######################################################################

enum class LinearSolverType
{
  SPARSE_LU = 0,
  SPARSE_CHOLESKY = 1,
  JACOBI = 2,
  GAUSS_SEIDEL = 3
};

//######################################################################

/**
 * Bitwise source flags.
 */
enum SourceFlags : int
{
  NO_FLAGS = 0,
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

//######################################################################

/**
 * A steady state solver for multigroup neutron diffusion applications.
 */
class SteadyStateSolver
{
protected:
  const std::string solver_string = "diffusion::SteadyStateSolver::";

protected:
  typedef Physics::Material Material;
  typedef Physics::MaterialPropertyType MaterialPropertyType;
  typedef Physics::CrossSections CrossSections;
  typedef Physics::IsotropicMultiGroupSource IsotropicMGSource;

  typedef std::vector<double> RobinBndryVals;
  typedef std::shared_ptr<Boundary> BndryPtr;

public:

  /*---------- Options ----------*/
  SolutionTechnique solution_technique = SolutionTechnique::GROUPSET_WISE;
  LinearSolverType linear_solver_type = LinearSolverType::SPARSE_LU;

  /*---------- General Information ----------*/
  size_t n_groups = 0;
  size_t n_precursors = 0;
  bool use_precursors = false;

  /*---------- Groupsets and Groups ----------*/
  std::vector<size_t> groups;
  std::vector<Groupset> groupsets;

  /*---------- Spatial Grid Information ----------*/
  std::shared_ptr<Mesh> mesh;
  std::shared_ptr<Discretization> discretization;

  /*---------- Physics Information ----------*/
public:
  std::vector<std::shared_ptr<Material>> materials;
  std::vector<std::shared_ptr<CrossSections>> material_xs;
  std::vector<std::shared_ptr<IsotropicMGSource>> material_src;

protected:
  /**
   * The maximum number of precursors that live on a material.
   * This is used to promote sparsity in the precursor vector for problems with
   * many materials and precursor sets, such as in burn-up applications.
   */
  size_t max_precursors_per_material = 0;

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

  /*---------- Boundary Information ----------*/
public:
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
  /*---------- Solutions ----------*/
  Vector phi;
  Vector phi_ell;

  Vector precursors;

public:

  /**
   * Initialize the solver.
   *
   * This routine ensures that the specified setup is valid. For example, the
   * mesh, groups, groupsets, materials, and boundaries are all checked and
   * relevant properties are initialized.
   */
  void
  initialize();


  /**
   * Run the steady state multigroup diffusion simulation.
   */
  void
  execute();

protected:

  /**
   * Solve the system for a single groupset. This routine iteratively converges
   * the scattering and fission term. In this particular case, matrices will
   * uniformly be SPD.
   */
  void
  solve_groupset(Groupset& groupset,
                 SourceFlags source_flags);

  /**
   * Solve the full multigroup system. This routine uses a direct solver to
   * solve the <tt>often</tt> asymmetric system.
   */
  void
  solve_full_system(SourceFlags source_flags);


protected:
  /**
   * Virtual function for assembling a groupset matrix.
   */
  virtual void assemble_matrix(Groupset& groupset) = 0;

  /**
   * Virtual function for setting a groupset source.
   */
  virtual void set_source(Groupset& groupset, Vector& b,
                          SourceFlags source_flags) = 0;

protected:
  /**
   * Validate the general setup of the simulation.
   */
  void input_checks();

  /**
 * Grab the appropriate material properties from the materials list.
 *
 * This routine performs checks to ensure that the number of materials matches
 * the number of unique material identifiers on the mesh, that CrossSections
 * objects exist on each Physics::properties list, and that the group
 * structures among the CrossSections and IsotropicMultiGroupSource objects are
 * compatible with the specified group structure. This routine also defines
 * a mapping between unique material IDs and the corresponding property's
 * location in the associated list. Lastly, the number of groups and precursors
 * are set. The number of groups is simply the size of the \p groups vector and
 * the number of precursors is the number of unique decay constants across all
 * materials.
 */
  void initialize_materials();


  /**
   * Create a boundary condition for each boundary and each group.
   */
  void initialize_boundaries();

  /**
   * Virtual function for creating a discretization.
   */
  virtual void initialize_discretization() = 0;

protected:
  /**
   * Transfer a groupset vector to a full multigroup vector.
   */
  virtual void scoped_transfer(const Groupset& groupset,
                               const Vector& x,
                               Vector& destination) = 0;

  /**
   * Copy the elements corresponding to the specified groupset from one full
   *  multigroup vector to another.
   */
  virtual void scoped_copy(const Groupset& groupset,
                           const Vector& x,
                           Vector& destination) = 0;

  /**
   * Return the \f$\ell_2\f$-norm between the last two iterates of the
   * multigroup flux vectors for elements belonging to the specified groupset.
   */
  virtual double compute_change(const Groupset& groupset) = 0;
};

}

#endif //STEADYSTATE_SOLVER_H
