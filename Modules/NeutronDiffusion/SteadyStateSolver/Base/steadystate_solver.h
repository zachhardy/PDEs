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


namespace NeutronDiffusion
{

/// Algorithms to solve the multigroup diffusion problem.
enum class SolutionTechnique
{
  FULL_SYSTEM = 0,   ///< Solve the full multigroup system.
  GROUPSET_WISE = 1  ///< Iteratively solve by groupset.
};

//######################################################################

/// Bitwise source flags
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

/// A steady state solver for multigroup neutron diffusion applications.
class SteadyStateSolver
{
protected:
  const std::string solver_string = "diffusion::SteadyStateSolver::";

protected:
  typedef Grid::Mesh Mesh;
  typedef Math::DiscretizationMethod DiscretizationMethod;
  typedef Math::Discretization Discretization;

  typedef Physics::Material Material;
  typedef Physics::MaterialPropertyType MaterialPropertyType;
  typedef Physics::CrossSections CrossSections;
  typedef Physics::IsotropicMultiGroupSource IsotropicMGSource;

  typedef std::vector<double> RobinBndryVals;
  typedef std::shared_ptr<Boundary> BndryPtr;

public:

  /*---------- Options ----------*/
  SolutionTechnique solution_technique = SolutionTechnique::GROUPSET_WISE;

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

  /*---------- Material Information ----------*/
public:
  std::vector<std::shared_ptr<Material>> materials;
  std::vector<std::shared_ptr<CrossSections>> material_xs;
  std::vector<std::shared_ptr<IsotropicMGSource>> material_src;

protected:
  /** The maximum number of precursors that live on a material.
   *  This is used to promote sparsity in the precursor vector for problems with
   *  many different materials and precursor sets, such as in burnup
   *  applications. */
  size_t max_precursors_per_material = 0;

  /** Map a material ID to a particular CrossSection object.
   *  This mapping alleviates the need to store multiple copies of the
   *  CrossSections objects when one property appears more than once. */
  std::vector<int> matid_to_xs_map;

  /** Map a material ID to a particular IsotropicMultiGroupSource object.
   *  This mapping alleviates the need to store multiple copies of the
   *  IsotropicMultiGroupSource objects when one property appears more than
   *  once. */
  std::vector<int> matid_to_src_map;

  /*---------- Boundary Information ----------*/
public:
  /** A list containing a pair with the boundary type and index corresponding
   *  to the location of the boundary values within the boundary values vector.
   *  This is similar to the matid_to_xs_map attribute. */
  std::vector<std::pair<BoundaryType, size_t>> boundary_info;

  /** The multigroup boundary values. The outer index corresponds to the
   *  boundary index, the middle to the group, and the last to the boundary
   *  value index. For non-Robin boundaries, this always has one entry at the
   *  innermost lever. For Robin boundaries, three entries in the order of
   *  <tt>(a, b, f)</tt> are used. */
  std::vector<std::vector<std::vector<double>>> boundary_values;

protected:
  /** The multigroup boundary conditions. This is a vector of vectors of
   *  pointers to Boundary objects. The outer indexing corresponds to the
   *  boundary index and the inner index to the group. These are created at
   *  solver initialization. */
  std::vector<std::vector<BndryPtr>> boundaries;

public:
  /*---------- Solutions ----------*/
  Math::Vector phi;
  Math::Vector phi_ell;

  Math::Vector precursors;

public:
  void
  initialize();

  void
  execute();

  void
  solve_groupset(Groupset& groupset, SourceFlags source_flags);

protected:
  /** Virtual function for assembling a groupset matrix. */
  virtual void assemble_matrix(Groupset& groupset) = 0;

  /** Virtual function for setting a groupset source. */
  virtual void set_source(Groupset& groupset, Math::Vector& b,
                          SourceFlags source_flags) = 0;

protected:
  void input_checks();

  void initialize_materials();
  void initialize_boundaries();

  /** Virtual function for creating a discretization. */
  virtual void initialize_discretization() = 0;

protected:
  virtual void scoped_transfer(const Groupset& groupset,
                               const Math::Vector& x,
                               Math::Vector& destination) = 0;

  virtual void scoped_copy(const Groupset& groupset,
                           const Math::Vector& x,
                           Math::Vector& destination) = 0;

  virtual double compute_change(const Groupset& groupset) = 0;
};

}

#endif //STEADYSTATE_SOLVER_H
