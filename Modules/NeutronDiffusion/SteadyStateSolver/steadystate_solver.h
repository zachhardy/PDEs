#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "../boundaries.h"

#include "Grid/mesh.h"
#include "Discretization/discretization.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "vector.h"
#include "matrix.h"

#include "LinearSolvers/linear_solver.h"


namespace neutron_diffusion
{

/// A steady state solver for multigroup neutron diffusion applications.
class SteadyStateSolver
{
private:
  const std::string solver_string = "diffusion::SteadyStateSolver::";

protected:
  typedef grid::Mesh Mesh;
  typedef math::DiscretizationMethod DiscretizationMethod;
  typedef math::Discretization Discretization;

  typedef physics::Material Material;
  typedef physics::MaterialPropertyType MaterialPropertyType;
  typedef physics::CrossSections CrossSections;
  typedef physics::IsotropicMultiGroupSource IsotropicMGSource;

  typedef std::vector<double> RobinBndryVals;
  typedef std::shared_ptr<Boundary> BndryPtr;

  typedef math::LinearSolverType LinearSolverType;
  typedef math::LinearSolver LinearSolver;

public:

  /*---------- Options ----------*/

  LinearSolverType linear_solver_type;

  /*---------- General information ----------*/

  size_t n_groups = 0;
  size_t n_precursors = 0;
  bool use_precursors = false;

  /*---------- Spatial grid information ----------*/

  std::shared_ptr<Mesh> mesh;
  std::shared_ptr<Discretization> discretization;

  /*---------- Material information ----------*/

public:
  std::vector<std::shared_ptr<Material>> materials;
  std::vector<std::shared_ptr<CrossSections>>  material_xs;
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

  /*---------- Boundary information ----------*/
protected:
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

  /** The multigroup boundary conditions. This is a vector of vectors of
   *  pointers to Boundary objects. The outer indexing corresponds to the
   *  boundary index and the inner index to the group. These are created at
   *  solver initialization. */
  std::vector<std::vector<BndryPtr>> boundaries;

public:
  void initialize();

  void add_boundary(BoundaryType bndry_type,
                    std::vector<std::vector<double>> mg_bndry_values = {});

protected:
  void input_checks();
  void initialize_materials();
  void initialize_boundaries();

};

}

#endif //STEADYSTATE_SOLVER_H
