#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "../diffusion.h"
#include "../boundaries.h"

#include "Mesh/mesh.h"
#include "SpatialDiscretization/spatial_discretization.h"

#include "Material/material.h"
#include "Material/CrossSections/cross_sections.h"

#include "Math/vector.h"
#include "Math/matrix.h"


namespace diffusion
{

/// A steady-state multigroup neutron diffusion solver.
class SteadyStateSolver
{
public:
  typedef std::vector<double> BoundaryValues;
  typedef std::vector<BoundaryValues> MultiGroupBoundaryValues;
  typedef std::pair<BoundaryType, MultiGroupBoundaryValues> BoundaryDescription;
  typedef std::shared_ptr<Boundary> BndryPtr;

private:
  std::string solver_string = "diffusion::SteadyStateSolver";

public:
  size_t n_groups = 0; ///< The number of energy groups.
  size_t n_precursors = 0; ///< The total number of delayed neutron precursors.
  bool use_precursors = false; ///< A flag for using precursors.

  /** The maximum number of precursors that live on any given material.
   *  This is used to promote sparsity in the precursor vector for problems
   *  with many materials with different precursor properties, such as burnup
   *  applications. */
  size_t max_precursors_per_material = 0;

  std::shared_ptr<grid::Mesh>  mesh;
  std::shared_ptr<math::SpatialDiscretization> discretization;
  std::vector<std::shared_ptr<material::Material>> materials;

  std::vector<std::shared_ptr<material::CrossSections>>  material_xs;
  std::vector<std::shared_ptr<material::IsotropicMultiGroupSource>> material_src;

  /** Map a material ID to a particular CrossSection object.
   *  See \ref initialize_materials. */
  std::vector<int> matid_to_xs_map;
  /** Map a material ID to a particular IsotropicMultiGroupSource object.
   * See \ref initialize_materials. */
  std::vector<int> matid_to_src_map;

  /** A mapping between the boundary ID and its description.
   *  The boundary description is defined by a pair with a BoundaryType and a
   *  MultiGroupBoundaryValues entry. This latter is a vector of
   *  vectors whose outer indexing corresponds to a group and inner to the
   *  scalar boundary values. The inner indexing will generally be a one-entry
   *  vector unless a fully specified Robin boundary is provided. */
  std::vector<BoundaryDescription> boundary_info;

  /** The multigroup boundary conditions. This is a vector of vectors of pointers
   *  to Boundary objects. The outer indexing corresponds to the boundary ID and
   *  the inner index to the group number. These are initialized at solver
   *  initialization. */
  std::vector<std::vector<BndryPtr>> boundaries;

  math::Vector phi;
  math::Vector precursors;

  math::Vector system_rhs;
  math::Matrix system_matrix;

public:

  void initialize(); ///< Initialize the diffusion solver.
  void execute(); ///< Execute the diffusion solver.

  void assemble_matrix(); ///< Assemble the matrix
  void assemble_rhs_vector(); ///< Assemble the right-hand side.

private:
  void fv_assemble_matrix(); ///< Assemble the matrix with the FV method.
  void fv_assemble_rhs_vector(); ///< Assemble the RHS with the FV method.

private:

  void check_inputs(); ///< Check the solver inputs for validity.
  void initialize_materials(); ///< Initialize the materials for simple access.
  void initialize_boundaries(); ///< Initialize the boundary conditions.
};

}

#endif //STEADYSTATE_SOLVER_H
