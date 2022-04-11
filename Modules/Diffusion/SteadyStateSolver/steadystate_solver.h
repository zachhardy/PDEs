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

/**
 * \brief A steady-state diffusion solver.
 *
 * The intended application for this solver is multigroup neutron diffusion
 * problems.
 */
class SteadyStateSolver
{
public:
  typedef std::vector<double> BoundaryValues;
  typedef std::vector<BoundaryValues> MultiGroupBoundaryValues;
  typedef std::pair<BoundaryType, MultiGroupBoundaryValues> BoundaryDescription;
  typedef std::shared_ptr<Boundary> BndryPtr;

private:
  /// The solver type as a string.
  std::string solver_string = "diffusion::SteadyStateSolver";

public:
  //========== GENERAL INFORMATION ==========

  /// The number of energy groups.
  int n_groups = 0;
  /// The number of delayed neutron precursors.
  int n_precursors = 0;
  /// The maximumum number of delayed neutron precursors on a single material.
  int max_precursors_per_material = 0;
  /// Flag for using delayed neutron precursors.
  bool use_precursors = false;

  //========== SPATIAL DOMAIN INFORMATION ==========

  /// A pointer to the Mesh the problem is defined on.
  std::shared_ptr<Mesh>  mesh;

  /// A pointer to the discretization of the problem on Mesh.
  std::shared_ptr<SpatialDiscretization> discretization;

  //========== MATERIAL INFORMATION ==========

  /// The Material objects that exist within the problem.
  std::vector<std::shared_ptr<Material>> materials;

  /// Shorthand accessors for CrossSections properties.
  std::vector<std::shared_ptr<CrossSections>>  material_xs;

  /// Shorthand accessors for IsotropicMultiGroupSource properties.
  std::vector<std::shared_ptr<IsotropicMultiGroupSource>> material_src;

  /// Mapping from Cell `material_id` to CrossSections property index.
  std::vector<int> matid_to_xs_map;

  /// Mapping from Cell `material_id` to IsotropicMultiGroupSource property index.
  std::vector<int> matid_to_src_map;

  //========== BOUNDARY CONDITIONS ==========
  /**
   * \brief A mapping between boundary ID and the boundary descriptions.
   *
   * The boundary description is given by a pair with a BoundaryType and a
   * MultiGroupBoundaryValues entry. This latter is given by a vector of
   * vectors whose outer indexing corresponds to a group and inner to the
   * scalar boundary values. The inner indexing will generally be a one-entry
   * vector unless a fully specified Robin boundary is provided.
   */
  std::vector<BoundaryDescription> boundary_info;
  /**
   * \brief The multigroup boundary conditions.
   *
   * This is a vector of vectors of pointers to Boundary objects. The outer
   * indexing corresponds to the boundary ID and the inner index to the group
   * number. These are initialized at solver initialization.
   */
  std::vector<std::vector<BndryPtr>> boundaries;


  //========== SYSTEM STORAGE ==========
  Vector phi;        ///< The scalar flux unknowns.
  Vector precursors; ///< The delayed neutron precursor unknwons.

  Vector system_rhs;
  Matrix system_matrix;

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
