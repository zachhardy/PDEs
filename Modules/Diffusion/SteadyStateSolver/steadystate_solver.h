#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "../diffusion.h"
#include "../boundaries.h"

#include "Mesh/mesh.h"
#include "spatial_discretization.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "vector.h"
#include "matrix.h"



namespace diffusion
{

/**
 *
 * \property SteadyStateSolver::matid_to_xs_map
 * \brief Map a material ID to a particular CrossSection object.
 *
 * See \ref initialize_materials implementation for details.
 *
 * \property SteadyStateSolver::matid_to_src_map
 * \brief Map a material ID to a particular IsotropicMultiGroupSource object.
 *
 * See \ref initialize_materials implementation for details.
 *
 * \property SteadyStateSolver::boundary_info
 * \brief A mapping between the boundary ID and its description.
 *
 * The boundary description is defined by a pair with a BoundaryType
 * and a MultiGroupBoundaryValues entry. This latter is a vector of
 * vectors whose outer indexing corresponds to a group and inner to
 * the scalar boundary values. The inner indexing will generally be a
 * one-entry vector unless a fully specified Robin boundary is
 * provided.
 *
 * \property SteadyStateSolver::boundaries
 * \brief The multigroup boundary conditions.
 *
 * This is a vector of vectors of pointers to Boundary objects. The outer
 * indexing corresponds to the boundary ID and the inner index to the group
 * number. These are created at solver initialization.
 */

/// A steady state solver for multigroup neutron diffusion applications.
class SteadyStateSolver
{
protected:
  typedef std::vector<double> BoundaryValues;
  typedef std::vector<BoundaryValues> MultiGroupBoundaryValues;
  typedef std::pair<BoundaryType, MultiGroupBoundaryValues> BoundaryDescription;
  typedef std::shared_ptr<Boundary> BndryPtr;

  typedef grid::Mesh Mesh;
  typedef discretization::DiscretizationMethod DiscretizationMethod;
  typedef discretization::SpatialDiscretization SpatialDiscretization;

  typedef material::Material Material;
  typedef material::MaterialPropertyType MaterialPropertyType;
  typedef material::CrossSections CrossSections;
  typedef material::IsotropicMultiGroupSource IsotropicMGSource;

private:
  const std::string solver_string = "diffusion::SteadyStateSolver";

public:
  size_t n_groups = 0;
  size_t n_precursors = 0;
  bool use_precursors = false;

  std::shared_ptr<grid::Mesh>  mesh;
  std::shared_ptr<SpatialDiscretization> discretization;

  std::vector<std::shared_ptr<Material>> materials;
  std::vector<std::shared_ptr<CrossSections>>  material_xs;
  std::vector<std::shared_ptr<IsotropicMGSource>> material_src;

  /** A mapping between the boundary ID and its description. The boundary
   * description is defined by a pair with a BoundaryType and a
   * MultiGroupBoundaryValues entry. This latter is a vector of vectors whose
   * outer indexing corresponds to a group and inner to the scalar boundary
   * values. The inner indexing will generally be a one-entry vector unless a
   * fully specified Robin boundary is provided.
   */
  std::vector<BoundaryDescription> boundary_info;

  math::Vector phi;
  math::Vector precursors;

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

  /** The multigroup boundary conditions.
   *  This is a vector of vectors of pointers to Boundary objects. The outer
   *  indexing corresponds to the boundary ID and the inner index to the group
   *  number. These are created at solver initialization.
   */
  std::vector<std::vector<BndryPtr>> boundaries;

  math::Vector system_rhs;
  math::Matrix system_matrix;

public:
  void initialize();
  void execute();

protected:
  void assemble_matrix();
  void assemble_rhs_vector();

  void fv_assemble_matrix();
  void fv_assemble_rhs_vector();

  void check_inputs();
  void initialize_materials();
  void initialize_boundaries();
};

}

#endif //STEADYSTATE_SOLVER_H
