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
 * \brief A steady-state multigroup neutron diffusion solver.
 *
 * \property max_precursors_per_material
 *           The maximum number of precursors that live on a material. This is
 *           used to promote sparsity in the precursor vector for problems with
 *           many different materials and precursor sets, such as in burnup
 *           applications.
 * \property matid_to_xs_map
 *           Map a material ID to a particular CrossSection object.
 *           See \ref initialize_materials implementation for details.
 * \property matid_to_src_map
 *           Map a material ID to a particular IsotropicMultiGroupSource object.
 *           See \ref initialize_materials implementation for details.
 * \property boundary_info
 *           A mapping between the boundary ID and its description.
 *           The boundary description is defined by a pair with a BoundaryType
 *           and a MultiGroupBoundaryValues entry. This latter is a vector of
 *           vectors whose outer indexing corresponds to a group and inner to
 *           the scalar boundary values. The inner indexing will generally be a
 *           one-entry vector unless a fully specified Robin boundary is
 *           provided.
 * \property boundaries
 *           The multigroup boundary conditions. This is a vector of vectors of
 *           pointers to Boundary objects. The outer indexing corresponds to the
 *           boundary ID and the inner index to the group number. These are
 *           created at solver initialization.
 */
class SteadyStateSolver
{
public:
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
  std::string solver_string = "diffusion::SteadyStateSolver";

public:
  size_t n_groups = 0; ///< The number of energy groups.
  size_t n_precursors = 0; ///< The total number of delayed neutron precursors.
  bool use_precursors = false; ///< A flag for using precursors.
  size_t max_precursors_per_material = 0;

  std::shared_ptr<grid::Mesh>  mesh;
  std::shared_ptr<SpatialDiscretization> discretization;
  std::vector<std::shared_ptr<Material>> materials;

  std::vector<std::shared_ptr<CrossSections>>  material_xs;
  std::vector<std::shared_ptr<IsotropicMGSource>> material_src;
  std::vector<int> matid_to_xs_map;
  std::vector<int> matid_to_src_map;

  std::vector<BoundaryDescription> boundary_info;
  std::vector<std::vector<BndryPtr>> boundaries;

  math::Vector phi;
  math::Vector precursors;

  math::Vector system_rhs;
  math::Matrix system_matrix;

public:

  void initialize();
  void execute();

private:
  void assemble_matrix();
  void assemble_rhs_vector();

  void fv_assemble_matrix();
  void fv_assemble_rhs_vector();

private:
  void check_inputs();
  void initialize_materials();
  void initialize_boundaries();
};

}

#endif //STEADYSTATE_SOLVER_H
