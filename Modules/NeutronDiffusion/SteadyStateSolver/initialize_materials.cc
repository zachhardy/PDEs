#include "steadystate_solver.h"

#include <set>
#include <cassert>


void
NeutronDiffusion::SteadyStateSolver::
initialize_materials()
{
  std::cout << "Initializing materials.\n";

  //============================================================
  // Check the material IDs on the mesh
  //============================================================

  // Determine the unique material IDs on the mesh
  std::set<unsigned int> unique_material_ids;
  std::vector<size_t> invalid_cells;
  for (const auto &cell: mesh->cells)
  {
    unique_material_ids.insert(cell.material_id);
    if (cell.material_id == -1)
      invalid_cells.emplace_back(cell.id);
  }
  assert(invalid_cells.size() == 0 ||
         invalid_cells.size() == mesh->cells.size());

  // If all cells are invalid, set to zero
  if (invalid_cells.size() == mesh->cells.size())
  {
    unique_material_ids.clear();
    unique_material_ids.insert(0);
    for (auto &cell: mesh->cells)
      cell.material_id = 0;
  }

  //============================================================
  // Clear the current material information
  //============================================================

  material_xs.clear();
  material_src.clear();

  size_t n_materials = materials.size();
  matid_to_xs_map.assign(n_materials, -1);
  matid_to_src_map.assign(n_materials, -1);

  //============================================================
  // Go through the materials and their properties to get
  // cross-sections and multi-group sources and perform checks
  // to ensure these are valid.
  //============================================================

  for (const int& material_id: unique_material_ids)
  {
    const auto material = materials[material_id];

    bool found_xs = false;
    for (const auto& property: material->properties)
    {
      // Get cross-section properties
      if (property->type() == MaterialPropertyType::CROSS_SECTIONS)
      {
        auto xs = std::static_pointer_cast<CrossSections>(property);
        if (n_groups == 0) n_groups = xs->n_groups;
        assert(xs->n_groups >= n_groups);

        material_xs.emplace_back(xs);
        matid_to_xs_map[material_id] = material_xs.size() - 1;
        found_xs = true;
      }

      // Get multigroup source properties
      else if (property->type() == MaterialPropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto src = std::static_pointer_cast<IsotropicMGSource>(property);
        assert(src->values.size() >= n_groups);

        material_src.emplace_back(src);
        matid_to_src_map[material_id] = material_src.size() - 1;
      }
    }//for properties
    assert(found_xs);
  }//for materials

  //============================================================
  // Define cell-wise cross sections
  //============================================================

  cellwise_xs.reserve(mesh->cells.size());
  for (const auto &cell: mesh->cells)
  {
    const auto xs = material_xs[matid_to_xs_map[cell.material_id]];
    cellwise_xs.emplace_back(xs);
  }

  //============================================================
  // Define precursor quantities
  //============================================================

  // Define the precursor properties
  if (use_precursors)
  {
    max_precursors = 0;
    std::set<double> unique_decay_constants;
    for (const auto& xs : material_xs)
    {
      for (unsigned int j = 0; j < xs->n_precursors; ++j)
        unique_decay_constants.insert(xs->precursor_lambda[j]);

      if (xs->n_precursors > max_precursors)
        max_precursors = xs->n_precursors;
    }
    n_precursors = unique_decay_constants.size();

    if (n_precursors == 0)
      use_precursors = false;
  }//if use precursors

  std::cout << "Materials initialized: " << materials.size() << std::endl;
}