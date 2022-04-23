#include "steadystate_solver.h"

#include <set>

/**
 * \brief Grab the appropriate material properties from the materials list.
 *
 * This routine performs checks to ensure that the number of materials matches
 * the number of unique material identifiers on the mesh, that CrossSections
 * objects exist on each Material::properties list, and that the group
 * structures among the CrossSections and IsotropicMultiGroupSource objects are
 * equivalent. Lastly, this routine defines a mapping between unique material
 * IDs and the corresponding property's location in the associated lists.
 */
void neutron_diffusion::SteadyStateSolver::initialize_materials()
{
  std::cout << "Initializing materials...\n";

  // Determine the unique material IDs on the mesh
  std::set<int> unique_material_ids;
  std::vector<size_t> invalid_cells;
  for (const auto& cell : mesh->cells)
  {
    unique_material_ids.insert(cell->material_id);
    if (cell->material_id < 0)
      invalid_cells.emplace_back(cell->id);
  }

  // If all cells are invalid, set to zero
  if (invalid_cells.size() == mesh->cells.size())
  {
    unique_material_ids.clear();
    unique_material_ids.insert(0);
    for (auto& cell : mesh->cells)
      cell->material_id = 0;
  }

  // If only some cells have no material property, throw an error
  else if (invalid_cells.size() > 0)
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << invalid_cells.size() << " cells with invalid or unset material IDs "
        << "were encountered.";
    throw std::runtime_error(err.str());
  }

  // Clear the current materials
  material_xs.clear();
  material_src.clear();

  size_t n_materials = materials.size();
  matid_to_xs_map.assign(n_materials, -1);
  matid_to_src_map.assign(n_materials, -1);

  // Go through the materials and obtain the relevant properties
  for (const int& material_id : unique_material_ids)
  {
    auto material = materials[material_id];

    // Go through material properties
    bool found_xs = false;
    for (const auto& property : material->properties)
    {
      // Handle cross sections
      if (property->type == MaterialPropertyType::CROSS_SECTIONS)
      {
        auto xs = std::static_pointer_cast<CrossSections>(property);
        material_xs.emplace_back(xs);
        matid_to_xs_map[material_id] = material_xs.size() - 1;

        found_xs = true;
        if (material_xs.size() == 1)
          n_groups = xs->n_groups;
      }

      // Handle multigroup sources
      else if (property->type == MaterialPropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto src = std::static_pointer_cast<IsotropicMGSource>(property);
        material_src.emplace_back(src);
        matid_to_src_map[material_id] = material_src.size() - 1;
      }

      // Make sure cross sections were found
      if (not found_xs)
      {
        std::stringstream err;
        err << solver_string << __FUNCTION__ << ": "
            << "Cross sections not found on material " << material_id << ".";
        throw std::runtime_error(err.str());
      }

      // Check group structures for cross sections and sources.
      if (material_xs.back()->n_groups != n_groups)
      {
        std::stringstream err;
        err << solver_string << __FUNCTION__ << ": "
            << "All cross sections must have the same group structure. "
            << "This was violated on material " << material_id << ".";
        throw std::runtime_error(err.str());
      }
      if (matid_to_src_map[material_id] > 0)
      {
        if (material_src.back()->values.size() != n_groups)
        {
          std::stringstream err;
          err << solver_string << __FUNCTION__ << ": "
              << "Multigroup sources must have the same number of entries as "
              << "the number of groups. This was violated on material"
              << material_id << ".";
          throw std::runtime_error(err.str());
        }
      }
    }

    /* Define the precursor properties. For the total number of precursors
     * tally the number of unique decay constants. */
    if (use_precursors)
    {
      max_precursors_per_material = 0;
      std::set<double> unique_decay_constants;
      for (const auto& xs : material_xs)
      {
        for (size_t j = 0; j < xs->n_precursors; ++j)
          unique_decay_constants.insert(xs->precursor_lambda[j]);

        if (xs->n_precursors > max_precursors_per_material)
          max_precursors_per_material = xs->n_precursors;
      }
      n_precursors = unique_decay_constants.size();

      if (n_precursors == 0)
        use_precursors = false;
    }
  }

  std::cout << "Materials Details:\n"
            << "\t# of Material:     " << materials.size() << "\n"
            << "\t# of Groups:       " << n_groups << "\n"
            << "\t# of Precursors:   " << n_precursors << "\n";
}