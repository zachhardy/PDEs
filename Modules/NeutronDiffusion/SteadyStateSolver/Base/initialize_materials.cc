#include "steadystate_solver.h"

#include <set>

/**
 * Grab the appropriate material properties from the materials list.
 *
 * This routine performs checks to ensure that the number of materials matches
 * the number of unique material identifiers on the mesh, that CrossSections
 * objects exist on each Material::properties list, and that the group
 * structures among the CrossSections and IsotropicMultiGroupSource objects are
 * compatible with the specified group structure. This routine also defines
 * a mapping between unique material IDs and the corresponding property's
 * location in the associated list. Lastly, the number of groups and precursors
 * are set. The number of groups is simply the size of the \p groups vector and
 * the number of precursors is the number of unique decay constants across all
 * materials.
 */
void NeutronDiffusion::SteadyStateSolver::initialize_materials()
{
  std::cout << "Initializing materials...\n";

  // Determine the unique material IDs on the mesh
  std::set<int> unique_material_ids;
  std::vector<size_t> invalid_cells;
  for (const auto& cell : mesh->cells)
  {
    unique_material_ids.insert(cell.material_id);
    if (cell.material_id < 0)
      invalid_cells.emplace_back(cell.id);
  }

  // If all cells are invalid, set to zero
  if (invalid_cells.size() == mesh->cells.size())
  {
    unique_material_ids.clear();
    unique_material_ids.insert(0);
    for (auto& cell : mesh->cells)
      cell.material_id = 0;
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

  //================================================== Loop over materials
  for (const int& material_id : unique_material_ids)
  {
    auto material = materials[material_id];

    //================================================== Loop over properties
    bool found_xs = false;
    for (const auto& property : material->properties)
    {
      //============================== Parse cross sections
      if (property->type == MaterialPropertyType::CROSS_SECTIONS)
      {
        auto xs = std::static_pointer_cast<CrossSections>(property);

        if (xs->n_groups < groups.size())
        {
          std::stringstream err;
          err << solver_string << __FUNCTION__ << ": "
              << "Cross sections encountered with fewer groups ("
              << xs->n_groups << ") than the simulation (" << groups.size()
              << "). All cross sections must have at least as many groups as "
              << "the simulation.";
          throw std::runtime_error(err.str());
        }

        material_xs.emplace_back(xs);
        matid_to_xs_map[material_id] = material_xs.size() - 1;
        found_xs = true;
      }//if CrossSections

      //============================== Parse multigroup sources
      else if (property->type == MaterialPropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto src = std::static_pointer_cast<IsotropicMGSource>(property);

        if (src->values.size() < groups.size())
        {
          std::stringstream err;
          err << solver_string << __FUNCTION__ << ": "
              << "Multigtoup source encountered with fewer groups ("
              << src->values.size() << ") than the simulation ("
              << groups.size() << "). All sources must have at least as many "
              << "groups as the simulation.";
          throw std::runtime_error(err.str());
        }

        material_src.emplace_back(src);
        matid_to_src_map[material_id] = material_src.size() - 1;
      }//if IsotropicMultiGroupSource
    }//for properties

    // Make sure cross sections were found
    if (not found_xs)
    {
      std::stringstream err;
      err << solver_string << __FUNCTION__ << ": "
          << "Cross sections not found on material " << material_id << ".";
      throw std::runtime_error(err.str());
    }
  }//for materials

  // Define the number of groups
  n_groups = groups.size();

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
  }//if use precursors

  std::cout << "Materials initialized: " << materials.size() << "\n";
  std::cout << "Done initializing materials.\n";
}