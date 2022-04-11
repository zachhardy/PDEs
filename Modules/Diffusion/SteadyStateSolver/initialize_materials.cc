#include "steadystate_solver.h"

#include <set>


void diffusion::SteadyStateSolver::initialize_materials()
{
  // Determine unique material ids
  std::set<int> unique_material_ids;
  for (const auto& cell : mesh->cells)
    unique_material_ids.insert(cell->material_id);

  // Checks
  if (materials.size() != unique_material_ids.size())
  {
    std::stringstream err;
    err << "diffusion::SteadyStateSolver::" << __FUNCTION__ << ": "
        << "Number of materials and unique material ids must agree.";
    throw std::runtime_error(err.str());
  }

  // Clear existing properties
  material_xs.clear();
  material_src.clear();
  matid_to_xs_map.assign(materials.size(), -1);
  matid_to_src_map.assign(materials.size(), -1);

  // Iterate over unique material ids
  typedef IsotropicMultiGroupSource IsotropicMGSource;
  for (const int& material_id : unique_material_ids)
  {
    auto material = materials[material_id];

    // Iterate over properties
    bool found_xs = false; bool found_src = false;
    for (const auto& property : material->properties)
    {
      // Handle cross sections
      if (property->type == MaterialPropertyType::CROSS_SECTIONS)
      {
        auto xs = std::static_pointer_cast<CrossSections>(property);
        material_xs.push_back(xs);
        matid_to_xs_map[material_id] = material_xs.size() - 1;
        if (n_groups == 0) n_groups = xs->n_groups;
        found_xs = true;
      }

      // Handle isotropic multigroup sources
      else if (property->type == MaterialPropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto src = std::static_pointer_cast<IsotropicMGSource>(property);
        material_src.push_back(src);
        matid_to_src_map[material_id] = material_src.size() - 1;
        found_src = true;
      }
    }//for property

    // Check validity of properties
    if (!found_xs)
    {
      std::stringstream err;
      err << "diffusion::SteadyStateSolver::" << __FUNCTION__ << ": "
          << "Cross section property not found for material "
          << material_id << ".";
      throw std::runtime_error(err.str());
    }
    if (material_xs.back()->n_groups != n_groups)
    {
      std::stringstream err;
      err << "diffusion::SteadyStateSolver::" << __FUNCTION__ << ": "
          << "All cross sections must have the same number of groups. "
          << "This was violated by material " << material_id << ".";
      throw std::runtime_error(err.str());
    }
    if (found_src and material_src.back()->values.size() != n_groups)
    {
      std::stringstream err;
      err << "diffusion::SteadyStateSolver::" << __FUNCTION__ << ": "
          << "All sources must agree with the group structure. "
          << "This was violated by material " << material_id << ".";
      throw std::runtime_error(err.str());
    }
  }//for material

  // Define precursor informaton
  if (use_precursors)
  {
    n_precursors = 0;
    max_precursors_per_material = 0;
    for (const auto& xs : material_xs)
    {
      n_precursors += xs->n_precursors;
      if (xs->n_precursors > max_precursors_per_material)
        max_precursors_per_material = xs->n_precursors;
    }
  }
  if (n_precursors == 0)
    use_precursors = false;
}
