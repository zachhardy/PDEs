#include "steadystate_solver.h"
#include "macros.h"

#include <set>


void
NeutronDiffusion::SteadyStateSolver::initialize_materials()
{
  std::cout << "Initializing materials.\n";

  // Determine the unique material IDs on the mesh
  std::set<int> unique_material_ids;
  std::vector<size_t> invalid_cells;
  for (const auto& cell : mesh->cells)
  {
    unique_material_ids.insert(cell.material_id);
    if (cell.material_id < 0)
      invalid_cells.emplace_back(cell.id);
  }
  Assert(invalid_cells.size() == 0 ||
         invalid_cells.size() == mesh->cells.size(),
         "Cell's with invalid or unset material IDs encounterd.");

  // If all cells are invalid, set to zero
  if (invalid_cells.size() == mesh->cells.size())
  {
    unique_material_ids.clear();
    unique_material_ids.insert(0);
    for (auto& cell : mesh->cells)
      cell.material_id = 0;
  }

  //============================== Clear the current materials
  material_xs.clear();
  material_src.clear();

  size_t n_materials = materials.size();
  matid_to_xs_map.assign(n_materials, -1);
  matid_to_src_map.assign(n_materials, -1);

  //============================== Loop over materials
  for (const int& material_id : unique_material_ids)
  {
    auto material = materials[material_id];

    //============================== Loop over properties
    bool found_xs = false;
    for (const auto& property : material->properties)
    {
      //============================== Parse cross sections
      if (property->type() == MaterialPropertyType::CROSS_SECTIONS)
      {
        auto xs = std::static_pointer_cast<CrossSections>(property);
        Assert(xs->n_groups >= groups.size(),
               "Cross sections with fewer groups than the simulation "
               "encountered.\nAll cross sections must have at least as many "
               "groups as the simulation.");

        material_xs.emplace_back(xs);
        matid_to_xs_map[material_id] = material_xs.size() - 1;
        found_xs = true;
      }

        //============================== Parse multigroup sources
      else if (property->type() == MaterialPropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto src = std::static_pointer_cast<IsotropicMGSource>(property);
        Assert(src->values.size() >= groups.size(),
               "Multigroup source encountered with fewer groups than the "
               "simulation.\nAll sources must have at least as many groups "
               "as the simulation.");

        material_src.emplace_back(src);
        matid_to_src_map[material_id] = material_src.size() - 1;
      }
    }//for properties

    Assert(found_xs, "Cross sections not found on a provided material.");
  }//for materials

  // Define the number of groups
  n_groups = groups.size();

  /* Define the precursor properties. For the total number of precursors
     * tally the number of unique decay constants. */
  if (use_precursors)
  {
    max_precursors = 0;
    std::set<double> unique_decay_constants;
    for (const auto& xs : material_xs)
    {
      for (size_t j = 0; j < xs->n_precursors; ++j)
        unique_decay_constants.insert(xs->precursor_lambda[j]);

      if (xs->n_precursors > max_precursors)
        max_precursors = xs->n_precursors;
    }
    n_precursors = unique_decay_constants.size();

    if (n_precursors == 0)
      use_precursors = false;
  }//if use precursors

  std::cout << "Materials initialized: " << materials.size() << "\n";
}