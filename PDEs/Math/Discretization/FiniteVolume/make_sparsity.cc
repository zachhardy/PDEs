#include "fv.h"


void
Math::FiniteVolume::
make_sparsity_pattern(std::vector<std::vector<size_t>> pattern,
                      const unsigned int n_components,
                      const bool is_coupled) const
{
  // Resive based on the number of DoFs
  pattern.resize(n_dofs(n_components));

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const size_t ir = cell.id * n_components;

    for (unsigned int c = 0; c < n_components; ++c)
    {
      if (is_coupled)
        for (unsigned int cp = 0; cp < n_components; ++cp)
          pattern[ir + c].emplace_back(ir + cp);
      else
        pattern[ir + c].emplace_back(ir + c);

    }

    // Loop over faces
    for (const auto& face : cell.faces)
    {
      if (face.has_neighbor)
      {
        const size_t jr = face.neighbor_id*n_components;
        for (unsigned int c = 0; c < n_components; ++c)
          pattern[ir + c].emplace_back(jr + c);
      }
    }//for face
  }//for cells
}
