#include "fv.h"

/**
 * Define the sparsity pattern.
 *
 * This routine defines the number of non-zero entries per row for a problem
 * with the specified number of components. If the \p is_coupled flag is set
 * to \p true, it is assumed that all components are coupled to one another,
 * otherwise, it is assumed that the system is uncoupled in all components.
 * For FV discretizations, each node is coupled to those on the neighboring
 * cells.
 *
 * \param prealloc The number of entries per row to preallocate.
 * \param n_components The number of components in the solution.
 * \param is_coupled A flag for allocating storage for coupling between
 *   solution components.
 */
void math::FiniteVolume::
make_sparsity_pattern(std::vector<size_t> prealloc,
                      const size_t n_components, const bool is_coupled) const
{
  // Resive based on the number of DoFs
  prealloc.resize(this->n_dofs(n_components));

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const size_t ir = cell->id * n_components;

    for (size_t c = 0; c < n_components; ++c)
      prealloc[ir + c] += is_coupled ? n_components : 1;

    // Loop over faces
    for (const auto& face : cell->faces)
    {
      if (face.has_neighbor)
      {
        const size_t jr = face.neighbor_id * n_components;
        for (size_t c = 0; c < n_components; ++c)
          prealloc[ir + c] += 1;
      }
    }//for face
  }//for cells
}
