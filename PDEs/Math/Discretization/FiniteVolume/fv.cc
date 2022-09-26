#include "fv.h"


using namespace PDEs;
using namespace Math;


FiniteVolume::
FiniteVolume(std::shared_ptr<Grid::Mesh> reference_mesh) :
    Discretization(reference_mesh, SpatialDiscretizationMethod::FINITE_VOLUME)
{
  std::cout << "Creating Finite Volume spatial discretization.\n";
  std::cout << "Finished creating Finite Volume spatial discretization.\n";
}


size_t
FiniteVolume::n_nodes() const
{
  return mesh->cells.size();
}


unsigned int
FiniteVolume::nodes_per_cell() const
{
  return 1;
}


size_t
FiniteVolume::n_dofs(const unsigned int n_components) const
{
  return n_components * mesh->cells.size();
}


unsigned int
FiniteVolume::dofs_per_cell(const unsigned int n_components) const
{
  return n_components;
}


std::vector<Grid::Node>
FiniteVolume::nodes(const Grid::Cell& cell) const
{
  return {cell.centroid};
}


void
FiniteVolume::
make_sparsity_pattern(std::vector<std::vector<size_t>> pattern,
                      const unsigned int n_components,
                      const bool is_coupled) const
{
  // Resive based on the number of DoFs
  pattern.resize(n_dofs(n_components));

  // Loop over cells
  for (const auto& cell: mesh->cells)
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
    for (const auto& face: cell.faces)
    {
      if (face.has_neighbor)
      {
        const size_t jr = face.neighbor_id * n_components;
        for (unsigned int c = 0; c < n_components; ++c)
          pattern[ir + c].emplace_back(jr + c);
      }
    }//for face
  }//for cells
}
