#include "fv.h"


using namespace PDEs;
using namespace Math;


FiniteVolume::
FiniteVolume(std::shared_ptr<Grid::Mesh> reference_mesh) :
  Discretization(reference_mesh, SpatialDiscretizationMethod::FINITE_VOLUME)
{}


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
  return n_components*mesh->cells.size();
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
