#include "fv.h"


using namespace Math;


FiniteVolume::
FiniteVolume(std::shared_ptr<Grid::Mesh> reference_mesh) :
  Discretization(reference_mesh, DiscretizationMethod::FINITE_VOLUME)
{}


size_t
FiniteVolume::nodes_per_cell() const
{ return 1; }


size_t
FiniteVolume::dofs_per_cell(const size_t n_components) const
{ return n_components; }


size_t
FiniteVolume::n_nodes() const
{ return mesh->cells.size(); }


size_t
FiniteVolume::n_dofs(const size_t n_components) const
{ return n_components * mesh->cells.size(); }


std::vector<Grid::Point>
FiniteVolume::nodes(const Grid::Cell& cell) const
{ return {cell.centroid}; }
