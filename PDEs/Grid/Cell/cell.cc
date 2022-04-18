#include "cell.h"

#include <cmath>


//######################################################################
/// Get the Cell type as a string.
std::string cell_type_name(const CellType cell_type)
{
  switch (cell_type)
  {
    case CellType::SLAB:      return "SLAB";
    case CellType::ANNULUS:   return "ANNULUS";
    case CellType::SHELL:     return "SHELL";
    default:                  return "NONE";
  }
}


//######################################################################
Cell::Cell(const Cell& other)
  : type(other.type),
    id(other.id),
    material_id(other.material_id),
    centroid(other.centroid),
    volume(other.volume),
    vertex_ids(other.vertex_ids),
    faces(other.faces)
{}


//######################################################################
Cell::Cell(Cell&& other)
    : type(other.type),
    id(other.id),
    material_id(other.material_id),
    centroid(other.centroid),
    volume(other.volume),
    vertex_ids(std::move(other.vertex_ids)),
    faces(std::move(other.faces))
{}


//######################################################################
Cell& Cell::operator=(const Cell& other)
{
  id = other.id;
  material_id = other.material_id;
  centroid = other.centroid;
  volume = other.volume;
  vertex_ids = other.vertex_ids;
  faces = other.faces;
  return *this;
}


//######################################################################
std::string Cell::to_string() const
{
  std::stringstream ss;
  ss << "***** Cell " << id << " *****\n";
  ss << "type: " << cell_type_name(type) << "\n";
  ss << "material_id: " << material_id << "\n";
  ss << "centroid: " << centroid << "\n";
  ss << "volume: " << volume << "\n";

  int v = 0;
  ss << "n_vertex_ids: " << vertex_ids.size() << "\n";
  for (auto& v_id : vertex_ids)
    ss << "Vertex " << v++ << ": " << v_id << "\n";

  int f = 0;
  ss << "n_faces: " << faces.size() << "\n";
  for (auto& face : faces)
    ss  << "--- Face " << f++ << " ---\n"
        << face.to_string() ;

  return ss.str();
}
