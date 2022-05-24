#include "mesh.h"
#include "cell.h"
#include "face.h"

using namespace pdes;

std::string
Grid::coordinate_system_str(const CoordinateSystem coord_sys)
{
  switch (coord_sys)
  {
    case CoordinateSystem::CARTESIAN:   return "CARTESIAN";
    case CoordinateSystem::CYLINDRICAL: return "CYLINDRICAL";
    case CoordinateSystem::SPHERICAL:   return "SPHERICAL";
    default:                            return "UNDEFINED";
  }
}


std::string
Grid::cell_type_str(const CellType cell_type)
{
  switch (cell_type)
  {
    case CellType::SLAB:    return "SLAB";
    case CellType::ANNULUS: return "ANNULUS";
    case CellType::SHELL:   return "SHELL";
    default:                return "UNDEFINED";
  }
}


std::string
Grid::Cell::str() const
{
  std::stringstream ss;
  ss << "***** Cell " << id << " *****\n";
  ss << "type: " << cell_type_str(type) << "\n";
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
        << face.str();

  return ss.str();
}


std::string
Grid::Face::str() const
{
  std::stringstream ss;
  size_t v = 0;
  ss << "n_vertex_ids: " << vertex_ids.size() << "\n";
  for (auto& v_id : vertex_ids)
    ss << "Vertex " << v++ << ": " << v_id << "\n";

  ss << "normal: " << normal << "\n";
  ss << "centroid: " << centroid << "\n";
  ss << "area: " << area << "\n";
  ss << "has_neighbor: " << has_neighbor << "\n";
  ss << "neighbor_id: " << neighbor_id << "\n";

  return ss.str();
}


std::ostream&
Grid::operator<<(std::ostream& os, const Cell& cell)
{ return os << cell.str(); }


std::ostream&
Grid::operator<<(std::ostream& os, const Face& face)
{ return os << face.str(); }
