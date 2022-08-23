#include "mesh.h"
#include "cell.h"
#include "face.h"

#include <iostream>
#include <sstream>


using namespace PDEs;


std::string
Grid::coordinate_system_str(const CoordinateSystemType coord_sys)
{
  switch (coord_sys)
  {
    case CoordinateSystemType::CARTESIAN:
      return "CARTESIAN";
    case CoordinateSystemType::CYLINDRICAL:
      return "CYLINDRICAL";
    case CoordinateSystemType::SPHERICAL:
      return "SPHERICAL";
    default:
      return "UNDEFINED";
  }
}


std::string
Grid::cell_type_str(const CellType cell_type)
{
  switch (cell_type)
  {
    case CellType::SLAB:
      return "SLAB";
    case CellType::ANNULUS:
      return "ANNULUS";
    case CellType::SHELL:
      return "SHELL";
    case CellType::QUADRILATERAL:
      return "QUADRILATERAL";
    default:
      return "UNDEFINED";
  }
}


Grid::Cell::Cell(const CellType cell_type) :
    type(cell_type)
{}


std::string
Grid::Cell::str() const
{
  std::stringstream ss;
  ss << "***** Cell " << id << " *****\n";
  ss << "type: " << cell_type_str(type) << "\n";
  ss << "material_id: " << material_id << "\n";
  ss << "centroid: " << centroid << "\n";
  ss << "volume: " << volume << "\n";

  ss << "n_vertex_ids: " << vertex_ids.size() << "\n";
  ss << "vertices: [ ";
  for (auto& v_id: vertex_ids)
    ss << v_id << " ";
  ss << "]\n";

  int f = 0;
  ss << "n_faces: " << faces.size() << "\n";
  for (auto& face: faces)
    ss << "--- Face " << f++ << " ---\n"
       << face.str();

  return ss.str();
}


std::string
Grid::Face::str() const
{
  std::stringstream ss;
  ss << "n_vertex_ids: " << vertex_ids.size() << "\n";
  ss << "vertices: [ ";
  for (auto& v_id: vertex_ids)
    ss << v_id << " ";
  ss << "]\n";

  ss << "normal: " << normal << "\n";
  ss << "centroid: " << centroid << "\n";
  ss << "area: " << area << "\n";
  ss << "has_neighbor: " << has_neighbor << "\n";
  ss << "neighbor_id: " << neighbor_id << "\n";
  return ss.str();
}


std::ostream&
operator<<(std::ostream& os, const Grid::Cell& cell)
{
  return os << cell.str();
}


std::ostream&
operator<<(std::ostream& os, const Grid::Face& face)
{
  return os << face.str();
}
