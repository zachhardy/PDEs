#include "mesh.h"
#include "cell.h"
#include "face.h"

#include <iostream>
#include <sstream>


using namespace PDEs;
using namespace Grid;


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

//######################################################################

Cell::Cell(const CellType cell_type) :
    type(cell_type)
{}


std::string
Cell::str() const
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



std::ostream&
Grid::operator<<(std::ostream& os, const Grid::Cell& cell)
{
  return os << cell.str();
}
