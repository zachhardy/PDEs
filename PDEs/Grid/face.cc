#include "mesh.h"
#include "cell.h"
#include "face.h"

#include <iostream>
#include <sstream>


using namespace PDEs;
using namespace Grid;


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
Grid::operator<<(std::ostream& os, const Grid::Face& face)
{
  return os << face.str();
}
