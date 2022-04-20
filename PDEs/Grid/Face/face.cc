#include "face.h"

Face::Face(const Face& other)
  : vertex_ids(other.vertex_ids),
    normal(other.normal),
    centroid(other.centroid),
    area(other.area),
    has_neighbor(other.has_neighbor),
    neighbor_id(other.neighbor_id)
{}


//######################################################################

Face::Face(Face&& other)
  : vertex_ids(std::move(other.vertex_ids)),
    normal(other.normal),
    centroid(other.centroid),
    area(other.area),
    has_neighbor(other.has_neighbor),
    neighbor_id(other.neighbor_id)
{}

//######################################################################

Face& Face::operator=(const Face& other)
{
  vertex_ids = other.vertex_ids;
  normal = other.normal;
  centroid = other.centroid;
  area = other.area;
  has_neighbor = other.has_neighbor;
  neighbor_id = other.neighbor_id;
  return *this;
}

//######################################################################

std::string Face::to_string() const
{
  std::stringstream ss;
  int v = 0;
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
