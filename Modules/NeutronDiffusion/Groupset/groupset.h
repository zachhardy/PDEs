#ifndef GROUPSET_H
#define GROUPSET_H

#include "vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

#include <memory>
#include <vector>
#include <cinttypes>


namespace neutron_diffusion
{

class Groupset
{
public:
  int id;
  std::vector<uint64_t> groups;

  /*---------- Options ----------*/
  uint64_t max_iterations = 100;
  double tolerance = 1.0e-8;

  /*---------- System Storage ----------*/
  math::Matrix<double> matrix;
  math::Vector<double> rhs;

  Groupset() : id(-1) {}
  explicit Groupset(int groupset_num) : id(groupset_num) {}
};

}

#endif //GROUPSET_H
