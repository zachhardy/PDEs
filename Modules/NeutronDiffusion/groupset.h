#ifndef GROUPSET_H
#define GROUPSET_H

#include "vector.h"
#include "matrix.h"
#include "Sparse/sparse_matrix.h"

#include <memory>
#include <vector>
#include <cinttypes>


using namespace Math;


namespace NeutronDiffusion
{

  class Groupset
  {
  public:
    int id = -1;
    std::vector<size_t> groups;

    /*---------- Options ----------*/
    size_t max_iterations = 100;
    double tolerance = 1.0e-8;

    /*---------- System Storage ----------*/
    SparseMatrix A;
    Vector b;


    explicit Groupset(int groupset_num) : id(groupset_num) {}
  };

}

#endif //GROUPSET_H
