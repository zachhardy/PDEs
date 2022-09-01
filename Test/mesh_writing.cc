#include "mesh.h"
#include "ortho_grids.h"

#include <vector>
#include <iostream>


using namespace PDEs;
using namespace Grid;


int main()
{
  size_t n = 5;
  double h = 1.0 / (n - 1);
  std::vector<double> x_verts(n);
  for (size_t i = 0; i < n + 1; ++i)
    x_verts[i] = i * h;

  auto mesh = create_2d_orthomesh(x_verts, x_verts);

  return 0;
}
