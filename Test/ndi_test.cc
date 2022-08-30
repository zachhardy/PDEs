#include "material.h"
#include "CrossSections/cross_sections.h"

#include <iostream>


using namespace PDEs;
using namespace Physics;


int main()
{
  auto xs = std::make_shared<CrossSections>();
  xs->read_ndi_file("xs_data/1001.ndi");

  return 0;
}