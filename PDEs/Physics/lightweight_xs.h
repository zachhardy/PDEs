#ifndef PDES_LIGHTWEIGHT_CROSS_SECTIONS_H
#define PDES_LIGHTWEIGHT_CROSS_SECTIONS_H

#include "CrossSections/cross_sections.h"

namespace Physics
{
  class LightWeightCrossSections
  {
  private:
    std::shared_ptr<CrossSections> ref_xs;

  public:
    std::vector<double> sigma_t;
    std::vector<double> diffusion_coeff;
    std::vector<double> buckling;

    LightWeightCrossSections(std::shared_ptr<CrossSections> xs);

    void update(const double time);
  };
}


#endif //PDES_LIGHTWEIGHT_CROSS_SECTIONS_H
