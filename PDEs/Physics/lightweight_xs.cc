#include "lightweight_xs.h"


using namespace Physics;

LightWeightCrossSections::
LightWeightCrossSections(std::shared_ptr<CrossSections> xs)
  : ref_xs(xs), sigma_t(xs->sigma_t)
{}


void LightWeightCrossSections::
update(const double time, const double temperature,
       const double reference_temperature)
{
  if (ref_xs->sigma_a_function)
  {
    const auto& f = ref_xs->sigma_a_function;
    for (unsigned int g = 0; g < ref_xs->n_groups; ++g)
    {
      const double sig_a = f(g, time, temperature,
                             reference_temperature, ref_xs->sigma_a[g]);
      sigma_t[g] = sig_a + ref_xs->sigma_s[g];
    }
  }
}

