#ifndef PDES_LIGHTWEIGHT_CROSS_SECTIONS_H
#define PDES_LIGHTWEIGHT_CROSS_SECTIONS_H

#include "CrossSections/cross_sections.h"


namespace PDEs
{
  namespace Physics
  {

    /**
     * A utility class for storing a limited amount of cross-section data. This
     * is primarily used for storing current functional cross-section values.
     */
    class LightWeightCrossSections
    {
    public:
      /**
       * A copy of potentially modified total cross-section values.
       */
      std::vector<double> sigma_t;

      /**
       * Default constructor.
       */
      LightWeightCrossSections(std::shared_ptr<CrossSections> xs);

      /**
       * Update the cross-section values stored in this object using functions
       * within the cross-sections this points to.
       */
      void
      update(const std::vector<double>& args);

    private:
      /**
       * A pointer to cross-section to obtain nominal values from.
       */
      std::shared_ptr<CrossSections> ref_xs;
    };
  }
}

#endif //PDES_LIGHTWEIGHT_CROSS_SECTIONS_H
