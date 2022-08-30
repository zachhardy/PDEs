#ifndef PDES_LIGHTWEIGHT_CROSS_SECTIONS_H
#define PDES_LIGHTWEIGHT_CROSS_SECTIONS_H

#include "CrossSections/cross_sections.h"


namespace PDEs
{
  namespace Physics
  {

    /**
     * A utility class for storing a limited amount of cross-section data.
     *
     * This class stores a pointer to a CrossSections object and is primarily
     * used for storing potentially functional cross-sections so that default
     * values are not changed. At present, this class only accommodates
     * total cross-sections
     */
    class LightWeightCrossSections
    {
    public:
      std::vector<double> sigma_t;

    private:
      std::shared_ptr<CrossSections> ref_xs;

    public:
      /** Default constructor. */
      LightWeightCrossSections(std::shared_ptr<CrossSections> xs);

      /**
       * Update the cross-section values.
       *
       * This routine calls the stored functions within the reference
       * CrossSections for the update. Said functions should have the structure
       * of the \p args argument encoded with in.
       */
      void update(const std::vector<double>& args);


    };
  }
}

#endif //PDES_LIGHTWEIGHT_CROSS_SECTIONS_H
