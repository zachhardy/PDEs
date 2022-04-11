#include "material.h"


//######################################################################
std::string material_property_type(const MaterialPropertyType property_type)
{
  switch (property_type)
  {
    case MaterialPropertyType::SCALAR:
      return "SCALAR";
    case MaterialPropertyType::CROSS_SECTIONS:
      return "CROSS_SECTION";
    case MaterialPropertyType::ISOTROPIC_MG_SOURCE:
      return "ISOTROPIC_MG_SOURCE";
    default:
      return "UNKNOWN";
  }
}
