#include "material.h"


namespace PDEs
{
  namespace Physics
  {
    MaterialProperty::
    MaterialProperty(const MaterialPropertyType type) :
        property_type(type)
    {}


    MaterialPropertyType MaterialProperty::type() const
    {
      return property_type;
    }


    Material::Material(const std::string name) :
        material_name(name)
    {}


    std::string Material::name() const
    {
      return material_name;
    }


    ScalarProperty::ScalarProperty(const double value) :
        MaterialProperty(MaterialPropertyType::SCALAR), value(value)
    {}


    IsotropicMultiGroupSource::
    IsotropicMultiGroupSource(const std::vector<double>& src) :
        MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
        values(src)
    {}
  }
}
