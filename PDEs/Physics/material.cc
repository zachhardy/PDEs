#include "material.h"


using namespace Physics;


MaterialProperty::
MaterialProperty(const MaterialPropertyType type) :
  property_type(type)
{}


MaterialProperty::
MaterialProperty(const MaterialPropertyType type, const std::string name) :
  property_type(type), property_name(name)
{}


MaterialPropertyType
MaterialProperty::type() const
{ return property_type; }


std::string
MaterialProperty::name() const
{ return property_name; }


//######################################################################


Material::Material(const std::string name) :
  material_name(name)
{}


std::string
Material::name() const
{ return material_name; }


//######################################################################


ScalarProperty::ScalarProperty() :
  MaterialProperty(MaterialPropertyType::SCALAR)
{}


ScalarProperty::ScalarProperty(const std::string name) :
  MaterialProperty(MaterialPropertyType::SCALAR, name)
{}


ScalarProperty::ScalarProperty(const double value) :
  MaterialProperty(MaterialPropertyType::SCALAR), value(value)
{}


ScalarProperty::
ScalarProperty(const double value, const std::string name) :
  MaterialProperty(MaterialPropertyType::SCALAR, name), value(value)
{}

//######################################################################

IsotropicMultiGroupSource::
IsotropicMultiGroupSource(const std::vector<double> src) :
  MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
  values(src)
{}


IsotropicMultiGroupSource::
IsotropicMultiGroupSource(const std::vector<double> src,
                          const std::string name) :
  MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE, name),
  values(src)
{}
