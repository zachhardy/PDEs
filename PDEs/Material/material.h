#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <vector>
#include <memory>

namespace material
{

// Forward declaration
class MaterialProperty;

enum class MaterialPropertyType
{
  SCALAR = 0,              ///< A scalar-valued property.
  CROSS_SECTIONS = 1,      ///< Neutronics cross sections.
  ISOTROPIC_MG_SOURCE = 2  ///< Isotropic neutron source.
};

inline std::string
material_property_name(const MaterialPropertyType property_type)
{
  switch (property_type)
  {
    case MaterialPropertyType::SCALAR: return "SCALAR";
    case MaterialPropertyType::CROSS_SECTIONS: return "CROSS_SECTIONS";
    case MaterialPropertyType::ISOTROPIC_MG_SOURCE:
    { return "ISOTROPIC_MG_SOURCE"; }
    default: return "UNDEFINED";
  }
}

//######################################################################

/**
 * \brief A class representing a material.
 *
 * In order to support future multi-physics applications. A Material holds
 * a collection a MaterialProperty objects which hold derived objects that
 * describe properties that are necessary for a particular type of physics.
 */
class Material
{
public:
  const std::string name = "Generic Material";
  std::vector<std::shared_ptr<MaterialProperty>> properties;

public:
  explicit Material(const std::string material_name)
      : name(material_name)
  {}
};

//######################################################################

/// An abstract class for representing a material property.
class MaterialProperty
{
public:
  const MaterialPropertyType type;
  const std::string name = "Generic Property";

public:
  explicit MaterialProperty(const MaterialPropertyType property_type)
      : type(property_type)
  {}

  explicit MaterialProperty(const std::string property_name,
                            const MaterialPropertyType property_type)
      : name(property_name), type(property_type)
  {}
};

//######################################################################

/// A simple scalar material property.
class ScalarProperty : public MaterialProperty
{
public:
  double value = 1.0;

public:
  ScalarProperty()
    : MaterialProperty(MaterialPropertyType::SCALAR)
  {}

  explicit ScalarProperty(const std::string property_name)
    : MaterialProperty(property_name, MaterialPropertyType::SCALAR)
  {}

  explicit ScalarProperty(const double scalar_value)
    : MaterialProperty(MaterialPropertyType::SCALAR), value(scalar_value)
  {}

  explicit ScalarProperty(const std::string property_name,
                          const double scalar_value)
    : MaterialProperty(property_name, MaterialPropertyType::SCALAR),
      value(scalar_value)
  {}
};

//######################################################################

/// An isotropic multigroup source for radiation transport problems.
class IsotropicMultiGroupSource : public MaterialProperty
{
public:
  std::vector<double> values;

public:
  IsotropicMultiGroupSource()
    : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE)
  {}

  IsotropicMultiGroupSource(const std::string property_name)
    : MaterialProperty(property_name, MaterialPropertyType::ISOTROPIC_MG_SOURCE)
  {}

  IsotropicMultiGroupSource(const std::vector<double>& mg_values)
    : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  IsotropicMultiGroupSource(const std::string property_name,
                            const std::vector<double>& mg_values)
    : MaterialProperty(property_name, MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  IsotropicMultiGroupSource(std::vector<double>&& mg_values)
    : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  IsotropicMultiGroupSource(const std::string property_name,
                            std::vector<double>&& mg_values)
    : MaterialProperty(property_name,
                       MaterialPropertyType::ISOTROPIC_MG_SOURCE),
        values(mg_values)
  {}

  IsotropicMultiGroupSource(std::initializer_list<double>& mg_values)
    : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  IsotropicMultiGroupSource(const std::string property_name,
                            std::initializer_list<double>& mg_values)
    : MaterialProperty(property_name,
                       MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}
};

}
#endif //MATERIAL_H
