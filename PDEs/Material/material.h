#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <vector>
#include <memory>

namespace pdes::Physics
{

/**
 * Available types of material properties.
 */
enum class MaterialPropertyType
{
  SCALAR = 0,              ///< A scalar-valued property.
  CROSS_SECTIONS = 1,      ///< Neutronics cross sections.
  ISOTROPIC_MG_SOURCE = 2  ///< Isotropic neutron source.
};


/**
 * Abstract base class for material properties.
 */
class MaterialProperty
{
public:
  const MaterialPropertyType type;
  const std::string name = "Generic Property";

public:
  /**
   * Construct a material property of the specified type.
   */
  MaterialProperty(const MaterialPropertyType property_type) :
      type(property_type)
  {}

  /**
   * Construct a material property of the specified type
   * with the specified name.
   */
  MaterialProperty(const MaterialPropertyType property_type,
                   const std::string property_name) :
      type(property_type), name(property_name)
  {}
};


//######################################################################

/**
 * A class representing a material.
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

  /**
   * Construct a material with unique name.
   */
  explicit
  Material(const std::string material_name) : name(material_name) {}
};

//######################################################################

/**
 * A simple scalar valued material property.
 */
class ScalarProperty : public MaterialProperty
{
private:
  double value = 1.0;

public:
  /**
   * Construct a scalar valued property with the default value of 1.
   */
  ScalarProperty() :
      MaterialProperty(MaterialPropertyType::SCALAR)
  {}

  /**
   * Construct a scalar valued property with the default value of 1 with
   * the specified name.
   */
  explicit
  ScalarProperty(const std::string property_name) :
      MaterialProperty(MaterialPropertyType::SCALAR, property_name)
  {}

  /**
   * Construct a scalar valued property with the specified value
   */
  ScalarProperty(const double scalar_value) :
      MaterialProperty(MaterialPropertyType::SCALAR),
      value(scalar_value)
  {}


  /**
   * Construct a scalar valued property with the specified name and value.
   */
  ScalarProperty(const double scalar_value,
                 const std::string property_name) :
      MaterialProperty(MaterialPropertyType::SCALAR),
      value(scalar_value)
  {}
};

//######################################################################

/**
 * An isotropic multigroup source material property. This is meant to be used
 * as a vector-valued property in multigroup radiation diffusion and transport
 * simulations.
 */
class IsotropicMultiGroupSource : public MaterialProperty
{
public:
  std::vector<double> values;

public:
  /**
   * Construct a multigroup source from an STL vector.
   */
  IsotropicMultiGroupSource(const std::vector<double> src) :
      MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(src)
  {}

  /**
   * Construct a multigroup source with the specified name from an STL vector.
   */
  IsotropicMultiGroupSource(const std::vector<double> src,
                            const std::string property_name) :
      MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE,
                       property_name),
      values(src)
  {}

};

}
#endif //MATERIAL_H
