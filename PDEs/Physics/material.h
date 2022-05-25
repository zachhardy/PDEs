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
private:
  const MaterialPropertyType property_type;
  const std::string property_name = "Generic Property";

public:
  /**
   * Construct a material property of the specified type.
   */
  MaterialProperty(const MaterialPropertyType type);

  /**
   * Construct a named material property of the specified type.
   */
  MaterialProperty(const MaterialPropertyType type, const std::string name);

  /**
   * Get the material property type.
   */
  MaterialPropertyType
  type() const;

  /**
   * Get the material property name.
   */
  std::string
  name() const;
};


//######################################################################

/**
 * A class representing a material.
 *
 * In order to support future multi-physics applications. A Physics holds
 * a collection a MaterialProperty objects which hold derived objects that
 * describe properties that are necessary for a particular type of physics.
 */
class Material
{
private:
  const std::string material_name = "Generic Material";

public:
  std::vector<std::shared_ptr<MaterialProperty>> properties;

public:
  /**
   * Default constructor.
   */
   Material() = default;

  /**
   * Construct a material with unique name.
   */
  explicit
  Material(const std::string name);

  /**
   * Return the material name.
   */
  std::string
  name() const;
};

//######################################################################

/**
 * A simple scalar valued material property.
 */
class ScalarProperty : public MaterialProperty
{
public:
  double value = 1.0;

public:
  /**
   * Construct a scalar valued property with the default value of 1.
   */
  ScalarProperty();

  /**
   * Construct a named scalar valued property with the default value of 1.
   */
  explicit
  ScalarProperty(const std::string name);

  /**
   * Construct a scalar valued property with the specified value.
   */
  ScalarProperty(const double value);

  /**
   * Construct a scalar valued property with the specified name and value.
   */
  ScalarProperty(const double value,
                 const std::string name);
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
  IsotropicMultiGroupSource(const std::vector<double> src);

  /**
   * Construct a named multigroup source from an STL vector.
   */
  IsotropicMultiGroupSource(const std::vector<double> src,
                            const std::string name);

};

}
#endif //MATERIAL_H
