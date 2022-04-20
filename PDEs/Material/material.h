#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <vector>
#include <memory>


// Forward declaration
class MaterialProperty;


/// Material property types.
enum class MaterialPropertyType
{
  SCALAR = 0,              ///< A scalar-valued property.
  CROSS_SECTIONS = 1,      ///< Neutronics cross sections.
  ISOTROPIC_MG_SOURCE= 2   ///< Isotropic neutron source.
};

/// Get the material property name.
std::string material_property_name(const MaterialPropertyType property_type);


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
  /// An optional name for identification.
  const std::string name = "Generic Material";

  /// A list of MaterialProperty objects.
  std::vector<std::shared_ptr<MaterialProperty>> properties;

public:
  /// Default constructor.
  Material() = default;
  /// Default constructor with a name.
  explicit Material(const std::string material_name) : name(material_name) {}
};


//######################################################################
/// An abstract class for representing a material property.
class MaterialProperty
{
public:
  /// An identifier for the material property type.
  const MaterialPropertyType type;
  /// An optional name for identification.
  const std::string name = "Generic Property";

public:
  /// Default constructor.
  explicit MaterialProperty(const MaterialPropertyType property_type)
      : type(property_type)
  {}

  /// Default constructor with a name.
  explicit MaterialProperty(const std::string property_name,
                            const MaterialPropertyType property_type)
      : name(property_name), type(property_type)
  {}
};


//######################################################################
/**
 * \brief A class for a scalar-valued MaterialProperty.
 *
 * A scalar-valued property is best used in simple simulations. For example,
 * in a heat conduction problem, a material may be defined by 2 ScalarProperty
 * objects: 1 for the conductivity coefficient and 1 for the heat source.
 * In this example, the solver would need to be able to identify which property
 * corresponds to the conductivity and to the source.
 */
class ScalarProperty : public MaterialProperty
{
public:
  double value = 1.0; ///< The material property value.

public:
  /// Default constructor.
  ScalarProperty() : MaterialProperty(MaterialPropertyType::SCALAR) {}

  /// Default constructor with a name.
  explicit ScalarProperty(const std::string property_name)
      : MaterialProperty(property_name, MaterialPropertyType::SCALAR)
  {}

  /// Construct with the specified value.
  explicit ScalarProperty(const double scalar_value)
      : MaterialProperty(MaterialPropertyType::SCALAR), value(scalar_value)
  {}

  /// Construct with a specified value and a name.
  explicit ScalarProperty(const std::string property_name,
                          const double scalar_value)
      : MaterialProperty(property_name, MaterialPropertyType::SCALAR),
        value(scalar_value)
  {}
};


//######################################################################
/**
 * \brief A class for an isotropic multigroup source.
 *
 * This is meant for use in multigroup neutronics solvers. In general, the
 * number of values specified for the MaterialProperty should match the number
 * of groups in the problem.
 */
class IsotropicMultiGroupSource : public MaterialProperty
{
public:
  std::vector<double> values; ///< The multigroup source values.

public:
  /// Default constructor.
  IsotropicMultiGroupSource()
    : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE)
  {}

  /// Default constructor with a name.
  IsotropicMultiGroupSource(const std::string property_name)
    : MaterialProperty(property_name, MaterialPropertyType::ISOTROPIC_MG_SOURCE)
  {}

  /// Copy constructor with an STL vector.
  IsotropicMultiGroupSource(const std::vector<double>& mg_values)
    : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  /// Copy constructor with a name and an STL vector
  IsotropicMultiGroupSource(const std::string property_name,
                            const std::vector<double>& mg_values)
      : MaterialProperty(property_name, MaterialPropertyType::ISOTROPIC_MG_SOURCE),
        values(mg_values)
  {}

  /// Move constructor with an STL vector
  IsotropicMultiGroupSource(std::vector<double>&& mg_values)
    : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  /// Move constructor with a name and an STL vector.
  IsotropicMultiGroupSource(const std::string property_name,
                            std::vector<double>&& mg_values)
      : MaterialProperty(property_name, MaterialPropertyType::ISOTROPIC_MG_SOURCE),
        values(mg_values)
  {}

  /// Construct with an initializer list.
  IsotropicMultiGroupSource(std::initializer_list<double>& mg_values)
      : MaterialProperty(MaterialPropertyType::ISOTROPIC_MG_SOURCE), values(mg_values)
  {}

  /// Construct a named source and an initializer list.
  IsotropicMultiGroupSource(const std::string property_name,
                            std::initializer_list<double>& mg_values)
      : MaterialProperty(property_name, MaterialPropertyType::ISOTROPIC_MG_SOURCE),
        values(mg_values)
  {}

};


#endif //MATERIAL_H
