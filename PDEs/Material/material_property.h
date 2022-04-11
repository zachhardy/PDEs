#ifndef MATERIAL_PROPERTY_H
#define MATERIAL_PROPERTY_H

#include <string>
#include <vector>

/// Material property types.
enum class PropertyType
{
  SCALAR = 0,                       ///< A scalar-valued property.
  CROSS_SECTIONS = 1,               ///< Neutronics cross sections.
  ISOTROPIC_MG_SOURCE= 2   ///< Isotropic neutron source.
};


/// An abstract class for representing a material property.
class MaterialProperty
{
public:
  const PropertyType type;
  const std::string name = "Generic Property";

public:
  /// Default constructor.
  explicit MaterialProperty(const PropertyType property_type)
    : type(property_type)
  {}

  /// Construct a property with a name.
  explicit MaterialProperty(const std::string property_name,
                            const PropertyType property_type)
    : name(property_name), type(property_type)
  {}
};


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
  double value = 1.0;

public:
  /// Default constructor.
  ScalarProperty() : MaterialProperty(PropertyType::SCALAR) {}

  /// Construct a scalar property with the specified value.
  explicit ScalarProperty(const double scalar_value)
    : MaterialProperty(PropertyType::SCALAR), value(scalar_value)
  {}

  /// Construct a named scalar property set to zero.
  explicit ScalarProperty(const std::string property_name)
    : MaterialProperty(property_name, PropertyType::SCALAR)
  {}

  /// Construct a named scalar property with the specified value.
  explicit ScalarProperty(const std::string& property_name,
                          const double scalar_value)
    : MaterialProperty(property_name, PropertyType::SCALAR),
      value(scalar_value)
  {}
};


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
  std::vector<double> values;

public:
  /// Default constructor.
  IsotropicMultiGroupSource()
    : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE)
  {}

  /// Construct a named, empty source.
  IsotropicMultiGroupSource(const std::string& property_name)
    : MaterialProperty(property_name, PropertyType::ISOTROPIC_MG_SOURCE)
  {}

  /// Construct the source based on the provided STL vector.
  IsotropicMultiGroupSource(const std::vector<double>& mg_values)
    : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE), values(mg_values)
  {}

  /// Construct a named source from the provided STL vector.
  IsotropicMultiGroupSource(const std::string& property_name,
                            const std::vector<double>& mg_values)
    : MaterialProperty(property_name, PropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  /// Move constructor with the provided STL vector.
  IsotropicMultiGroupSource(std::vector<double>&& mg_values)
    : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE), values(mg_values)
  {}

  /// Move constructor with a name and the provided STL vector.
  IsotropicMultiGroupSource(const std::string& property_name,
                            std::vector<double>&& mg_values)
    : MaterialProperty(property_name, PropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

  /// Construct a source fram the provided initializer list.
  IsotropicMultiGroupSource(std::initializer_list<double>& mg_values)
    : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE), values(mg_values)
  {}

  /// Construct a named source from the provided iinitializer list.
  IsotropicMultiGroupSource(const std::string& property_name,
                            std::initializer_list<double>& mg_values)
    : MaterialProperty(property_name, PropertyType::ISOTROPIC_MG_SOURCE),
      values(mg_values)
  {}

};

#endif //MATERIAL_PROPERTY_H
