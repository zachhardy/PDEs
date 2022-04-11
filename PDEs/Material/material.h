#ifndef MATERIAL_H
#define MATERIAL_H

#include "material_property.h"

#include <string>
#include <vector>
#include <memory>


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

  /**
   * \brief Construct a named material.
   * \param material_name A material_name for identification purposes.
   */
  explicit Material(std::string material_name) : name(material_name) {}
};

#endif //MATERIAL_H
