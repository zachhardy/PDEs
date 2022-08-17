#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <vector>
#include <memory>


namespace Physics
{
  enum class MaterialPropertyType
  {
    SCALAR = 0,              ///< A scalar-valued property.
    CROSS_SECTIONS = 1,      ///< Neutronics cross sections.
    ISOTROPIC_MG_SOURCE = 2  ///< Isotropic neutron source.
  };


  /** Abstract base class for material properties. */
  class MaterialProperty
  {
  private:
    const MaterialPropertyType property_type;
    const std::string property_name = "Generic Property";

  public:
    MaterialProperty(const MaterialPropertyType type);
    MaterialProperty(const MaterialPropertyType type, const std::string name);

    MaterialPropertyType type() const;
    std::string name() const;
  };


  //######################################################################

  /**
   * A class representing a material. In order to support future multi-physics
   * applications. A Physics holds a collection a MaterialProperty objects which
   * hold derived objects that describe properties that are necessary for a
   * particular type of physics.
   */
  class Material
  {
  private:
    const std::string material_name = "Generic Material";

  public:
    std::vector<std::shared_ptr<MaterialProperty>> properties;

  public:
    Material() = default;
    explicit Material(const std::string name);

    std::string name() const;
  };

  //######################################################################

  /** A simple scalar valued material property. */
  class ScalarProperty : public MaterialProperty
  {
  public:
    double value = 1.0;

  public:
    ScalarProperty();
    ScalarProperty(const std::string name);
    ScalarProperty(const double value);
    ScalarProperty(const double value, const std::string name);
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
    IsotropicMultiGroupSource(const std::vector<double> src);

    IsotropicMultiGroupSource(const std::vector<double> src,
                              const std::string name);

  };

}
#endif //MATERIAL_H
