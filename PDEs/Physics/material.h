#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <vector>
#include <memory>


namespace PDEs
{
  namespace Physics
  {
    /**
     * The available types of material properties.
     */
    enum class MaterialPropertyType
    {
      SCALAR = 0,  ///< A scalar-valued property.
      CROSS_SECTIONS = 1, ///< Neutronics cross sections.
      ISOTROPIC_MG_SOURCE = 2  ///< Isotropic neutron source.
    };


    /**
     * Abstract base class for material properties.
     */
    class MaterialProperty
    {
    public:
      /**
       * Construct a material property of the specified \p type.
       */
      explicit MaterialProperty(const MaterialPropertyType type);

      /**
       * Return the material property type.
       */
      MaterialPropertyType
      type() const;

    protected:
      /**
       * The material property type. This is used to distinguish between
       * property types when stored as a pointer to this base class.
       */
      const MaterialPropertyType property_type;
    };


    /**
     * A class representing a material. In order to support future multi-physics
     * applications, materials hold a collection a material properties which
     * are then parsed by the relevant physics solver for access.
     */
    class Material
    {
    public:
      /**
       * A list of material properties that belong to this material.
       */
      std::vector<std::shared_ptr<MaterialProperty>> properties;

      /**
       * Default constructor. Construct a generic empty material.
       */
      Material() = default;

      /**
       * Construct a named empty material.
       */
      explicit Material(const std::string name);

      /**
       * Return the name of the material.
       */
      std::string
      name() const;

    private:
      /**
       * The name of the material.
       */
      const std::string material_name = "Generic Material";
    };


    /**
     * A scalar valued material property. This can be used to represent scalar
     * coefficients in generic equations.
     */
    class ScalarProperty : public MaterialProperty
    {
    public:
      /**
       * The scalar value of the property.
       */
      double value;

      /**
       * Default constructor. Construct a scalar property set to \p value. The
       * default behavior is to set the property value to 1.0.
       */
      ScalarProperty(const double value = 1.0);
    };


    /**
     * An isotropic multi-group source material property. This is meant to be used
     * as a vector-valued property in multi-group radiation diffusion and
     * transport simulations.
     */
    class IsotropicMultiGroupSource : public MaterialProperty
    {
    public:
      /**
       * The multi-group source values.
       */
      std::vector<double> values;

      /**
       * Default constructor. Construct an isotropic multi-group source from a
       * vector of values.
       */
      IsotropicMultiGroupSource(const std::vector<double>& src);

    };

  }
}
#endif //MATERIAL_H
