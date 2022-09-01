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
    protected:
      const MaterialPropertyType property_type;

    public:
      /** Construct a material property of the specified \p type. */
      explicit MaterialProperty(const MaterialPropertyType type);

      /** Return the material property type. */
      MaterialPropertyType type() const;
    };


    /**
     * A class representing a material. In order to support future multi-physics
     * applications, materials hold a collection a material properties which
     * are then parsed by the relevant physics solver for access.
     */
    class Material
    {
    public:
      std::vector<std::shared_ptr<MaterialProperty>> properties;

    private:
      const std::string material_name = "Generic Material";

    public:
      /** Default constructor. Construct an empty material. */
      Material() = default;

      /** Construct a named empty material. */
      explicit Material(const std::string name);

      /** Return the name of the material. */
      std::string name() const;
    };


    /**
     * A scalar valued material property.
     *
     * This can be used to represent scalar coefficients in generic equations.
     * By default, the scalar value is set to unity.
     */
    class ScalarProperty : public MaterialProperty
    {
    public:
      double value;

    public:
      /** Default constructor. Construct a scalar property set to \p value. */
      ScalarProperty(const double value = 1.0);
    };


    /**
     * An isotropic multi-group source material property. This is meant to be
     * used as a vector-valued property in multi-group radiation diffusion and
     * transport simulations.
     */
    class IsotropicMultiGroupSource : public MaterialProperty
    {
    public:
      std::vector<double> values;

    public:
      /** Construct an isotropic multigroup source from an STL vector. */
      IsotropicMultiGroupSource(const std::vector<double>& src);
    };
  }// Physics
}// PDEs
#endif //MATERIAL_H
