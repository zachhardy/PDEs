#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include "../material.h"

#include <string>
#include <unordered_map>

/**
 * A class for material cross sections for use in neutron diffusion
 * and neutron transport applications.
 */
class CrossSections : public MaterialProperty
{
protected:
  typedef std::vector<std::vector<double>> TransferMatrix;
  typedef std::vector<std::vector<double>> EmmissionSpectra;

public:
  unsigned int n_groups; ///< The number of energy groups.
  unsigned int scattering_order; ///< The maximum scattering order.

  /// The total number of delayed neutron precursors.
  unsigned int n_precursors;

  /// The atom-density \f$ \rho \f$ of the material in atoms/b-cm.
  double density = 1.0;
  /// A flag for fissile materials.
  bool is_fissile = false;
  /// A flag for the existence of delayed neutron precursors.
  bool has_precursors = false;

  std::vector<double> sigma_t;  ///< The total cross sections (b).
  std::vector<double> sigma_a;  ///< The absorption cross sections (b).
  std::vector<double> sigma_s;  ///< The scattering cross sections (b).
  std::vector<double> sigma_f;  ///< The fission cross sections (b).

  /**
   * \brief The removal cross section (b).
   *
   * \f[ \sigma_r = \sigma_t - \sigma_{s, g \rightarrow g} \f]
   */
  std::vector<double> sigma_r;

   /// The group-to-group transfer cross sections for each scattering moment,
   /// given by \f$ \sigma_{\ell, g^\prime \rightarrow g} \f$ (b).
  std::vector<TransferMatrix> transfer_matrices;

  /**
   * \brief The total fission spectrum.
   *
   * \f[ \chi = (1 - \beta) \chi_p + \beta \sum_j \gamma_j \chi_{d,j}
   *          = \frac{\nu_p}{\nu}) \chi_p
   *            + \frac{\nu_d}{\nu} \sum_j \gamma_j \chi_{d,j} \f]
   */
  std::vector<double> chi;
  std::vector<double> chi_prompt;   ///< The prompt fission spectrum.
  EmmissionSpectra chi_delayed;     ///< The delayed fission spectrum.

  /**
   * \brief The total neutrons per fission.
   *
   * \f[ \nu = \nu_p + \nu_d \f]
   */
  std::vector<double> nu;
  std::vector<double> nu_prompt;  ///< The prompt neutrons per fission.
  std::vector<double> nu_delayed; ///< The delayed neutrons per fission.

  std::vector<double> nu_sigma_f;
  std::vector<double> nu_prompt_sigma_f;
  std::vector<double> nu_delayed_sigma_f;

  /// The decay constants of each delayed neutron precursor (s\f$^{-1}\f$).
  std::vector<double> precursor_lambda;
  /**
   * \brief The yield fraction for the delayed neutron precursors.
   *
   * \f[ \sum_j \gamma_j = 1 \f]
   */
  std::vector<double> precursor_yield;

  /// The inverse of the neutron speed (s/cm).
  std::vector<double> inv_velocity;

  /// The diffusion coefficient (cm)
  std::vector<double> diffusion_coeff;

public:
  /// Default constructor.
  CrossSections() : MaterialProperty(MaterialPropertyType::CROSS_SECTIONS) {}

  /**
   * \brief Construct a named property.
   * \param property_name A name for identification purposes.
   */
  explicit CrossSections(const std::string property_name)
    : MaterialProperty(property_name, MaterialPropertyType::CROSS_SECTIONS)
  {}

  /// Clear all properties.
  void reset();

  /**
   * \brief Read a file to set the cross section values.
   * \param file_name The path to the file to read.
   */
  void read_xs_file(const std::string& file_name);

private:

  /// Compute \f$ \sigma_s \f$ from the zeroth scattering moment.
  void compute_scattering_from_transfers();

  /// Enforce the relationship \f$ \sigma_t = \sigma_a + \sigma_s \f$.
  void reconcile_cross_sections();

  /// Validate and process fission related properties.
  void reconcile_fission_properties();

  /// Compute the macroscopic cross sections (cm\f$^{-1}\f$).
  void compute_macroscopic_cross_sections();

private:
  /// Read a cross section block from the cross section file.
  void read_cross_section(const std::string& keyword,
                          std::vector<double>& destination,
                          std::ifstream& file,
                          std::istringstream& line_stream,
                          unsigned int& line_number);

  /// Read the transfer matrix block of the cross section file.
  void read_transfer_matrices(const std::string& keyword,
                              std::vector<TransferMatrix>& destination,
                              std::ifstream& file,
                              std::istringstream& line_stream,
                              unsigned int& line_number);

  /// Read a precursor property from the cross section file.
  void read_precursor_property(const std::string& keyword,
                               std::vector<double>& destination,
                               std::ifstream& file,
                               std::istringstream& line_stream,
                               unsigned int& line_number);

  /// Read the delayed neutron spectra from the cross section file.
  void read_delayed_spectra(const std::string& keyword,
                            EmmissionSpectra& destination,
                            std::ifstream& file,
                            std::istringstream& line_stream,
                            unsigned int& line_number);
};


#endif //CROSS_SECTIONS_H
