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
  unsigned int n_precursors; ///< The number of delayed neutron precursors.

  double density = 1.0; ///< The atom-density \f$ \rho \f$ in atoms/b-cm.
  bool is_fissile = false; ///< A flag to indicate fissile properties.
  bool has_precursors = false; ///< A flag to indicate precursor presence.

  std::vector<double> sigma_t;  ///< The total cross sections (b).
  std::vector<double> sigma_a;  ///< The absorption cross sections (b).
  std::vector<double> sigma_s;  ///< The scattering cross sections (b).
  std::vector<double> sigma_f;  ///< The fission cross sections (b).

   /** The removal cross sections (b) given by \f$ \sigma_r = \sigma_t -
    * \sigma_{s, g \rightarrow g} \f$. */
  std::vector<double> sigma_r;

  /** The group-to-group transfer cross sections for each scattering moment,
   *  given by \f$ \sigma_{\ell, g^\prime \rightarrow g} \f$ (b). */
  std::vector<TransferMatrix> transfer_matrices;

  std::vector<double> chi;        ///< The total fission spectrum.
  std::vector<double> chi_prompt; ///< The prompt fission spectrum.
  EmmissionSpectra chi_delayed;   ///< The delayed fission spectrum.

  std::vector<double> nu;         ///< The total neutrons per fission.
  std::vector<double> nu_prompt;  ///< The prompt neutrons per fission.
  std::vector<double> nu_delayed; ///< The delayed neutrons per fission.

  std::vector<double> nu_sigma_f;
  std::vector<double> nu_prompt_sigma_f;
  std::vector<double> nu_delayed_sigma_f;

  /// The delayed neutron precursor decay constants (s\f$^{-1}\f$).
  std::vector<double> precursor_lambda;
  /// The delayed neutron precursor yield fractions.
  std::vector<double> precursor_yield;

  /// The inverse of the neutron speed (s/cm).
  std::vector<double> inv_velocity;
  /// The diffusion coefficient (cm)
  std::vector<double> diffusion_coeff;

public:
  /// Default constructor.
  CrossSections() : MaterialProperty(MaterialPropertyType::CROSS_SECTIONS) {}

  /// Default constructor with a name.
  explicit CrossSections(const std::string property_name)
    : MaterialProperty(property_name, MaterialPropertyType::CROSS_SECTIONS)
  {}

public:
  /// Clear all of the cross section data.
  void reset();

  /**
   * \brief Set the cross section data by reading a file.
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
  /**
   * \brief Read a cross section block from the cross section file.
   * \param keyword The identifier for the current property block.
   * \param destination The cross section vector to store the results in.
   * \param file The file being parsed.
   * \param line_stream Storage for a line in the file.
   * \param line_number The current line number in the file.
   */
  void read_cross_section(const std::string& keyword,
                          std::vector<double>& destination,
                          std::ifstream& file,
                          std::istringstream& line_stream,
                          unsigned int& line_number);

  /**
   * \brief Read the transfer matrix block of the cross section file.
   * \param keyword The identifier for the current property block.
   * \param destination The vector of transfer matrices to store the result in.
   * \param file The file being parsed.
   * \param line_stream Storage for a line in the file.
   * \param line_number The current line number in the file.
   */
  void read_transfer_matrices(const std::string& keyword,
                              std::vector<TransferMatrix>& destination,
                              std::ifstream& file,
                              std::istringstream& line_stream,
                              unsigned int& line_number);

  /**
   * \brief Read a precursor property from the cross section file.
   * \param keyword The identifier for the current property block.
   * \param destination The precursor property vector to store the result in.
   * \param file The file being parsed.
   * \param line_stream Storage for a line in the file.
   * \param line_number The current line number in the file.
   */
  void read_precursor_property(const std::string& keyword,
                               std::vector<double>& destination,
                               std::ifstream& file,
                               std::istringstream& line_stream,
                               unsigned int& line_number);

  /**
   * \brief Read the delayed neutron spectra from the cross section file.
   * \param keyword The identifier for the current property block.
   * \param destination The vector of emmission spectra to store the results in.
   * \param file The file being parsed.
   * \param line_stream Storage for a line in the file.
   * \param line_number The current line number in the file.
   */
  void read_delayed_spectra(const std::string& keyword,
                            EmmissionSpectra& destination,
                            std::ifstream& file,
                            std::istringstream& line_stream,
                            unsigned int& line_number);
};


#endif //CROSS_SECTIONS_H
