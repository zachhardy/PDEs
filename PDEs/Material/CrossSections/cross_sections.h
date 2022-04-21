#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include "../material.h"

#include <string>
#include <unordered_map>

namespace material
{

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
  unsigned int n_groups;
  unsigned int scattering_order;
  unsigned int n_precursors;

  double density = 1.0; ///< Atom density in atoms/b-cm.
  bool is_fissile = false;

  std::vector<double> sigma_t;  ///< Total cross section
  std::vector<double> sigma_a;  ///< Absorption cross section
  std::vector<double> sigma_s;  ///< Scattering cross section
  std::vector<double> sigma_f;  ///< Fission cross section
  std::vector<double> sigma_r;  ///< Removal cross section

  /// Moment-wise group-to-group transfer matrices
  std::vector<TransferMatrix> transfer_matrices;

  std::vector<double> chi;        ///< Total fission spectrum.
  std::vector<double> chi_prompt; ///< Prompt fission spectrum.
  EmmissionSpectra chi_delayed;   ///< Delayed fission spectrum.

  std::vector<double> nu;         ///< Total neutrons per fission.
  std::vector<double> nu_prompt;  ///< Prompt neutrons per fission.
  std::vector<double> nu_delayed; ///< Delayed neutrons per fission.

  std::vector<double> nu_sigma_f;
  std::vector<double> nu_prompt_sigma_f;
  std::vector<double> nu_delayed_sigma_f;

  std::vector<double> precursor_lambda; ///< Decay constants in (s\f$^{-1}\f$).
  std::vector<double> precursor_yield;  ///< Precursor yield fractions.

  std::vector<double> inv_velocity; ///< Inverse speed (s/cm)
  std::vector<double> diffusion_coeff; ///< Diffusion coefficient

public:
  CrossSections() : MaterialProperty(MaterialPropertyType::CROSS_SECTIONS)
  {}

  explicit CrossSections(const std::string property_name)
    : MaterialProperty(property_name, MaterialPropertyType::CROSS_SECTIONS)
  {}

public:
  void reset();
  void read_xs_file(const std::string& file_name);

private:
  void compute_scattering_from_transfers();
  void reconcile_cross_sections();
  void reconcile_fission_properties();
  void compute_macroscopic_cross_sections();

private:
  void read_cross_section(const std::string& keyword,
                          std::vector<double>& destination,
                          std::ifstream& file,
                          std::istringstream& line_stream,
                          size_t & line_number);

  void read_transfer_matrices(const std::string& keyword,
                              std::vector<TransferMatrix>& destination,
                              std::ifstream& file,
                              std::istringstream& line_stream,
                              size_t & line_number);

  void read_precursor_property(const std::string& keyword,
                               std::vector<double>& destination,
                               std::ifstream& file,
                               std::istringstream& line_stream,
                               size_t& line_number);

  void read_delayed_spectra(const std::string& keyword,
                            EmmissionSpectra& destination,
                            std::ifstream& file,
                            std::istringstream& line_stream,
                            size_t& line_number);
};

}

#endif //CROSS_SECTIONS_H
