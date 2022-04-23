#include "cross_sections.h"

#include <iostream>
#include <sstream>
#include <fstream>


/// Clear all of the cross section data.
void physics::CrossSections::reset()
{
  n_groups = 0;
  n_precursors = 0;
  scattering_order = 0;

  is_fissile = false;
  density = 1.0;

  sigma_t.clear();
  sigma_a.clear();
  sigma_s.clear();
  sigma_r.clear();
  sigma_f.clear();

  chi.clear();
  chi_prompt.clear();
  chi_delayed.clear();

  nu.clear();
  nu_prompt.clear();
  nu_delayed.clear();

  nu_sigma_f.clear();
  nu_prompt_sigma_f.clear();
  nu_delayed_sigma_f.clear();

  precursor_lambda.clear();
  precursor_yield.clear();

  inv_velocity.clear();
  diffusion_coeff.clear();

  transfer_matrices.clear();
}

//######################################################################

/**
 * \brief Read a cross section block from the cross section file.
 * \param keyword The identifier for the current property block.
 * \param destination The cross section vector to store the results in.
 * \param file The file being parsed.
 * \param line_stream Storage for a line in the file.
 * \param line_number The current line number in the file.
 */
void physics::CrossSections::read_cross_section(
    const std::string& keyword,
    std::vector<double>& destination,
    std::ifstream& file,
    std::istringstream& line_stream,
    size_t& line_number)
{
  std::string line;

  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;

  //========== Go through entries
  size_t g = 0;
  int group; double value;
  while (line != keyword + "_END")
  {
    line_stream >> group >> value;

    if (g >= n_groups)
    {
      std::stringstream err;
      err << "CrossSections::" << __FUNCTION__ << ": "
          << "Too man entries in " << keyword << " block. "
          << "There can only be " << n_groups << " entries.";
      throw std::runtime_error(err.str());
    }
    if (group >= n_groups or group < 0)
    {
      std::stringstream err;
      err << "CrossSections::" << __FUNCTION__ << ": "
          << "Invalid group number encountered on line "
          << line_number << ".";
      throw std::runtime_error(err.str());
    }

    destination[group] = value;

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number; ++g;
  }
}

//######################################################################

/**
 * \brief Read the transfer matrix block of the cross section file.
 * \param keyword The identifier for the current property block.
 * \param destination The vector of transfer matrices to store the result in.
 * \param file The file being parsed.
 * \param line_stream Storage for a line in the file.
 * \param line_number The current line number in the file.
 */
void physics::CrossSections::read_transfer_matrices(
    const std::string& keyword,
    std::vector<TransferMatrix>& destination,
    std::ifstream& file,
    std::istringstream& line_stream,
    size_t& line_number)
{
  std::string word, line;

  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;

  //========== Go through entries
  double value;
  int moment, group, gprime;
  while (line != keyword + "_END" and !file.eof())
  {
    line_stream >> word;
    if (word == "M_GPRIME_G_VAL")
    {
      line_stream >> moment >> gprime >> group >> value;

      if (group >= n_groups or group < 0) {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Invalid incident group encountered on line "
            << line_number << ".";
        throw std::runtime_error(err.str());
      }
      if (gprime >= n_groups or group < 0) {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Invalid destination group encountered on line "
            << line_number << ".";
        throw std::runtime_error(err.str());
      }

      if (moment < destination.size())
        destination[moment][group][gprime] = value;
    }

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}

//######################################################################

/**
 * \brief Read a precursor property from the cross section file.
 * \param keyword The identifier for the current property block.
 * \param destination The precursor property vector to store the result in.
 * \param file The file being parsed.
 * \param line_stream Storage for a line in the file.
 * \param line_number The current line number in the file.
 */
void physics::CrossSections::read_precursor_property(
    const std::string& keyword,
    std::vector<double>& destination,
    std::ifstream& file,
    std::istringstream& line_stream,
    size_t& line_number)
{
  std::string line;

  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;

  //========== Go through entries
  size_t j = 0;
  int precursor_num; double value;
  while (line != keyword + "_END")
  {
    line_stream >> precursor_num >> value;

    if (j >= n_precursors)
    {
      std::stringstream err;
      err << "CrossSections::" << __FUNCTION__ << ": "
          << "Too man entries in " << keyword << " block. "
          << "There can only be " << n_precursors << " entries.";
      throw std::runtime_error(err.str());
    }
    if (precursor_num >= n_precursors or precursor_num < 0)
    {
      std::stringstream err;
      err << "CrossSections::" << __FUNCTION__ << ": "
          << "Too man entries in " << keyword << " block. "
          << "Invalid destination precursor number.";
      throw std::runtime_error(err.str());
    }

    destination[precursor_num] = value;

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number; ++j;
  }
}


//######################################################################

/**
 * \brief Read the delayed neutron spectra from the cross section file.
 * \param keyword The identifier for the current property block.
 * \param destination The vector of emmission spectra to store the results in.
 * \param file The file being parsed.
 * \param line_stream Storage for a line in the file.
 * \param line_number The current line number in the file.
 */
void physics::CrossSections::read_delayed_spectra(
    const std::string& keyword,
    EmmissionSpectra& destination,
    std::ifstream& file,
    std::istringstream& line_stream,
    size_t& line_number)
{
  std::string word, line;

  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;

  line_stream >> word;

  //========== Go through entries
  double value;
  int group, precursor_num;
  while (line != keyword + "_END")
  {
    if (word == "G_PRECURSORJ_VAL")
    {
      line_stream >> group >> precursor_num >> value;
      if (group >= n_groups or group < 0 or
          precursor_num >= n_precursors or precursor_num < 0) {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "The group or precursor number exceeded the maximum "
            << "on line " << line_number << ".";
        throw std::runtime_error(err.str());
      }
      destination[group][precursor_num] = value;
    }

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}
