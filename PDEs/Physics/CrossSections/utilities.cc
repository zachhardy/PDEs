#include "cross_sections.h"

#include <iostream>
#include <fstream>

#include "macros.h"


void
Physics::CrossSections::reset()
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


void
Physics::CrossSections::read_cross_section(
  const std::string keyword,
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
  int group;
  double value;
  while (line != keyword + "_END")
  {
    line_stream >> group >> value;
    Assert(g < n_groups,
           "The number of entries in the " + keyword +
           " block exceeds the total number of groups.");
    Assert(group < n_groups && group >= 0,
           "Invalid group number encountered on line " +
           std::to_string(line_number) + ".");

    destination[group] = value;

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    ++g;
  }
}


//######################################################################


void
Physics::CrossSections::read_transfer_matrices(
  const std::string keyword,
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
      Assert(moment >= 0,
             "Invalid scattering moment encountered on line " +
             std::to_string(line_number))
      Assert(group < n_groups && group >= 0,
             "Invalid incident group encountered on line " +
             std::to_string(line_number) + ".");
      Assert(gprime < n_groups && gprime >= 0,
             "Invalid destination group encountered on line " +
             std::to_string(line_number) + ".");

      if (moment < destination.size())
        destination[moment][group][gprime] = value;
    }

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}


//######################################################################


void
Physics::CrossSections::read_precursor_property(
  const std::string keyword,
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
  int precursor_num;
  double value;
  while (line != keyword + "_END")
  {
    line_stream >> precursor_num >> value;
    Assert(j < n_precursors,
           "The number of entries in the " + keyword +
           " block exceeds the total number of precursor species.");
    Assert(precursor_num < n_precursors && precursor_num >= 0,
           "Invalid precursor species number encountered on line " +
           std::to_string(line_number) + ".");

    destination[precursor_num] = value;

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    ++j;
  }
}


//######################################################################


void
Physics::CrossSections::read_delayed_spectra(
  const std::string keyword,
  EmissionSpectra& destination,
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
  int group, precursor_num;
  while (line != keyword + "_END")
  {
    line_stream >> word;
    if (word == "G_PRECURSORJ_VAL")
    {
      line_stream >> group >> precursor_num >> value;
      Assert(group < n_groups && group >= 0,
             "Invalid group encountered on line " +
             std::to_sting(line_number) + ".");
      Assert(precursor_num < n_precursors && precursor_num >= 0,
             "Invalid precursor number encountered on line " +
             std::to_string(line_number) + ".");

      destination[group][precursor_num] = value;
    }

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}
