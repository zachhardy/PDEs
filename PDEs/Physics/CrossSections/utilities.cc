#include "cross_sections.h"

#include <iostream>
#include <fstream>

#include "macros.h"


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
    assert(g < n_groups);
    assert(group < n_groups && group >= 0);

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
      assert(moment >= 0);
      assert(group < n_groups && group >= 0);
      assert(gprime < n_groups && gprime >= 0);

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
    assert(j < n_precursors);
    assert(precursor_num < n_precursors && precursor_num >= 0);

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
      assert(group < n_groups && group >= 0);
      assert(precursor_num < n_precursors && precursor_num >= 0);

      destination[group][precursor_num] = value;
    }

    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}
