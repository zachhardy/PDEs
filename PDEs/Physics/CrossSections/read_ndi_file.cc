#include "cross_sections.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>
#include <cassert>


using namespace PDEs;
using namespace Physics;


void
CrossSections::
read_ndi_file(const std::string file_name,
             const double atom_density,
             const bool verbose)
{
  std::cout << "Reading cross-section file \"" << file_name << "\".\n";

  // Clear old cross-sections
  reset();

  density = atom_density;

  // Open the file
  std::ifstream file;
  file.open(file_name);
  assert(file.is_open());

  //############################################################
  /** Lambda for reading vector data. */
  auto read_1d_data =
      [](std::vector<double>& destination,
         const unsigned int n,
         std::ifstream& file,
         std::istringstream& line_stream,
         size_t& line_number)
      {
        std::string line;
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;

        double value;
        unsigned int g = 0;
        while (true)
        {
          while (line_stream >> value)
            destination.at(g++) = value;

          if (g == n) break;

          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //############################################################
  /** Lambda for reading transfer matrix data. */
  auto read_transfer_matrix =
      [](std::vector<TransferMatrix>& destination,
         const unsigned int M,
         const unsigned int G,
         std::ifstream& file,
         std::istringstream& line_stream,
         size_t& line_number)
      {
        std::string line;
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;

        // Skip the 0 moment identifier
        std::getline(file, line);
        line_stream = std::istringstream (line);
        ++line_number;

        double value;
        unsigned int ell = 0, g = 0, gp = 0;
        while (true)
        {
          // Parse all values of the line
          std::vector<double> vals;
          while (line_stream >> value)
            vals.push_back(value);

          // If more than one entry, this is real data. The NDI data is
          // contiguously written and is of a different format than the normal
          // ".xs" file in that the innermost index if the "to" group instead
          // of the "from" group. For this reason, elements are added down
          // columns instead of across rows. When the end of the column is
          // reached, the index counters are moved to the top of the subsequent
          // column.
          if (vals.size() > 1)
            for (const auto& v : vals)
            {
              destination.at(ell).at(g++).at(gp) = v;
              if (g == G) { ++gp; g = 0; }
            }

          // If only one entry, this corresponds to the start of the next
          // scattering moment. When this occurs, ensure innermost index is at
          // its maximum value and that the moment index is one less than the
          // designated moment index from the file. This moves the counters to
          // point to the start of the next scattering moment data.
          else
          {
            assert(gp == G && ell + 1 == (unsigned int)vals.front());
            ++ell; g = gp = 0;
          }

          // When the maximum of the innermost index is reached for the last
          // scattering moment, leaving the while loop. This is the end of
          // the scattering data.
          if (ell == M - 1 && gp == G) break;

          // Otherwise, move to the next line to keep parsing
          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };


  // Read the file
  size_t line_number = 0;
  std::string line, word;
  while (std::getline(file, line))
  {
    std::istringstream line_stream(line);
    line_stream >> word;

    auto& f = file;
    auto& ls = line_stream;
    auto& ln = line_number;

    if (word == "num_grps")
    {
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;

      line_stream >> n_groups;
      reinit();
    }

    if (word == "vel")
    {
      read_1d_data(inv_velocity, n_groups, f, ls, ln);
      for (auto& v: inv_velocity) v = 1.0 / v;
    }
    if (word == "sig_tot")
      read_1d_data(sigma_t, n_groups, f, ls, ln);
    if (word == "pn_order")
    {
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;

      line_stream >> n_moments;
      transfer_matrices.resize(n_moments);
      for (auto& transfer_matrix : transfer_matrices)
        transfer_matrix.resize(n_groups, std::vector<double>(n_groups, 0.0));

      // Move to the "pn_full" line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;

      read_transfer_matrix(transfer_matrices, n_moments, n_groups, f, ls, ln);
    }

    word = "";
    ++line_number;
  }// while open
  file.close();

  compute_scattering_from_transfers();
  reconcile_cross_sections();
  reconcile_fission_properties();
  compute_macroscopic_cross_sections();

  if (verbose)
    std::cout << "Cross Sections Details:\n"
              << "\t# of Groups    : " << n_groups << "\n"
              << "\t# of Precursors: " << n_precursors << "\n"
              << "\tFissile?       : " << is_fissile << "\n";

  std::cout << "Done reading cross-sections file \"" << file_name << "\".\n";
}
