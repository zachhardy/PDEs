#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include "material.h"

#include <string>
#include <unordered_map>
#include <cstddef>
#include <functional>


namespace PDEs
{
  namespace Physics
  {

    /**
     * A class holding nuclear data for neutron interactions.
     */
    class CrossSections : public MaterialProperty
    {
    protected:
      using TransferMatrix = std::vector<std::vector<double>>;
      using EmissionSpectra = std::vector<std::vector<double>>;

    public:
      /**
       * A convenient typedef for functional cross-sections. This function takes
       * as input a group number, a vector of arguments, and a reference
       * cross-section value. The function must have knowledge of how to parse
       * the vector of arguments.
       */
      using XSFunction = std::function<double(const unsigned int group_num,
                                              const std::vector<double>& args,
                                              const double reference)>;

    public:
      unsigned int n_groups = 0;
      unsigned int n_moments = 0; ///< The number of scattering moments.
      unsigned int n_precursors = 0;

      double density = 1.0; ///< Atom density in <tt>atoms/b-cm</tt>

      bool is_fissile = false;

      std::vector<double> E_bounds; ///< Energy bin boundaries.

      std::vector<double> sigma_t;  ///< Total cross section
      std::vector<double> sigma_a;  ///< Absorption cross section
      std::vector<double> sigma_s;  ///< Scattering cross section
      std::vector<double> sigma_f;  ///< Fission cross section
      std::vector<double> sigma_r;  ///< Removal cross section

      /// Moment-wise group-to-group transfer matrices
      std::vector<TransferMatrix> transfer_matrices;

      std::vector<double> chi;  ///< Total fission spectrum.
      std::vector<double> chi_prompt; ///< Prompt fission spectrum.
      EmissionSpectra chi_delayed;  ///< Delayed fission spectrum.

      std::vector<double> nu;  ///< Total neutrons per fission.
      std::vector<double> nu_prompt;  ///< Prompt neutrons per fission.
      std::vector<double> nu_delayed;  ///< Delayed neutrons per fission.
      std::vector<double> beta;  ///< Delayed neutron fraction.

      std::vector<double> nu_sigma_f;
      std::vector<double> nu_prompt_sigma_f;
      std::vector<double> nu_delayed_sigma_f;

      std::vector<double> precursor_lambda; ///< Decay constants in (s\f$^{-1}\f$).
      std::vector<double> precursor_yield;  ///< Precursor yield fractions.

      std::vector<double> inv_velocity; ///< Inverse speed (s/cm)
      std::vector<double> diffusion_coeff; ///< Diffusion coefficient (cm)
      std::vector<double> buckling; ///< Material buckling term

      /// A modifier function for the absorption cross-sections.
      XSFunction sigma_a_function;

    public:
      //################################################## Constructors

      /** Default constructor. */
      CrossSections();

      /** Delete the current cross-section data. */
      void reset();

      /** Reinitialize the cross-section data. */
      void reinit();

      /** Make pure scatterer. */
      void make_pure_scatterer();

      /**
       * Read a ".xs" file containing the cross-section information. Once
       * the cross-sections are parsed, multiply the relevant quantities by the
       * specified atom density \p rho.
       */
      void
      read_xs_file(const std::string file_name,
                   const double rho = 1.0,
                   const bool verbose = false);

      /**
       * Read a ".ndi" file containing the cross-section information. Once
       * the cross-sections are parsed, multiply the relevant quantities by the
       * specified atom density \p rho.
       */
      void
      read_ndi_file(const std::string file_name,
                    const double rho = 1.0,
                    const bool verbose = false);

      /** Write the energy bin boundaries to a file. */
      void write_group_structure(
          const std::string directory = ".",
          const std::string file_prefix = "e_bounds") const;

    private:
      //################################################## Operations

      /**
       * Compute \f$ \sigma_s \f$ from the zeroth scattering moment.
       *
       * Compute the group-wise scattering cross sections from the zeroth
       * transfer matrices. This is defined as the sum of all transfers from a
       * fixed group  to any other. Mathematically, this is given by \f$
       * \sigma_{s,g} =  \sum_{g^\prime} \sigma_{0, g \rightarrow g^\prime} \f$,
       * which is obtained via column-wise sums.
       */
      void compute_scattering_from_transfers();

      /**
       * Enforce the relationship \f$ \sigma_t = \sigma_a + \sigma_s \f$.
       *
       * If the absorption cross section was not specified, then compute it
       * via \f$ \sigma_a = \sigma_t - \sigma_s \f$, where \f$ \sigma_s \f$ is
       * obtained from the transfer matrix. Otherwise, modify \f$ \sigma_t \f$
       * using the provided absorption cross-section and scattering
       * cross-section obtained from the transfer matrix. If both \f$ \sigma_t
       * \f$ and \f$ \sigma_a \f$ are provided, and they do not agree with the
       * transfer matrix, the \f$ \sigma_a \f$ values are taken as true.
       */
      void reconcile_cross_sections();

      /**
       * Validate and process fission related properties.
       *
       * This routine does a number of things.
       * 1. If the cross sections are not fissile but delayed neutron precursor
       *    properties were specified, they are cleared.
       * 2. If fissile and delayed neutron precursor properties were specified,
       *    checks for prompt and delayed \f$ \nu \f$ and \f$ \chi \f$ are
       *    performed and fission/emission spectra as well as precursor yields
       *    \f$ \gamma \f$ are normalized to unity. Lastly, the total \f$ \nu
       *    \f$ and \f$ \chi \f$ quantities are computed from their prompt and
       *    delayed counterparts. The total neutrons per fission is computed via
       *    \f$ \nu = \nu_p + \nu_d \f$ and the total spectra via
       *    \f[ \chi = (1 - \beta) \chi_p + \beta \sum_j \gamma_j \chi_{d_j}
       *             = \frac{\nu_p}{\nu} \chi_p +
       *               \frac{\nu_d}{\nu} \sum_j \gamma_j \chi_{d,j}.
       *    \f]
       * 3. If fissile and no delayed neutron precursor properties were
       *    specified, checks for total \f$ \nu \f$ and \f$ \chi \f$ are
       *    performed and the fission spectrum is normalized to unity.
       */
      void reconcile_fission_properties();

      /**
       * Compute the macroscopic cross-sections.
       *
       * Compute the macroscopic cross-sections via \f$ \Sigma_x = \rho
       * \sigma_x \f$. If the \p diffusion_coeff was unspecified, it is
       * computed via its standard definition, given by \f$ D =
       * \frac{1}{3 \Sigma_t} \f$.
       */
      void compute_macroscopic_cross_sections();
    };
  }
}
#endif //CROSS_SECTIONS_H
