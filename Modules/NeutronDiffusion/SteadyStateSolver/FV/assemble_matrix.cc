#include "steadystate_solver_fv.h"


using namespace pdes;
using namespace Math;


void
NeutronDiffusion::SteadyStateSolver_FV::
assemble_matrix(Groupset& groupset)
{
  SparseMatrix& A = groupset.matrix;
  A = 0.0;

  // Get groupset range
  const auto n_gsg = groupset.groups.size();
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    const size_t i = cell.id * n_gsg;

    //==================== Total interaction terms
    const double* sig_t = xs->sigma_t.data();
    for (size_t g = gs_i; g <= gs_f; ++g)
      A.add(i + g, i + g, sig_t[g] * volume);

    //============================== Cross-group coupling
    if (solution_technique == SolutionTechnique::FULL_SYSTEM)
    {
      for (size_t g = gs_i; g <= gs_f; ++g)
      {
        //==================== Scattering term
        const double* sig_s = xs->transfer_matrices[0][g].data();
        for (size_t gp = gs_i; gp <= gs_f; ++gp)
          A.add(i + g, i + gp, -sig_s[gp] * volume);

        //==================== Fission term
        if (xs->is_fissile)
        {
          //========== Total fission
          if (not use_precursors)
          {
            const double chi = xs->chi[g];
            const double* nu_sigf = xs->nu_sigma_f.data();
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
              A.add(i + g, i + gp, -chi * nu_sigf[gp] * volume);
          }

          //========== Prompt + delayed fission
          else
          {
            //========== Prompt
            const double chi_p = xs->chi_prompt[g];
            const double* nup_sigf = xs->nu_prompt_sigma_f.data();
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
              A.add(i + g, i + gp, -chi_p * nup_sigf[gp] * volume);

            //========== Delayed
            const double* chi_d = xs->chi_delayed[g].data();
            const double* gamma = xs->precursor_yield.data();
            const double* nud_sigf = xs->nu_delayed_sigma_f.data();
            for (size_t j = 0; j < xs->n_precursors; ++j)
              for (size_t gp = gs_i; gp <= gs_f; ++gp)
              {
                double fission_d = chi_d[j] * gamma[j] * nud_sigf[gp];
                A.add(i + g, i + gp, -fission_d * volume);
              }
          }
        }//if fissile
      }//for group
    }//if full system

    //================================================== Loop over faces
    for (const auto& face : cell.faces)
    {
      //============================== Interior faces
      if (face.has_neighbor)
      {
        // Get neighbor info
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const size_t j = nbr_cell.id * n_gsg;

        const int nbr_xs_id = matid_to_xs_map[nbr_cell.material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];

        // Geometric quantities
        const double d_pf = cell.centroid.distance(face.centroid);
        const double d_pn = cell.centroid.distance(nbr_cell.centroid);
        const double w = d_pf / d_pn; // harmonic mean weighting factor

        //==================== Diffusion term
        const double* D = xs->diffusion_coeff.data();
        const double* D_nbr = nbr_xs->diffusion_coeff.data();
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          const double D_eff = 1.0 / (w/D[g] + (1.0 - w)/D_nbr[g]);
          const double value = D_eff / d_pn * face.area;

          A.add(i + g, i + g, value);
          A.add(i + g, j + g, -value);
        }
      }//if interior face

      //============================== Boundary faces
      else
      {
        const auto bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        //==================== Dirichlet boundary term
        if (bndry_type == BoundaryType::ZERO_FLUX or
            bndry_type == BoundaryType::DIRICHLET)
        {
          const double* D = xs->diffusion_coeff.data();
          const double d_pf = cell.centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
            A.add(i + g, i + g, D[g]/d_pf * face.area);

        }

        //==================== Robin boundary term
        else if (bndry_type == BoundaryType::VACUUM or
                 bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          const double* D = xs->diffusion_coeff.data();
          const double d_pf = cell.centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            double value = bc->a*D[g] / (bc->b*D[g] + bc->a*d_pf);
            A.add(i + g, i + g, value * face.area);
          }
        }
      }//if boundary face
    }//for face
  }//for cell
}


