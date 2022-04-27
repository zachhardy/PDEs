#include "steadystate_solver_fv.h"

/**
 * \brief Transfer a groupset vector to a full multigroup vector.
 * \param groupset The groupset the vector \p x corresponds to.
 * \param x The groupset vector to be transferred.
 * \param destination The destination full multigroup vector.
 */
void neutron_diffusion::SteadyStateSolver_FV::
scoped_transfer(const Groupset& groupset, const math::Vector& x,
                math::Vector& destination)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();
  const auto n_gsg = groupset.groups.size();

  for (const auto& cell : mesh->cells)
  {
    const size_t gs_uk_map = cell->id * n_gsg;
    const size_t full_uk_map = cell->id * n_groups;
    for (size_t g = gs_i; g <= gs_f; ++g)
      destination[full_uk_map + g] = x[gs_uk_map + g];
  }
}

/**
 * \brief Copy the elements corresponding to the specified groupset from one
 *        full multigroup vector to another.
 * \param groupset The groupset being copied.
 * \param x The vector to be copied from.
 * \param destination The vector copied into.
 */
void neutron_diffusion::SteadyStateSolver_FV::
scoped_copy(const Groupset& groupset, const math::Vector& x,
            math::Vector& destination)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  for (const auto& cell : mesh->cells)
  {
    const size_t uk_map = cell->id * n_groups;
    for (size_t g = gs_i; g <= gs_f; ++g)
      destination[uk_map + g] = x[uk_map + g];
  }
}

/**
 * \brief Return the \f$\ell_2\f$-norm between the last two iterates of the
 *        multigroup flux elements that belong to the specified groupset.
 * \param groupset The groupset to compute the change for.
 */
double neutron_diffusion::SteadyStateSolver_FV::
compute_change(const Groupset& groupset)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  double norm = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const size_t uk_map = cell->id * n_groups;
    for (size_t g = gs_i; g <= gs_f; ++g)
    {
      const size_t dof = uk_map + g;
      double delta = std::fabs(phi[dof] - phi_ell[dof]);
      norm += delta * delta;
    }
  }
  return std::sqrt(norm);
}
