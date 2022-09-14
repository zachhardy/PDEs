#include <vector>
#include <set>
#include <map>
#include <cmath>

#include <iostream>
#include <iomanip>

#include "Grid/cartesian_vector.h"

using namespace PDEs::Grid;


int main(int argc, char** argv)
{
  int n = 2;

  // Define velocities
  std::vector<CartesianVector> V;
  for (int i = -n; i <= n; ++i)
    for (int j = -n; j <= n; ++j)
      for (int k = -n; k <= n; ++k)
        if (i != 0 && j != 0 && k != 0)
          V.emplace_back(i, j, k);
  std::cout << "Max velocity: " << V.back().length() << "\n"
            << "# of discrete velocities: " << V.size() << "\n";

  for (unsigned int i = 0; i < V.size(); ++i)
    std::cout << i << "  " << V[i].str() << "\n";

  // Define pairs
  std::vector<std::pair<unsigned int, unsigned int>> pairs;
  for (unsigned int i = 0; i < V.size(); ++i)
    for (unsigned int j = 0; j < V.size(); ++j)
        pairs.emplace_back(i, j);

  std::cout << "\n # of Pairs: " << pairs.size() << "\n";

  // Define admissible collisions
  using CollisionPair = std::pair<unsigned int, unsigned int>;
  using AdmissiblePairs = std::vector<CollisionPair>;
  std::vector<std::pair<AdmissiblePairs, double>> a_ijkl(pairs.size());

  std::cout << "\nTransition Probabilities:\n";


  for (unsigned int ij = 0; ij < pairs.size(); ++ij)
  {
    const auto& pair_ij = pairs[ij];
    const auto& v_i = V[pair_ij.first];
    const auto& v_j = V[pair_ij.second];


    unsigned int count = 0;
    for (unsigned int kl = 0; kl < pairs.size(); ++kl)
    {
      const auto& pair_kl = pairs[kl];
      const auto& v_k = V[pair_kl.first];
      const auto& v_l = V[pair_kl.second];

      const auto dp = v_i + v_j - v_k - v_l;
      const auto dE = v_i.dot(v_i) + v_j.dot(v_j) -
                      v_k.dot(v_k) - v_l.dot(v_l);

      if (std::fabs(dp.length()) < 1.0e-12 && std::fabs(dE) < 1.0e-12)
      {
        count += 1;
        a_ijkl[ij].first.push_back(pair_kl);
      }
    }
    a_ijkl[ij].second = 1.0/count;

    if (count > 2)
      std::cout
        << "Pair(" << std::setw(3) << pair_ij.first << ", "
                   << std::setw(3) << pair_ij.second << "): "
        << "# of Pairs " << std::setw(3) << a_ijkl[ij].first.size()
        << ",  Transition Probability: " << a_ijkl[ij].second
        << std::endl;
  }

  return 0;
}
