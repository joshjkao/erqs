#include "hamiltonian.hpp"
#include "validation.hpp"
#include <ranges>

PauliHamiltonian RandomHamiltonian(size_t n_terms, QSpace space,
                                   std::mt19937 &gen) {
  auto random_bitset = [&]() -> BitString {
    std::bernoulli_distribution dist{0.5};
    BitString bits;
    for (size_t i = 0; i < NQUBITS; ++i) {
      bits[i] = dist(gen);
    }
    return bits & space;
  };

  std::vector<double> coeffs;
  std::vector<PauliOperator> ops;

  for (auto _ : std::views::iota(0uz, n_terms)) {
    BitString x = random_bitset();
    BitString y = random_bitset() & ~x;
    BitString z = random_bitset() & ~x & ~y;
    coeffs.push_back(1);
    ops.emplace_back(x, y, z);
  }

  PauliHamiltonian H{coeffs, ops};
  if (!Validate(H, CHECK_ALL)) {
    throw std::runtime_error("invalid hamiltonian");
  }

  return H;
}
