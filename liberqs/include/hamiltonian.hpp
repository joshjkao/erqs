#pragma once
#include "common.hpp"
#include <random>
#include <vector>

struct PauliOperator {
  // assume these are always disjoint
  // a 1 means that gate acts on that bit
  BitString x;
  BitString y;
  BitString z;
};

struct PauliHamiltonian {
  std::vector<double> coeffs;
  std::vector<PauliOperator> operators;
};

PauliHamiltonian RandomHamiltonian(size_t n_terms, QSpace space,
                                   std::mt19937 &gen);
