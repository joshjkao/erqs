#pragma once
#include "common.hpp"
#include <random>
#include <vector>

struct PauliOperator {
  // assume these are always disjoint
  // this is checked upon construction
  // a 1 means that gate acts on that bit
  BitString x;
  BitString y;
  BitString z;
};

struct PauliHamiltonian {
  std::vector<double> coeffs;
  std::vector<PauliOperator> operators;
};

PauliHamiltonian ReadFromFile(const std::string& filename);

PauliHamiltonian RandomHamiltonian(size_t n_terms, QSpace space,
                                   std::mt19937 &gen);

void Print(const PauliHamiltonian& h);
