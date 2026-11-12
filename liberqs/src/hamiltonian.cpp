#include "hamiltonian.hpp"
#include "validation.hpp"
#include <fstream>
#include <iostream>
#include <ranges>

PauliHamiltonian ReadFromFile(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  PauliHamiltonian hamiltonian;
  std::string line;

  while (std::getline(file, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    std::stringstream ss(line);
    double coeff;
    if (!(ss >> coeff))
      continue;

    PauliOperator op;
    size_t qubit;
    char pauli;

    // Read pairs of (qubit_index, pauli_gate)
    while (ss >> qubit >> pauli) {
      // Note: If BitString uses 0-based indexing and your input
      // file is 1-based, replace `qubit` with `qubit - 1`.
      switch (pauli) {
      case 'X':
      case 'x':
        op.x.set(qubit - 1);
        break;
      case 'Y':
      case 'y':
        op.y.set(qubit - 1);
        break;
      case 'Z':
      case 'z':
        op.z.set(qubit - 1);
        break;
      default:
        break;
      }
    }

    hamiltonian.coeffs.push_back(coeff);
    hamiltonian.operators.push_back(op);
  }

  return hamiltonian;
}

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

#include <iostream>
#include <vector>

void Print(const PauliHamiltonian& h) {
  // check to ensure vectors are aligned
  if (h.coeffs.size() != h.operators.size()) {
    std::cerr << "Error: Number of coefficients and operators do not match.\n";
    return;
  }

  for (size_t i = 0; i < h.coeffs.size(); ++i) {
    std::cout << h.coeffs[i];
    const auto& op = h.operators[i];
    
    size_t num_qubits = op.x.size(); 
    
    for (size_t j = 0; j < num_qubits; ++j) {
      if (op.x[j]) {
        std::cout << " " << (j + 1) << " X";
      } else if (op.y[j]) {
        std::cout << " " << (j + 1) << " Y";
      } else if (op.z[j]) {
        std::cout << " " << (j + 1) << " Z";
      }
    }
    std::cout << "\n";
  }
}
