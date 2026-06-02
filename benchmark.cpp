#include "operations.hpp"
#include "optimization.hpp"
#include "quantumstate.hpp"
#include "validation.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <sys/resource.h>

long time_in_ms(auto func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  return duration.count();
}

int main(int argc, char **argv) {
  size_t max_depth, n_terms, n_factors, n_bits, h_terms;

  if (argc != 6) {
    std::cout << "not enough arguments\n";
    std::cout << "expected 6 arguments, got " << argc << "\n";
    return -1;
  }

  max_depth = static_cast<size_t>(atoi(argv[1]));
  n_terms = static_cast<size_t>(atoi(argv[2]));
  n_factors = static_cast<size_t>(atoi(argv[3]));
  n_bits = static_cast<size_t>(atoi(argv[4]));
  h_terms = static_cast<size_t>(atoi(argv[5]));

  size_t bit_mask = 0;
  for (size_t i = 0; i < n_bits; ++i) {
    bit_mask = bit_mask << 1 | 1;
  }

  std::random_device rd{};
  std::mt19937 gen{rd()};
  auto ket =
      RandomSumState(max_depth, n_terms, n_factors, QSpace{bit_mask}, gen);

  Simplify(ket);

  bool ket_valid = Validate(ket, ValidationArgs{.log_to_stdout = true});

  if (!ket_valid) {
    std::cout << "state is invalid after simplification";
    return -1;
  }

  SkewOperator inner;
  double exp_val;

  PauliHamiltonian H = RandomHamiltonian(5, bit_mask, gen);

  auto do_inner = [&]() { inner = Inner(ket, ket); };
  auto do_exp_val = [&]() { exp_val = ExpectedValue(H, ket); };

  auto inner_time = time_in_ms(do_inner);
  auto exp_val_time = time_in_ms(do_exp_val);

  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cout << "{\n";
  std::cout << "\"max_depth\": " << max_depth << ",\n";
  std::cout << "\"n_terms\": " << n_terms << ",\n";
  std::cout << "\"n_factors\": " << n_factors << ",\n";
  std::cout << "\"n_bits\": " << n_bits << ",\n";
  std::cout << "\"n_h_terms\": " << h_terms << ",\n";
  std::cout << "\"n_nodes\": " << CountNodes(ket) << ",\n";
  std::cout << "\"c_time\": " << inner_time << ",\n";
  std::cout << "\"exp_val_time\": " << exp_val_time << ",\n";
  std::cout << "\"mem_max\": " << usage.ru_maxrss << ",\n";
  std::cout << "\"resulting_terms\": " << inner.size() << ",\n";
  std::cout << "\"inner\": " << inner[0].coeff.real() << "\n";
  std::cout << "}\n";
}
