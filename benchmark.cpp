#include "operations.h"
#include "optimization.h"
#include "quantumstate.h"
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
  size_t max_depth, n_terms, n_factors, n_bits;

  if (argc < 5) {
    std::cout << "not enough arguments\n";
    return -1;
  }

  max_depth = static_cast<size_t>(atoi(argv[1]));
  n_terms = static_cast<size_t>(atoi(argv[2]));
  n_factors = static_cast<size_t>(atoi(argv[3]));
  n_bits = static_cast<size_t>(atoi(argv[4]));

  size_t bit_mask = 0;
  for (size_t i = 0; i < n_bits; ++i) {
    bit_mask = bit_mask << 1 | 1;
  }

  std::random_device rd{};
  std::mt19937 gen{rd()};
  auto root =
      RandomSumState(max_depth, n_terms, n_factors, QSpace{bit_mask}, gen);
  // pay attention to the order here
  RemoveSingles(root);
  SquashPures(root);
  // root = Simplify(root);
  // root = Simplify(root);
  Normalize(root);

  // Print(root);

  SkewOperator inner;

  auto do_inner = [&]() { inner = Inner(root, root); };

  auto time = time_in_ms(do_inner);

  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cout << "{\n";
  std::cout << "\"max_depth\": " << max_depth << ",\n";
  std::cout << "\"n_terms\": " << n_terms << ",\n";
  std::cout << "\"n_factors\": " << n_factors << ",\n";
  std::cout << "\"n_bits\": " << n_bits << ",\n";
  std::cout << "\"n_nodes\": " << CountNodes(root) << ",\n";
  std::cout << "\"c_time\": " << time << ",\n";
  std::cout << "\"mem_max\": " << usage.ru_maxrss << ",\n";
  std::cout << "\"resulting_terms\": " << inner.ketbras.size() << ",\n";
  std::cout << "\"inner\": " << inner.ketbras[0].coeff.real() << "\n";
  std::cout << "}\n";

  // std::cout << "before adjustment\n";
  // Print(root);
  // std::cout << "now add random term\n";
  // AddRandomTerm(root, gen);
  // Print(root);
  // std::cout << "now unify random product\n";
  // UnifyRandom(root, gen);
  // Print(root);
  // std::cout << "now normalize it\n";
  // Normalize(root);
  // Print(root);
  // std::cout << "now prune\n";
  // Prune(root, 0.1);
  // Print(root);
  // std::cout << "now remove singles\n";
  // RemoveSingles(root);
  // Print(root);
  // std::cout << "squash pures\n";
  // SquashPures(root);
  // Print(root);

  // PauliOperator op1{
  //     .x = BitString{"0000"}, .y = BitString{"0000"}, .z =
  //     BitString{"1000"}};
  // PauliOperator op2{
  //     .x = BitString{"0000"}, .y = BitString{"0000"}, .z =
  //     BitString{"0100"}};
  // PauliOperator op3{
  //     .x = BitString{"0000"}, .y = BitString{"0000"}, .z =
  //     BitString{"0010"}};
  // PauliOperator op4{
  //     .x = BitString{"0000"}, .y = BitString{"0000"}, .z =
  //     BitString{"0001"}};

  // PauliHamiltonian H{{1, 1, 1, 1}, {op1, op2, op3, op4}};

  // OptimizeCoefficients(H, root);
}
