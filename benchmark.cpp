#include "operations.h"
#include "optimization.h"
#include "quantumstate.h"
#include "validation.h"
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
  auto clone1 = Clone(root);
  auto clone2 = Clone(root);

  // Print(root);

  bool valid = Validate(root, ValidationArgs{.log_to_stdout = false});
  if (!valid)
    std::cout << "root state is invalid\n";

  RemoveSingles(clone1);
  SquashPures(clone1);
  bool valid_1 = Validate(clone1, ValidationArgs{});
  if (!valid_1)
    std::cout << "removed squashed state is invalid\n";
  if (!Equals_slow(root, clone1)) {
    std::cout << "error: clone1 doesn't equal root\n";
  }

  Simplify(clone2);
  bool valid_2 = Validate(clone2, ValidationArgs{});
  if (!valid_2)
    std::cout << "removed squashed state is invalid\n";
  if (!Equals_slow(root, clone2)) {
    std::cout << "error: clone2 doesn't equal root\n";
  }

  // SkewOperator inner;

  // auto do_inner = [&]() { inner = Inner(root, root); };

  // auto time = time_in_ms(do_inner);

  // struct rusage usage;
  // getrusage(RUSAGE_SELF, &usage);

  // std::cout << "{\n";
  // std::cout << "\"max_depth\": " << max_depth << ",\n";
  // std::cout << "\"n_terms\": " << n_terms << ",\n";
  // std::cout << "\"n_factors\": " << n_factors << ",\n";
  // std::cout << "\"n_bits\": " << n_bits << ",\n";
  // std::cout << "\"n_nodes\": " << CountNodes(root) << ",\n";
  // std::cout << "\"c_time\": " << time << ",\n";
  // std::cout << "\"mem_max\": " << usage.ru_maxrss << ",\n";
  // std::cout << "\"resulting_terms\": " << inner.ketbras.size() << ",\n";
  // std::cout << "\"inner\": " << inner.ketbras[0].coeff.real() << "\n";
  // std::cout << "}\n";

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
