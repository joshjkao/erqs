#include "operations.hpp"
#include "optimization.hpp"
#include "quantumstate.hpp"
#include "validation.hpp"
#include <chrono>
#include <iostream>
#include <memory>
#include <random>

auto time_in_ms(auto func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  return duration.count();
}

int main() {

  std::random_device rd{};
  std::mt19937 gen{rd()};
  auto root = RandomSumState(5, 5, 5, ~QSpace{0}, gen);

  auto clone1 = Clone(root);
  auto clone2 = Clone(root);

  RemoveSingles(clone1);
  SquashPures(clone1);

  Simplify(clone2);

  size_t coeffs_original = CountCoefficients(root);
  size_t coeffs1 = CountCoefficients(clone1);
  size_t coeffs2 = CountCoefficients(clone2);
  size_t nodes_original = CountNodes(root);
  size_t nodes1 = CountNodes(clone1);
  size_t nodes2 = CountNodes(clone2);

  std::cout << nodes_original << " " << nodes1 << " " << nodes2 << "\n";
  std::cout << coeffs_original << " " << coeffs1 << " " << coeffs2 << "\n";
  // std::cout << Equals_slow(clone1, root) << "\n";
  // std::cout << Equals_slow(clone2, root) << "\n";
  std::cout << Validate(clone1, ValidationArgs{}) << "\n";
  std::cout << Validate(clone2, ValidationArgs{}) << "\n";
  // std::cout << "first\n";
  // Print(clone1);
  // std::cout << "second\n";
  // Print(clone2);

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

  // std::cout << "normalize it again\n";
  // Normalize(root);
  // Print(root);

  PauliOperator op1{
      .x = BitString{"1100"}, .y = BitString{"0000"}, .z = BitString{"1000"}};
  PauliOperator op2{
      .x = BitString{"0000"}, .y = BitString{"1000"}, .z = BitString{"0100"}};
  PauliOperator op3{
      .x = BitString{"0000"}, .y = BitString{"0001"}, .z = BitString{"0010"}};
  PauliOperator op4{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0001"}};

  PauliHamiltonian H{{1, 1, 1, 1}, {op1, op2, op3, op4}};

  // auto conj = Operate(H, root);

  // auto inner = Inner(root, conj);

  // std::cout << "original\n";
  // Print(inner);
  Flatten(root);
  std::cout << ExpectedValue(H, root) << "\n";
  std::cout << ExpectedValue(H, clone1) << "\n";
  std::cout << ExpectedValue(H, clone2) << "\n";

  // OptimizeCoefficients(H, root);
}
