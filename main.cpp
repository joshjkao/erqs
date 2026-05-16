#include "operations.h"
#include "optimization.h"
#include "quantumstate.h"
#include <chrono>
#include <iostream>
#include <random>

uint64_t time_in_ms(auto func) {
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
  auto root = RandomSumState(4, 2, 2, ~QSpace{0}, gen);

  std::cout << "before adjustment\n";
  Print(root);

  std::cout << "now add random term\n";
  AddRandomTerm(root, gen);
  Print(root);

  std::cout << "now unify random product\n";
  UnifyRandom(root, gen);
  Print(root);

  std::cout << "now normalize it\n";
  Normalize(root);
  Print(root);

  std::cout << "now prune\n";
  Prune(root, 0.1);
  Print(root);

  std::cout << "now remove singles\n";
  RemoveSingles(root);
  Print(root);

  std::cout << "squash pures\n";
  SquashPures(root);
  Print(root);

  std::cout << "normalize it again\n";
  Normalize(root);
  Print(root);

  PauliOperator op1{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"1000"}};
  PauliOperator op2{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0100"}};
  PauliOperator op3{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0010"}};
  PauliOperator op4{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0001"}};

  PauliHamiltonian H{{1, 1, 1, 1}, {op1, op2, op3, op4}};

  auto conj = Operate(H, root);

  auto inner = Inner(root, conj);

  std::cout << "original\n";
  Print(inner);

  // OptimizeCoefficients(H, root);
}
