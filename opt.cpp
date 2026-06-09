#include "operations.h"
#include "optimization.h"
#include "quantumstate.h"
#include "validation.h"
#include <chrono>
// #include <iostream>
#include <cmath>
#include <iostream>
#include <memory>
#include <print>
#include <random>
#include <ranges>
#include <variant>

long time_in_ms(auto func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  return duration.count();
}

void random_update(std::shared_ptr<SumState> root, std::mt19937 &gen) {
  std::uniform_int_distribution<int> dist(0, 1);
  int rand = dist(gen);
  if (rand) {
    std::println("attempting a random term addition");
    AddRandomTerm(root, gen);
  } else {
    std::println("attempting a random unification");
    UnifyRandom(root, gen);
  }
}

int main() {
  std::random_device rd;
  auto seed = rd();
  std::mt19937 gen{seed};
  // auto root = RandomSumState(3, 2, 2, QSpace{"1111"}, gen);
  auto root = ZeroOneTensor(QSpace{"1111"});

  Simplify(root);
  Validate(root, CHECK_ALL);

  Print(root);
  std::cout << std::endl;

  auto clone = Clone(root);

  // Z local terms
  PauliOperator op1{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"1000"}};
  PauliOperator op2{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0100"}};
  PauliOperator op3{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0010"}};
  PauliOperator op4{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0001"}};

  // XX interaction terms
  PauliOperator int1{
      .x = BitString{"1100"}, .y = BitString{"0000"}, .z = BitString{"0000"}};
  PauliOperator int2{
      .x = BitString{"0110"}, .y = BitString{"0000"}, .z = BitString{"0000"}};
  PauliOperator int3{
      .x = BitString{"0011"}, .y = BitString{"0000"}, .z = BitString{"0000"}};
  PauliOperator int4{
      .x = BitString{"1001"}, .y = BitString{"0000"}, .z = BitString{"0000"}};

  PauliHamiltonian H{{2, 2, 2, 2, 1, 1, 1, 1},
                     {op1, op2, op3, op4, int1, int2, int3, int4}};

  Visitor visitor{
      .sum_visitor =
          [&](auto state) { OptimizeCoefficients_local(H, root, state); },
      .prod_visitor = [](auto) {},
      .pure_visitor = [](auto) {}};

  InvokeBottomUp(visitor, root);
  double minf_local = ExpectedValue(H, root);
  std::cout << "localized min " << minf_local << "\n";

  double minf_global = OptimizeCoefficients(H, clone);
  std::cout << "global min " << minf_global << " ";

  Normalize(clone);
  Print(clone);

  Normalize(root);
  Print(root);

  auto inner_self1 = Inner_slow(root, root);
  std::cout << inner_self1 << "\n";
  auto inner_self2 = Inner_slow(root, root);
  std::cout << inner_self2 << "\n";

  auto inner = Inner_slow(clone, root);
  std::cout << "inner between the two optimized states: " << inner << "\n";
  std::cout << inner * std::conj(inner);
}
