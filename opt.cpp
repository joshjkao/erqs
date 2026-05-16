#include "operations.h"
#include "optimization.h"
#include "quantumstate.h"
#include <chrono>
// #include <iostream>
#include <cmath>
#include <print>
#include <random>
#include <ranges>

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

int main(int argc, char **argv) {
  std::random_device rd;
  auto seed = rd();
  std::mt19937 gen{seed};
  auto root = RandomSumState(2, 2, 2, ~QSpace{0}, gen);

  SquashPures(root);
  RemoveSingles(root);

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

  PauliHamiltonian H{{0, 0, 0, 0, -0.5, -0.5, -0.5, -0.5},
                     {op1, op2, op3, op4, int1, int2, int3, int4}};

  double lambda = 1.0;

  if (argc > 1)
    lambda = atof(argv[1]);

  std::println("seed: {}, lambda: {}", seed, lambda);
  std::println(
      "iteration num_coeffs_proposed E_internal_proposed E_total_proposed "
      "accepted num_coeffs E_internal_min E_total_min time");

  double curr_min_cost = INFINITY;
  double curr_energy = INFINITY;

  for (auto i : std::views::iota(0, 50)) {
    auto root_cpy = Clone(root);
    random_update(root_cpy, gen);
    Normalize(root_cpy);
    Prune(root_cpy, 1e-3);
    RemoveSingles(root_cpy);
    SquashPures(root_cpy);

    auto time = time_in_ms([&]() { OptimizeCoefficients(H, root_cpy); });

    double energy = ExpectedValue(H, root_cpy);
    double cost = energy + lambda * CountCoefficients(root_cpy);
    bool accepted = false;
    if (cost < curr_min_cost) {
      root = root_cpy;
      curr_min_cost = cost;
      curr_energy = energy;
      accepted = true;
    }
    std::println("{} {} {} {} {} {} {} {} {}", i, CountCoefficients(root_cpy),
                 energy, cost, accepted, CountCoefficients(root), curr_energy,
                 curr_min_cost, time);
    Print(root);
    std::println();
  }
}
