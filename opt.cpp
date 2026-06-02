#include "operations.hpp"
#include "optimization.hpp"
#include "quantumstate.hpp"
#include "validation.hpp"
#include <chrono>
#include <iostream>
#include <memory>
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

int main() {
  std::random_device rd;
  auto seed = rd();
  std::mt19937 gen{seed};
  // auto root = RandomSumState(3, 2, 2, QSpace{"1111"}, gen);
  QSpace space{"1111"};
  auto root = ZeroOneTensor(space);

  Simplify(root);
  Validate(root, CHECK_ALL);

  Print(root);
  std::cout << std::endl;

  auto clone = Clone(root);

  std::vector<double> coeffs;
  std::vector<PauliOperator> z_locals;

  QSpace mask{1};
  QSpace zero{0};
  for (auto _ : std::views::iota(0uz, NQUBITS)) {
    coeffs.push_back(1.0);
    PauliOperator op{.x = zero, .y = mask, .z = zero};
    z_locals.push_back(op);
    mask = mask << 1;
  }

  PauliHamiltonian H{.coeffs = coeffs, .operators = z_locals};

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
