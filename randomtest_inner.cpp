#include "cxxopts.hpp"
#include "operations.h"
// #include "optimization.h"
#include "quantumstate.h"
#include "validation.h"
#include <iostream>
#include <random>

auto main(int argc, char **argv) -> int {
  cxxopts::Options options("randomtest_inner",
                           "attempt to validate inner products");
  options.add_options()("max_depth", "maximum tree depth",
                        cxxopts::value<size_t>())(
      "n_terms", "terms per sum state", cxxopts::value<size_t>())(
      "n_factors", "factors per product state",
      cxxopts::value<size_t>())("p,print_state", "print the state");

  options.parse_positional({"max_depth", "n_terms", "n_factors"});

  auto result = options.parse(argc, argv);

  // just error out
  size_t max_depth = result["max_depth"].as<size_t>();
  size_t n_terms = result["n_terms"].as<size_t>();
  size_t n_factors = result["n_factors"].as<size_t>();
  // bool print_state = result["print_state"].as<bool>();
  QSpace space = ~QSpace{0};

  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto ket = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto bra = RandomSumState(max_depth, n_terms, n_factors, space, gen);

  auto inner_slow = Inner_slow(bra, ket);
  auto inner_fast = Inner(bra, ket);

  std::cout << "{\n";
  std::cout << "\t\"inner_slow\": " << inner_slow << ",\n";
  std::cout << "\t\"inner_fast\": "
            << (inner_fast.ketbras.empty() ? 0.0 : inner_fast.ketbras[0].coeff)
            << "\n";
  std::cout << "}";
}
