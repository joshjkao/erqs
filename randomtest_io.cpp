
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
      "n_factors", "factors per product state", cxxopts::value<size_t>());

  options.parse_positional({"max_depth", "n_terms", "n_factors"});

  auto result = options.parse(argc, argv);

  // just error out
  size_t max_depth = result["max_depth"].as<size_t>();
  size_t n_terms = result["n_terms"].as<size_t>();
  size_t n_factors = result["n_factors"].as<size_t>();
  QSpace space = ~QSpace{0};

  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto rand1 = RandomSumState(max_depth, n_terms, n_factors, space, gen);

  std::string str = Stringify(rand1);

  auto copy = FromString(str);

  std::string copy_str = Stringify(copy);

  if (str == copy_str) {
    std::cout << "string comparison passed\n";
  } else {
    std::cout << "string comparison failed\n";
  }

  if (Equals_literal(rand1, copy)) {
    std::cout << "literal comparison passed\n";
  } else {
    std::cout << "literal comparison failed\n";
  }
}
