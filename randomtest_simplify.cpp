#include "cxxopts.hpp"
#include "optimization.h"
#include "quantumstate.h"
#include "validation.h"
#include <iostream>
#include <random>

auto main(int argc, char **argv) -> int {
  cxxopts::Options options(
      "randomtest_simplify",
      "generate a random tree, simplify it, and check if it's valid");
  options.add_options()("max_depth", "maximum tree depth",
                        cxxopts::value<size_t>())(
      "n_terms", "terms per sum state", cxxopts::value<size_t>())(
      "n_factors", "factors per product state", cxxopts::value<size_t>())(
      "p,print_state", "print the state")("e,test_equality",
                                          "test equality before and after");

  options.parse_positional({"max_depth", "n_terms", "n_factors"});

  auto result = options.parse(argc, argv);

  // just error out
  size_t max_depth = result["max_depth"].as<size_t>();
  size_t n_terms = result["n_terms"].as<size_t>();
  size_t n_factors = result["n_factors"].as<size_t>();
  bool print_state = result["print_state"].as<bool>();
  bool test_equality = result["test_equality"].as<bool>();
  QSpace space = ~QSpace{0};

  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto original = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto simplified = Clone(original);

  if (print_state) {
    Print(original);
    std::cout << "\n";
  }

  simplified = Simplify(simplified);
  bool simplified_is_valid = Validate(simplified, CHECK_ALL);

  // auto squash_removed = Clone(original);
  // SquashPures(squash_removed);
  // RemoveSingles(squash_removed);
  // SquashPures(squash_removed);
  // bool squash_removed_is_valid = Validate(squash_removed, CHECK_ALL);

  bool simplified_equals_original = false;
  // bool sr_equals_original = false;
  if (test_equality) {
    simplified_equals_original = Equals_slow(simplified, original);
    // sr_equals_original = Equals_slow(squash_removed, original);
  }

  std::cout << "{\n";
  std::cout << "\t\"simplified_is_valid\": " << simplified_is_valid << ",\n";
  std::cout << "\t\"simplified_matches_original\": "
            << simplified_equals_original << "\n";
  // std::cout << "\t\"squash_removed_is_valid\": " << squash_removed_is_valid
  // << ",\n";
  // std::cout << "\t\"squash_removed_matches_original\": " <<
  // sr_equals_original
  // << "\n";
  std::cout << "}";
}
