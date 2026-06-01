#include "cxxopts.hpp"
#include "operations.h"
#include "quantumstate.h"
#include "validation.h"
#include <iostream>
#include <random>

struct TrialArgs {
  size_t max_depth;
  size_t n_terms;
  size_t n_factors;
};

struct TrialResult {
  complex inner_fast;
  complex inner_slow;
  double real_error;
  double imag_error;
};

auto do_test(const TrialArgs &args) -> TrialResult {
  auto [max_depth, n_terms, n_factors] = args;
  QSpace space = ~QSpace{0};

  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto ket = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto bra = RandomSumState(max_depth, n_terms, n_factors, space, gen);

  auto inner_fast_op = Inner(bra, ket);
  auto inner_slow = Inner_slow(bra, ket);

  Validate(inner_fast_op, CHECK_ALL);

  complex inner_fast;
  if (inner_fast_op.size() > 2)
    throw std::runtime_error("skew op isn't a single term");
  else if (inner_fast_op.empty())
    inner_fast = 0;
  else if (GetSpace(inner_fast_op[0].ket) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else if (GetSpace(inner_fast_op[0].bra) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else
    inner_fast = inner_fast_op[0].coeff;

  double real_error =
      std::abs((inner_slow.real() - inner_fast.real()) / inner_slow.real());
  double imag_error =
      std::abs((inner_slow.imag() - inner_fast.imag()) / inner_slow.imag());

  TrialResult result{.inner_fast = inner_fast,
                     .inner_slow = inner_slow,
                     .real_error = real_error,
                     .imag_error = imag_error};

  return result;
}
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

  TrialArgs args{
      .max_depth = max_depth, .n_terms = n_terms, .n_factors = n_factors};

  auto res = do_test(args);
  std::cout << res.real_error << " " << res.imag_error << "\n";
}
