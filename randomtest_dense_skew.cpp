#include "cxxopts.hpp"
#include "operations.h"
#include "quantumstate.h"
#include "validation.h"
#include <iostream>
#include <pthread.h>
#include <random>
#include <ranges>
#include <stdexcept>

struct TrialArgs {
  size_t max_depth;
  size_t n_terms;
  size_t n_factors;
  size_t h_terms;
};

struct TrialResult {
  complex inner_fast;
  complex inner_slow;
  double real_error;
  double imag_error;
  long time_fast;
  long time_slow;
};

long time_in_ms(auto func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  return duration.count();
}

auto do_test(const TrialArgs &args) -> TrialResult {
  auto [max_depth, n_terms, n_factors, h_terms] = args;
  QSpace space = ~QSpace{0};
  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto random_ket = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto random_bra = RandomSumState(max_depth, n_terms, n_factors, space, gen);

  auto random_bitset = [&]() -> BitString {
    std::bernoulli_distribution dist{0.5};
    BitString bits;
    for (size_t i = 0; i < NQUBITS; ++i) {
      bits[i] = dist(gen);
    }
    return bits;
  };

  std::vector<double> coeffs;
  std::vector<PauliOperator> ops;

  for (auto _ : std::views::iota(0uz, h_terms)) {
    BitString x = random_bitset();
    BitString y = random_bitset() & ~x;
    BitString z = random_bitset() & ~x & ~y;
    coeffs.push_back(1);
    ops.emplace_back(x, y, z);
  }

  PauliHamiltonian H{coeffs, ops};
  if (!Validate(H, CHECK_ALL)) {
    throw std::runtime_error("invalid hamiltonian");
  }

  auto conj = Operate(H, random_ket);
  auto conj_clone = Clone(conj);
  // auto inner_fast_op = Inner(random_bra, conj);
  SkewOperator inner_fast_op;
  auto time_fast =
      time_in_ms([&]() { inner_fast_op = Inner(random_bra, conj); });

  // auto inner_slow = Inner_slow(random_bra, conj_clone);
  complex inner_slow;
  auto time_slow =
      time_in_ms([&]() { inner_slow = Inner_slow(random_bra, conj_clone); });

  complex inner_fast;
  if (inner_fast_op.ketbras.size() > 2)
    throw std::runtime_error("skew op isn't a single term");
  else if (inner_fast_op.ketbras.empty())
    inner_fast = 0;
  else if (GetSpace(inner_fast_op.ketbras[0].ket) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else if (GetSpace(inner_fast_op.ketbras[0].bra) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else
    inner_fast = inner_fast_op.ketbras[0].coeff;

  double real_error =
      std::abs((inner_slow.real() - inner_fast.real()) / inner_slow.real());
  double imag_error =
      std::abs((inner_slow.imag() - inner_fast.imag()) / inner_slow.imag());

  TrialResult result{.inner_fast = inner_fast,
                     .inner_slow = inner_slow,
                     .real_error = real_error,
                     .imag_error = imag_error,
                     .time_fast = time_fast,
                     .time_slow = time_slow};
  return result;
}

auto main(int argc, char **argv) -> int {
  cxxopts::Options options("randomtest_expectation",
                           "attempt to validate expectation values");
  options.add_options()("max_depth", "maximum tree depth",
                        cxxopts::value<size_t>())(
      "n_terms", "terms per sum state", cxxopts::value<size_t>())(
      "n_factors", "factors per product state", cxxopts::value<size_t>())(
      "h_terms", "number of terms to make the hamiltonian",
      cxxopts::value<size_t>())("p,print_state", "print the state");

  constexpr auto n_trials = 1000;

  options.parse_positional({"max_depth", "n_terms", "n_factors", "h_terms"});

  auto result = options.parse(argc, argv);

  // just error out
  size_t max_depth = result["max_depth"].as<size_t>();
  size_t n_terms = result["n_terms"].as<size_t>();
  size_t n_factors = result["n_factors"].as<size_t>();
  size_t h_terms = result["h_terms"].as<size_t>();
  // bool print_state = result["print_state"].as<bool>();

  TrialArgs args{.max_depth = max_depth,
                 .n_terms = n_terms,
                 .n_factors = n_factors,
                 .h_terms = h_terms};

#pragma omp parallel for
  for (auto _ : std::views::iota(0, n_trials)) {
    auto res = do_test(args);
#pragma omp critical
    {
      std::cout << res.real_error << " " << res.imag_error << " "
                << res.time_fast << " " << res.time_slow << "\n";
    }
  }
}
