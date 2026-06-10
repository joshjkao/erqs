#include "cxxopts.hpp"
#include "operations.hpp"
#include "quantumstate.hpp"
#include "skewoperator.hpp"
#include "validation.hpp"
#include <cstddef>
#include <iostream>
#include <pthread.h>
#include <random>
#include <ranges>
#include <stdexcept>
#include <type_traits>

struct TrialArgs {
  size_t max_depth;
  size_t n_terms;
  size_t n_factors;
  size_t h_terms;
};

struct TrialResult {
  complex inner_fast_1;
  complex inner_fast_2;
  complex inner_slow;
  double real_error_1;
  double imag_error_1;
  long time_fast_1;
  double real_error_2;
  double imag_error_2;
  long time_fast_2;
  long time_slow;
};

// constexpr static auto policy_naive = [](auto &bpool, const auto &,
//                                         auto &kpool) -> SkewOperator {
//   if (!bpool.empty()) {
//     auto ret = bpool.back();
//     bpool.pop_back();
//     return ret;
//   }
//   auto ret = kpool.back();
//   kpool.pop_back();
//   return ret;
// };

constexpr static auto policy_random = [](auto &bpool, const auto &,
                                         auto &kpool) -> SkewOperator {
  assert(!bpool.empty() || !kpool.empty());
  static thread_local std::random_device rand{};
  static thread_local std::mt19937 gen{rand()};
  static thread_local std::uniform_int_distribution<> dist(0, 1);
  std::ranges::shuffle(bpool, gen);
  std::ranges::shuffle(kpool, gen);
  if (bpool.empty()) {
    auto ret = kpool.back();
    kpool.pop_back();
    return ret;
  } else if (kpool.empty()) {
    auto ret = bpool.back();
    bpool.pop_back();
    return ret;
  }
  if (dist(gen)) {
    auto ret = kpool.back();
    kpool.pop_back();
    return ret;
  }
  auto ret = bpool.back();
  bpool.pop_back();
  return ret;
};

constexpr static auto policy_overlap = [](auto &bpool, const auto &ret,
                                          auto &kpool) -> SkewOperator {
  QSpace bspace = ret.GetBraSpace();
  QSpace kspace = ret.GetKetSpace();
  size_t b_max_overlap{0uz};
  size_t b_curr_best{0uz};
  for (auto i{0uz}; i < bpool.size(); ++i) {
    size_t overlap = (bpool[i].GetBraSpace() & bspace).count();
    if (overlap > b_max_overlap) {
      b_max_overlap = overlap;
      b_curr_best = i;
    }
  }
  size_t k_max_overlap{0uz};
  size_t k_curr_best{0uz};
  for (auto i{0uz}; i < kpool.size(); ++i) {
    size_t overlap = (kpool[i].GetKetSpace() & kspace).count();
    if (overlap > k_max_overlap) {
      k_max_overlap = overlap;
      k_curr_best = i;
    }
  }
  if (b_curr_best > k_curr_best || kpool.empty()) {
    SkewOperator ret1 = bpool[b_curr_best];
    bpool.erase(bpool.begin() + static_cast<std::ptrdiff_t>(b_curr_best));
    return ret1;
  }
  SkewOperator ret1 = kpool[k_curr_best];
  kpool.erase(kpool.begin() + static_cast<std::ptrdiff_t>(k_curr_best));
  return ret1;
};

template <typename ResultType> struct TimedResult {
  long time;
  ResultType ret;
};
template <> struct TimedResult<void> {
  long time;
};
auto time_in_ms(auto &&func) {
  using ResultType = std::invoke_result_t<decltype(func)>;
  auto start = std::chrono::high_resolution_clock::now();
  if constexpr (std::is_void_v<ResultType>) {
    std::invoke(std::forward<decltype(func)>(func));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return TimedResult<ResultType>{.time = duration.count()};
  } else {
    auto ret = std::invoke(std::forward<decltype(func)>(func));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return TimedResult<ResultType>{.time = duration.count(),
                                   .ret = std::move(ret)};
  }
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

  auto [time_fast_1, inner_fast_op_1] =
      time_in_ms([&]() { return Inner(random_bra, conj, policy_random); });

  auto [time_fast_2, inner_fast_op_2] =
      time_in_ms([&]() { return Inner(random_bra, conj, policy_overlap); });

  complex inner_slow;
  auto [time_slow] =
      time_in_ms([&]() { inner_slow = Inner_slow(random_bra, conj_clone); });

  complex inner_fast_1;
  if (inner_fast_op_1.size() > 2)
    throw std::runtime_error("skew op isn't a single term");
  else if (inner_fast_op_1.empty())
    inner_fast_1 = 0;
  else if (GetSpace(inner_fast_op_1[0].ket) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else if (GetSpace(inner_fast_op_1[0].bra) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else
    inner_fast_1 = inner_fast_op_1[0].coeff;

  complex inner_fast_2;
  if (inner_fast_op_2.size() > 2)
    throw std::runtime_error("skew op isn't a single term");
  else if (inner_fast_op_2.empty())
    inner_fast_2 = 0;
  else if (GetSpace(inner_fast_op_2[0].ket) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else if (GetSpace(inner_fast_op_2[0].bra) != QSpace{0})
    throw std::runtime_error("skewop isn't a number");
  else
    inner_fast_2 = inner_fast_op_2[0].coeff;

  double real_error_1 =
      std::abs((inner_slow.real() - inner_fast_1.real()) / inner_slow.real());
  double imag_error_1 =
      std::abs((inner_slow.imag() - inner_fast_1.imag()) / inner_slow.imag());

  double real_error_2 =
      std::abs((inner_slow.real() - inner_fast_2.real()) / inner_slow.real());
  double imag_error_2 =
      std::abs((inner_slow.imag() - inner_fast_2.imag()) / inner_slow.imag());

  TrialResult result{.inner_fast_1 = inner_fast_1,
                     .inner_fast_2 = inner_fast_2,
                     .inner_slow = inner_slow,
                     .real_error_1 = real_error_1,
                     .imag_error_1 = imag_error_1,
                     .time_fast_1 = time_fast_1,
                     .real_error_2 = real_error_2,
                     .imag_error_2 = imag_error_2,
                     .time_fast_2 = time_fast_2,
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

  constexpr auto n_trials = 100;

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
      std::cout << res.real_error_1 << " " << res.imag_error_1 << " "
                << res.real_error_2 << " " << res.imag_error_2 << " "
                << res.time_fast_1 << " " << res.time_fast_2 << " "
                << res.time_slow << "\n";
    }
  }
}
