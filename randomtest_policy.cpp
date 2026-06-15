#include "cxxopts.hpp"
#include "hamiltonian.hpp"
#include "operations.hpp"
#include "quantumstate.hpp"
#include "skewoperator.hpp"
#include "utility.hpp"
#include <complex>
#include <cstddef>
#include <iostream>
#include <nlohmann/detail/macro_scope.hpp>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
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
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(TrialArgs, max_depth, n_terms, n_factors,
                                   h_terms)
struct FastInnerMetrics {
  std::string policy;
  double real_inner;
  double imag_inner;
  double real_error;
  double imag_error;
  long time_fast;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FastInnerMetrics, policy, real_inner,
                                   imag_inner, real_error, imag_error,
                                   time_fast)
struct TrialResult {
  TrialArgs args;
  double real_inner;
  double imag_inner;
  long time_slow;
  std::vector<FastInnerMetrics> fast_metrics;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(TrialResult, args, real_inner, imag_inner,
                                   time_slow, fast_metrics)

constexpr static auto policy_naive = [](auto &bpool, const auto &,
                                        auto &kpool) -> SkewOperator {
  if (!bpool.empty()) {
    auto ret = bpool.back();
    bpool.pop_back();
    return ret;
  }
  auto ret = kpool.back();
  kpool.pop_back();
  return ret;
};

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
  if (b_curr_best > k_curr_best || (kpool.empty())) {
    SkewOperator ret1 = bpool[b_curr_best];
    bpool.erase(bpool.begin() + static_cast<std::ptrdiff_t>(b_curr_best));
    return ret1;
  }
  SkewOperator ret1 = kpool[k_curr_best];
  kpool.erase(kpool.begin() + static_cast<std::ptrdiff_t>(k_curr_best));
  return ret1;
};

std::array<std::string, 3> policy_names{"naive", "random", "overlap"};
std::array<std::function<SkewOperator(std::vector<SkewOperator> &,
                                      const SkewOperator &,
                                      std::vector<SkewOperator> &)>,
           3>
    policies{policy_naive, policy_random, policy_overlap};

auto do_test(const TrialArgs &args) -> TrialResult {
  auto [max_depth, n_terms, n_factors, h_terms] = args;
  QSpace space = ~QSpace{0};
  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto random_ket = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto random_bra = RandomSumState(max_depth, n_terms, n_factors, space, gen);

  auto H = RandomHamiltonian(h_terms, ~QSpace{0}, gen);

  auto conj = Operate(H, random_ket);

  complex inner_slow;
  auto [time_slow] =
      time_in_ms([&]() { inner_slow = Inner_slow(random_bra, conj); });

  std::vector<FastInnerMetrics> fast_metrics;

  for (const auto &&[name, policy] : std::views::zip(policy_names, policies)) {
    auto [time_fast, inner_fast_op] =
        time_in_ms([&]() { return Inner(random_bra, conj, policy); });

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

    FastInnerMetrics metrics{
        .policy = name,
        .real_inner = inner_fast.real(),
        .imag_inner = inner_fast.imag(),
        .real_error = real_error,
        .imag_error = imag_error,
        .time_fast = time_fast,
    };

    fast_metrics.push_back(metrics);
  }

  TrialResult result{.args = args,
                     .real_inner = inner_slow.real(),
                     .imag_inner = inner_slow.imag(),
                     .time_slow = time_slow,
                     .fast_metrics = fast_metrics};

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

  constexpr auto n_trials = 5;

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

  std::vector<TrialResult> trial_results;
  for (auto _ : std::views::iota(0, n_trials)) {
    auto res = do_test(args);
    trial_results.push_back(res);
  }

  nlohmann::json j{trial_results};
  std::cout << j;
}
