#include "cxxopts.hpp"
#include "hamiltonian.hpp"
#include "operations.hpp"
#include "quantumstate.hpp"
#include "skewoperator.hpp"
#include "utility.hpp"
#include <atomic>
#include <complex>
#include <cstddef>
#include <expected>
#include <iostream>
#include <memory>
#include <pthread.h>
#include <random>
#include <stdexcept>

struct TrialArgs {
  size_t max_depth;
  size_t n_terms;
  size_t n_factors;
  size_t h_terms;
};

struct ImplMetrics {
  std::string policy;
  double real_inner;
  double imag_inner;
  long time;
};
struct TrialResult {
  TrialArgs args;
  std::vector<ImplMetrics> metrics;
};
auto parse_skewop_ret(const SkewOperator &op)
    -> std::expected<complex, std::string> {
  complex ret;
  if (op.size() > 2)
    return std::unexpected<std::string>{"skewop isn't a single term"};
  else if (op.empty())
    ret = 0;
  else if (GetSpace(op[0].ket) != QSpace{0})
    return std::unexpected<std::string>{"skewop isn't a number"};
  else if (GetSpace(op[0].bra) != QSpace{0})
    return std::unexpected<std::string>{"skewop isn't a number"};
  else
    ret = op[0].coeff;
  return ret;
};

// constexpr static auto inner_slow =
//     [](std::shared_ptr<SumState> bra,
//        std::shared_ptr<SumState> ket) -> std::expected<complex, std::string>
//        {
//   return Inner_slow(bra, ket);
// };

constexpr static auto inner_double =
    [](std::shared_ptr<SumState> bra,
       std::shared_ptr<SumState> ket) -> std::expected<complex, std::string> {
  auto inner = Inner_double(bra, ket);
  return parse_skewop_ret(inner);
};

// constexpr static auto inner_naive =
//     [](std::shared_ptr<SumState> bra,
//        std::shared_ptr<SumState> ket) -> std::expected<complex, std::string>
//        {
//   constexpr static auto policy_naive = [](auto &bpool, const auto &,
//                                           auto &kpool) -> SkewOperator {
//     if (!bpool.empty()) {
//       auto ret = bpool.back();
//       bpool.pop_back();
//       return ret;
//     }
//     auto ret = kpool.back();
//     kpool.pop_back();
//     return ret;
//   };
//   auto inner = Inner(bra, ket, policy_naive);
//   return parse_skewop_ret(inner);
// };

constexpr static auto inner_random =
    [](std::shared_ptr<SumState> bra,
       std::shared_ptr<SumState> ket) -> std::expected<complex, std::string> {
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
  auto inner = Inner(bra, ket, policy_random);
  return parse_skewop_ret(inner);
};

constexpr static auto inner_overlap =
    [](std::shared_ptr<SumState> bra,
       std::shared_ptr<SumState> ket) -> std::expected<complex, std::string> {
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
  auto inner = Inner(bra, ket, policy_overlap);
  return parse_skewop_ret(inner);
};

struct InnerImpl {
  std::string name;
  std::function<std::expected<complex, std::string>(std::shared_ptr<SumState>,
                                                    std::shared_ptr<SumState>)>
      impl;
};

std::array<InnerImpl, 3> inner_impls{{{"doublecon", inner_double},
                                      {"random", inner_random},
                                      {"overlap", inner_overlap}}};

auto do_test(const TrialArgs &args) -> TrialResult {
  auto [max_depth, n_terms, n_factors, h_terms] = args;
  QSpace space = ~QSpace{0};
  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto random_ket = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto random_bra = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto H = RandomHamiltonian(h_terms, ~QSpace{0}, gen);
  auto conj = Operate(H, random_ket);
  std::vector<ImplMetrics> metrics;
  for (const auto &[name, inner_impl] : inner_impls) {
    auto [time, inner_op] =
        time_in_ms([&]() { return inner_impl(random_bra, conj); });

    auto inner = inner_op
                     .or_else([](const auto &err) -> decltype(inner_op) {
                       throw std::runtime_error(err);
                     })
                     .value();
    metrics.emplace_back(name, inner.real(), inner.imag(), time);
  }

  return {args, metrics};
}

auto main(int argc, char **argv) -> int {
  cxxopts::Options options("randomtest_expectation",
                           "attempt to validate expectation values");
  options.add_options()("max_depth", "maximum tree depth",
                        cxxopts::value<size_t>())(
      "n_terms", "terms per sum state", cxxopts::value<size_t>())(
      "n_factors", "factors per product state", cxxopts::value<size_t>())(
      "h_terms", "number of terms to make the hamiltonian",
      cxxopts::value<size_t>())("p,print_state", "print the state")(
      "s,skip_slow", "skip the slow inner product");

  constexpr auto n_trials = 100;

  options.parse_positional({"max_depth", "n_terms", "n_factors", "h_terms"});

  auto result = options.parse(argc, argv);

  // just error out
  size_t max_depth = result["max_depth"].as<size_t>();
  size_t n_terms = result["n_terms"].as<size_t>();
  size_t n_factors = result["n_factors"].as<size_t>();
  size_t h_terms = result["h_terms"].as<size_t>();
  // bool print_state = result["print_state"].as<bool>();
  // bool skip_slow = result["skip_slow"].as<bool>();

  TrialArgs args{.max_depth = max_depth,
                 .n_terms = n_terms,
                 .n_factors = n_factors,
                 .h_terms = h_terms};

  std::vector<TrialResult> trial_results(n_trials);

  std::atomic<size_t> completed_trials{0};

  size_t update_interval = n_trials / 20;
  if (update_interval == 0)
    update_interval = 1; // Prevent modulo by zero

#pragma omp parallel for
  for (size_t i = 0; i < n_trials; ++i) {
    auto res = do_test(args);

    for (size_t j = 0; j < inner_impls.size(); ++j) {
      trial_results[i] = res;
    }

    size_t current = ++completed_trials;

    if (current % update_interval == 0 || current == n_trials) {

#pragma omp critical
      {
        std::cout << "Progress: " << current << " / " << n_trials << " ("
                  << (current * 100 / n_trials) << "%)\r" << std::flush;
      }
    }
  }
  std::cout << std::endl;

  std::cout << "trial,";
  for (const auto &[name, _] : inner_impls) {
    std::cout << name << "_real_inner," << name << "_imag_inner," << name
              << "_time,";
  }
  std::cout << std::endl;
  for (const auto &trial_result : trial_results) {
    auto [res_args, res_metrics] = trial_result;
    for (const auto &metric : res_metrics) {
      std::cout << metric.real_inner << "," << metric.imag_inner << ","
                << metric.time << ",";
    }
    std::cout << "\n";
  }
}
