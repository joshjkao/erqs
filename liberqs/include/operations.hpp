#pragma once
#include "common.hpp"
#include "hamiltonian.hpp"
#include "optimization.hpp"
#include "skewoperator.hpp"
#include <cassert>
#include <concepts>
#include <ranges>
#include <type_traits>

// ---- INNER PRODUCTS ---- //
SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<PureState> &p2);
SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<ProductState> &p2);
SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<PureState> &p2);
SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<ProductState> &p2);
SkewOperator Inner(const ptr_variant &p1, const ptr_variant &p2);
SkewOperator Inner(const std::shared_ptr<SumState> &p1,
                   const std::shared_ptr<SumState> &p2);

// ---- INNER PRODUCTS ---- //
SkewOperator Inner_double(const std::shared_ptr<PureState> &p1,
                          const std::shared_ptr<PureState> &p2);
SkewOperator Inner_double(const std::shared_ptr<ProductState> &p1,
                          const std::shared_ptr<ProductState> &p2);
SkewOperator Inner_double(const std::shared_ptr<ProductState> &p1,
                          const std::shared_ptr<PureState> &p2);
SkewOperator Inner_double(const std::shared_ptr<PureState> &p1,
                          const std::shared_ptr<ProductState> &p2);
SkewOperator Inner_double(const ptr_variant &p1, const ptr_variant &p2);
SkewOperator Inner_double(const std::shared_ptr<SumState> &p1,
                          const std::shared_ptr<SumState> &p2);

using skewop_pool = std::vector<SkewOperator>;
template <typename Func>
concept OrderPolicy =
    std::invocable<Func, skewop_pool &, const SkewOperator &, skewop_pool &> &&
    std::convertible_to<
        std::invoke_result_t<Func, skewop_pool &, const SkewOperator &,
                             skewop_pool &>,
        SkewOperator>;
SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<ProductState> &p2,
                   OrderPolicy auto policy);

ptr_variant Tensor(const ptr_variant &p1, const ptr_variant &p2);
std::shared_ptr<SumState> Tensor(const std::shared_ptr<SumState> &p1,
                                 const std::shared_ptr<SumState> &p2);

// ---- OPERATIONS ---- //
std::shared_ptr<SumState> Operate(const PauliOperator &pauli,
                                  const std::shared_ptr<SumState> &state);
std::shared_ptr<SumState> Operate(const PauliHamiltonian &h,
                                  const std::shared_ptr<SumState> &state);
double Norm(const std::shared_ptr<SumState> &state);
double ExpectedValue(const PauliOperator &pauli,
                     const std::shared_ptr<SumState> &state);
double ExpectedValue(const PauliHamiltonian &h,
                     const std::shared_ptr<SumState> &state);
double ExpectedValue_slow(const PauliHamiltonian &h,
                          const std::shared_ptr<SumState> &state);

// templated implementation for arbitrary contraction orderings
SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<PureState> &p2, OrderPolicy auto) {
  return Inner(p1, p2);
}
SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<ProductState> &p2,
                   OrderPolicy auto policy) {
  if ((GetSpace(p1) & GetSpace(p2)).none())
    return SkewOperator{{KetBra{.coeff = 1.0, .ket = p2, .bra = p1}}};
  // QSpace shared_space = GetSpace(p1) & GetSpace(p2);
  auto bra = MakePure(QSpace{0}, QSpace{0});
  auto ket = MakePure(QSpace{0}, QSpace{0});
  SkewOperator ret{{KetBra{1.0, ket, bra}}};
  std::vector<std::shared_ptr<SumState>> unused_kets;
  std::vector<std::shared_ptr<SumState>> unused_bras;
  std::vector<SkewOperator> bra_pool;
  std::vector<SkewOperator> ket_pool;
  for (const auto &bra_sum : p1->states) {
    if ((bra_sum->space & p2->space).none()) {
      unused_bras.push_back(bra_sum);
      continue;
    }
    SkewOperator this_bra = FromBra(bra_sum);
    bra_pool.push_back(this_bra);
  }
  for (const auto &ket_sum : p2->states) {
    if ((ket_sum->space & p1->space).none()) {
      unused_kets.push_back(ket_sum);
      continue;
    }
    SkewOperator this_ket = FromKet(ket_sum);
    ket_pool.push_back(this_ket);
  }
  size_t bra_pool_size = bra_pool.size();
  size_t ket_pool_size = ket_pool.size();
  while (!bra_pool.empty() || !ket_pool.empty()) {
    SkewOperator op = policy(bra_pool, ret, ket_pool);
    if (bra_pool.size() != bra_pool_size) {
      ret = Multiply(op, ret);
      ret = Simplify(ret);
      bra_pool_size = bra_pool.size();
    } else if (ket_pool.size() != ket_pool_size) {
      ret = Multiply(ret, op);
      ret = Simplify(ret);
      ket_pool_size = ket_pool.size();
    } else {
      assert(false);
    }
  }
  if (!unused_bras.empty()) {
    auto unused_bras_p = MakeProduct(unused_bras);
    auto unused_bras_op = FromBra(unused_bras_p);
    ret = Multiply(unused_bras_op, ret);
    ret = Simplify(ret);
  }
  if (!unused_kets.empty()) {
    auto unused_kets_p = MakeProduct(unused_kets);
    auto unused_kets_op = FromKet(unused_kets_p);
    ret = Multiply(ret, unused_kets_op);
    ret = Simplify(ret);
  }
  // assert((shared_space & ret.GetBraSpace()).none());
  // assert((shared_space & ret.GetKetSpace()).none());
  return ret;
}
SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<PureState> &p2,
                   OrderPolicy auto policy) {
  std::shared_ptr<SumState> sum = MakeSum({1.0}, {p2});
  std::shared_ptr<ProductState> prod = MakeProduct({sum});
  return Inner(p1, prod, policy);
}
SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<ProductState> &p2,
                   OrderPolicy auto policy) {
  auto sum = MakeSum({1.0}, {p1});
  auto prod = MakeProduct({sum});
  return Inner(prod, p2, policy);
}
SkewOperator Inner(const ptr_variant &p1, const ptr_variant &p2,
                   OrderPolicy auto policy) {
  return std::visit([&](auto &s1, auto &s2) { return Inner(s1, s2, policy); },
                    p1, p2);
}
SkewOperator Inner(const std::shared_ptr<SumState> &p1,
                   const std::shared_ptr<SumState> &p2,
                   OrderPolicy auto policy) {
  SkewOperator ret;
  for (auto &&[c1, s1] : std::views::zip(p1->coeffs, p1->states)) {
    for (auto &&[c2, s2] : std::views::zip(p2->coeffs, p2->states)) {
      SkewOperator term = Inner(s1, s2, policy);
      Multiply(term, std::conj(c1));
      Multiply(term, c2);
      Add(ret, term);
    }
  }
  return ret;
}
