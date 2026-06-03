#include "operations.hpp"
#include "optimization.hpp"
#include "quantumstate.hpp"
#include <cassert>
#include <complex>
#include <ranges>
#include <unordered_map>
#include <variant>

SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<PureState> &p2) {
  QSpace shared = p1->space & p2->space;
  if ((p1->bits & shared) == (p2->bits & shared)) {
    auto bra = MakePure(p1->space ^ shared, p1->bits & ~shared);
    auto ket = MakePure(p2->space ^ shared, p2->bits & ~shared);
    return SkewOperator{{KetBra{1.0, ket, bra}}};
  } else {
    auto bra = MakePure(QSpace{0}, QSpace{0});
    auto ket = MakePure(QSpace{0}, QSpace{0});
    return SkewOperator{{KetBra{0.0, ket, bra}}};
  }
}

SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<ProductState> &p2) {
  if ((GetSpace(p1) & GetSpace(p2)).none())
    return SkewOperator{{KetBra{.coeff = 1.0, .ket = p2, .bra = p1}}};
  auto bra = MakePure(QSpace{0}, QSpace{0});
  auto ket = MakePure(QSpace{0}, QSpace{0});
  SkewOperator ret{{KetBra{1.0, ket, bra}}};
  std::vector<std::shared_ptr<SumState>> unused_kets;
  std::vector<std::shared_ptr<SumState>> unused_bras;
  for (const auto &bra_sum : p1->states) {
    if ((bra_sum->space & p2->space).none()) {
      unused_bras.push_back(bra_sum);
      continue;
    }
    SkewOperator this_bra = FromBra(bra_sum);
    ret = Multiply(this_bra, ret);
    ret = Simplify(ret);
  }
  for (const auto &ket_sum : p2->states) {
    if ((ket_sum->space & p1->space).none()) {
      unused_kets.push_back(ket_sum);
      continue;
    }
    SkewOperator this_ket = FromKet(ket_sum);
    ret = Multiply(ret, this_ket);
    ret = Simplify(ret);
  }

  // unused state skipping
  if (unused_bras.empty() || unused_kets.empty()) {
    ptr_variant unused_bra_v;
    if (unused_bras.empty())
      unused_bra_v = MakePure(QSpace{0}, QSpace{0});
    else
      unused_bra_v = MakeProduct(unused_bras);
    ptr_variant unused_ket_v;
    if (unused_kets.empty())
      unused_ket_v = MakePure(QSpace{0}, QSpace{0});
    else
      unused_ket_v = MakeProduct(unused_kets);
    SkewOperator unused;
    unused.AddKetBra(
        KetBra{.coeff = 1.0, .ket = unused_ket_v, .bra = unused_bra_v});
    ret = Multiply(ret, unused);
  }
  return Simplify(ret);
}

SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<PureState> &p2) {
  std::shared_ptr<SumState> sum = MakeSum({1.0}, {p2});
  std::shared_ptr<ProductState> prod = MakeProduct({sum});
  return Inner(p1, prod);
}
SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<ProductState> &p2) {
  auto sum = MakeSum({1.0}, {p1});
  auto prod = MakeProduct({sum});
  return Inner(prod, p2);
}
SkewOperator Inner(const ptr_variant &p1, const ptr_variant &p2) {
  return std::visit([](auto &s1, auto &s2) { return Inner(s1, s2); }, p1, p2);
}
SkewOperator Inner(const std::shared_ptr<SumState> &p1,
                   const std::shared_ptr<SumState> &p2) {
  SkewOperator ret;
  for (auto &&[c1, s1] : std::views::zip(p1->coeffs, p1->states)) {
    for (auto &&[c2, s2] : std::views::zip(p2->coeffs, p2->states)) {
      SkewOperator term = Inner(s1, s2);
      Multiply(term, std::conj(c1));
      Multiply(term, c2);
      Add(ret, term);
    }
  }
  return ret;
}

static ptr_variant tensor(const std::shared_ptr<PureState> &p1,
                          const std::shared_ptr<PureState> &p2) {
  assert((p1->space & p2->space).none());
  return MakePure(QSpace{p1->space | p2->space}, QSpace{p1->bits | p2->bits});
}
static ptr_variant tensor(const std::shared_ptr<PureState> &p1,
                          const std::shared_ptr<ProductState> &p2) {
  std::vector<std::shared_ptr<SumState>> states = p2->states;
  states.push_back(MakeSum({1.0}, {p1}));
  return MakeProduct(states);
}
static ptr_variant tensor(const std::shared_ptr<ProductState> &p1,
                          const std::shared_ptr<PureState> &p2) {
  return tensor(p2, p1);
}
static ptr_variant tensor(const std::shared_ptr<ProductState> &p1,
                          const std::shared_ptr<ProductState> &p2) {
  std::vector<std::shared_ptr<SumState>> states;
  for (const auto &s1 : p1->states) {
    states.push_back(s1);
  }
  for (const auto &s2 : p2->states) {
    states.push_back(s2);
  }
  return MakeProduct(states);
}

ptr_variant Tensor(const ptr_variant &p1, const ptr_variant &p2) {
  return std::visit(
      [](const auto &e1, const auto &e2) { return tensor(e1, e2); }, p1, p2);
}
std::shared_ptr<SumState> Tensor(const std::shared_ptr<SumState> &p1,
                                 const std::shared_ptr<SumState> &p2) {
  std::vector<complex> coeffs;
  std::vector<ptr_variant> states;
  for (size_t i1 = 0; i1 < p1->states.size(); ++i1) {
    for (size_t i2 = 0; i2 < p2->states.size(); ++i2) {
      coeffs.push_back(p1->coeffs[i1] * p2->coeffs[i2]);
      states.push_back(Tensor(p1->states[i1], p2->states[i2]));
    }
  }
  return MakeSum(coeffs, states);
}
static complex operate(const PauliOperator &pauli,
                       std::shared_ptr<PureState> &state) {
  // flip x bits
  state->bits ^= ~BitString{0} & pauli.x & state->space;

  // count y overlapping bits to calculate phase
  // then flip y bits
  BitString zero_where_y = pauli.y & state->space & ~state->bits;
  auto num_zero_y = zero_where_y.count();
  BitString one_where_y = pauli.y & state->space & state->bits;
  auto num_one_y = one_where_y.count();
  state->bits ^= ~BitString{0} & pauli.y & state->space;

  // count z overlapping bits with one
  BitString one_where_z = pauli.z & state->space & state->bits;
  auto num_one_z = one_where_z.count();

  return std::pow(complex{0., 1.}, num_zero_y) *
         std::pow(complex{0., -1.}, num_one_y) *
         std::pow(complex{-1., 0.}, num_one_z);
}
static complex operate(const PauliOperator &pauli,
                       std::shared_ptr<SumState> &state);
static complex operate(const PauliOperator &pauli,
                       std::shared_ptr<ProductState> &state) {
  for (auto &p : state->states) {
    operate(pauli, p);
  }
  return 1.;
}
static complex operate(const PauliOperator &pauli,
                       std::shared_ptr<SumState> &state) {
  for (auto &&[c, s] : std::views::zip(state->coeffs, state->states)) {
    c *= std::visit([&](auto &p) { return operate(pauli, p); }, s);
  }
  return 1.;
}
std::shared_ptr<SumState> Operate(const PauliOperator &pauli,
                                  const std::shared_ptr<SumState> &state) {
  std::shared_ptr<SumState> ret = Clone(state);
  operate(pauli, ret);
  return ret;
}
std::shared_ptr<SumState> Operate(const PauliHamiltonian &h,
                                  const std::shared_ptr<SumState> &state) {
  if (h.operators.empty())
    return Clone(state);
  auto ret = Operate(h.operators[0], state);
  for (auto &c : ret->coeffs)
    c *= h.coeffs[0];
  for (auto &&[termc, term] :
       std::views::zip(h.coeffs, h.operators) | std::views::drop(1)) {
    auto r = Operate(term, state);
    for (auto &&[c, p] : std::views::zip(r->coeffs, r->states)) {
      ret->coeffs.push_back(c * termc);
      ret->states.push_back(p);
    }
  }
  return ret;
}
double Norm(const std::shared_ptr<SumState> &state) {
  SkewOperator inner = Inner(state, state);
  assert(inner.size() == 1);
  assert(GetSpace(inner[0].bra) == QSpace{0});
  assert(GetSpace(inner[0].ket) == QSpace{0});
  assert(inner[0].coeff.imag() <= 1e-6);
  return inner[0].coeff.real();
}

double ExpectedValue(const PauliOperator &pauli,
                     const std::shared_ptr<SumState> &state) {
  auto ket = Operate(pauli, state);
  SkewOperator inner = Inner(state, ket);
  assert(inner.size() == 1);
  assert(GetSpace(inner[0].bra) == QSpace{0});
  assert(GetSpace(inner[0].ket) == QSpace{0});
  return inner[0].coeff.real();
}
double ExpectedValue(const PauliHamiltonian &h,
                     const std::shared_ptr<SumState> &state) {
  auto ket = Operate(h, state);
  SkewOperator inner = Inner(state, ket);
  double self_inner = Norm(state);
  assert(inner.size() == 1 || inner.empty());
  if (inner.empty())
    return 0.0;
  assert(GetSpace(inner[0].bra) == QSpace{0});
  assert(GetSpace(inner[0].ket) == QSpace{0});
  assert(inner[0].coeff.imag() <= 1e-6);
  assert(self_inner != 0.0);
  return inner[0].coeff.real() / self_inner;
}

double ExpectedValue_slow(const PauliHamiltonian &h,
                          const std::shared_ptr<SumState> &state) {

  auto clone = Clone(state);
  Flatten(clone);

  std::unordered_map<PureState, complex, PureStateHash> conjugated_coeffs;

  double norm_sq{0.0};
  for (const auto &coeff : clone->coeffs) {
    // this is fine because Flatten combines redundant pure states
    norm_sq += (std::conj(coeff) * coeff).real();
  }

  for (const auto &&[coeff, var] :
       std::views::zip(clone->coeffs, clone->states)) {
    assert(std::holds_alternative<std::shared_ptr<PureState>>(var));
    auto pure = std::get<std::shared_ptr<PureState>>(var);
    for (const auto &&[h_coeff, op] : std::views::zip(h.coeffs, h.operators)) {
      BitString y_mask = op.y & pure->bits;
      BitString not_y_mask = op.y & ~pure->bits;
      BitString z_mask = op.z & pure->bits;
      BitString flip_mask = op.x | op.y;
      BitString new_bits = pure->bits ^ flip_mask;
      int k =
          (2 * z_mask.count() + 1 * not_y_mask.count() + 3 * y_mask.count()) %
          4;
      static const complex phases[] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
      complex phase = phases[k];
      complex new_coeff = coeff * h_coeff * phase;
      conjugated_coeffs[PureState{pure->space, new_bits}] += new_coeff;
    }
  }
  complex accumulator{0.0};
  for (const auto &&[coeff, var] :
       std::views::zip(clone->coeffs, clone->states)) {
    auto pure = std::get<std::shared_ptr<PureState>>(var);
    auto it = conjugated_coeffs.find(*pure);
    if (it == conjugated_coeffs.end())
      continue;
    accumulator += std::conj(coeff) * it->second;
  }
  return accumulator.real() / norm_sq;
}
