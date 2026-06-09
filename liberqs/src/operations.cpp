#include "operations.h"
#include "optimization.h"
#include "quantumstate.h"
#include <cassert>
#include <complex>
#include <iostream>
#include <optional>
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

  QSpace unused_bra = p1->space;
  QSpace unused_ket = p2->space;
  SkewOperator ret;
  ret.ketbras.push_back(KetBra{.coeff = 1.0,
                               .ket = MakePure(QSpace{0}, QSpace{0}),
                               .bra = MakePure(QSpace{0}, QSpace{0})});
  for (auto &s1 : p1->states) {
    for (auto &s2 : p2->states) {
      if ((s1->space & unused_bra & s2->space & unused_ket).any()) {
        unused_bra ^= s1->space;
        unused_ket ^= s2->space;
        SkewOperator term = Inner(s1, s2);
        ret = Multiply(ret, term);
      }
    }
  }
  std::vector<std::shared_ptr<SumState>> unused_bras;
  std::vector<std::shared_ptr<SumState>> unused_kets;
  bool has_unused = false;
  for (auto &s1 : p1->states) {
    if ((s1->space & unused_bra).any()) {
      unused_bras.push_back(s1);
      has_unused = true;
    }
  }
  for (auto &s2 : p2->states) {
    if ((s2->space & unused_ket).any()) {
      unused_kets.push_back(s2);
      has_unused = true;
    }
  }
  if (has_unused) {
    ptr_variant ket;
    if (unused_kets.empty())
      ket = MakePure(QSpace{0}, QSpace{0});
    else
      ket = MakeProduct(unused_kets);
    ptr_variant bra;
    if (unused_bras.empty())
      bra = MakePure(QSpace{0}, QSpace{0});
    else
      bra = MakeProduct(unused_bras);
    KetBra unused{.coeff = 1.0, .ket = ket, .bra = bra};
    SkewOperator extra{{unused}};
    ret = Multiply(ret, extra);
  }
  return ret;
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

// needed to perform compression
using PureKB = std::tuple<BitString, BitString, BitString, BitString>;
struct KBTupleHash {
  std::size_t operator()(const PureKB &t) const {
    std::size_t seed = 0;

    // Boost's hash_combine algorithm
    auto hash_combine = [&seed](std::size_t hash_value) {
      seed ^= hash_value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    // 2. Hash the four std::bitset<N> elements
    std::hash<BitString> bitset_hasher;
    hash_combine(bitset_hasher(std::get<0>(t)));
    hash_combine(bitset_hasher(std::get<1>(t)));
    hash_combine(bitset_hasher(std::get<2>(t)));
    hash_combine(bitset_hasher(std::get<3>(t)));

    return seed;
  }
};

using pkb_result =
    std::tuple<complex, BitString, BitString, BitString, BitString>;
static std::optional<pkb_result> get_pkb_if_pure(const KetBra &kb) {
  auto ket_ptr_ptr = std::get_if<std::shared_ptr<PureState>>(&kb.ket);
  auto bra_ptr_ptr = std::get_if<std::shared_ptr<PureState>>(&kb.bra);
  if (ket_ptr_ptr && bra_ptr_ptr) {
    pkb_result res{kb.coeff, (*ket_ptr_ptr)->space, (*ket_ptr_ptr)->bits,
                   (*bra_ptr_ptr)->space, (*bra_ptr_ptr)->bits};
    return res;
  }
  return std::nullopt;
}

SkewOperator Simplify(const SkewOperator &op) {
  SkewOperator ret;
  std::unordered_map<PureKB, complex, KBTupleHash> pures;
  for (const auto &kb : op.ketbras) {
    auto pure_opt = get_pkb_if_pure(kb);
    if (pure_opt) {
      auto &[coeff, kspace, kbits, bspace, bbits] = pure_opt.value();
      PureKB pkb{kspace, kbits, bspace, bbits};
      auto [iterator, inserted] = pures.insert({pkb, coeff});
      if (!inserted) {
        iterator->second += coeff;
      }
    } else {
      ret.ketbras.push_back(kb);
    }
  }
  for (const auto &[pkb, coeff] : pures) {
    auto &[kspace, kbits, bspace, bbits] = pkb;
    auto ket = MakePure(kspace, kbits);
    auto bra = MakePure(bspace, bbits);
    KetBra kb{.coeff = coeff, .ket = ket, .bra = bra};
    ret.ketbras.push_back(kb);
  }
  for (auto &[coeff, ket, bra] : ret.ketbras) {
    Simplify(ket);
    Simplify(bra);
  }
  return ret;
}

SkewOperator CompressConstants(const SkewOperator &op) {
  SkewOperator ret;
  KetBra constants{
      .coeff = 0.0, .ket = MakePure(0ull, 0ull), .bra = MakePure(0ull, 0ull)};
  for (const auto &kb : op.ketbras) {
    auto pure_opt = get_pkb_if_pure(kb);
    if (!pure_opt) {
      ret.ketbras.push_back(kb);
      continue;
    }
    auto &[coeff, kspace, kbits, bspace, bbits] = pure_opt.value();
    if (kspace.any() || bspace.any()) {
      ret.ketbras.push_back(kb);
      continue;
    }
    constants.coeff += coeff;
  }
  ret.ketbras.push_back(constants);
  return ret;
}

void Add(SkewOperator &o1, const SkewOperator &o2) {
  for (const auto &kb : o2.ketbras) {
    o1.ketbras.push_back(kb);
  }
  o1 = Simplify(o1);
}
void Multiply(SkewOperator &o, const complex &c) {
  for (auto &kb : o.ketbras) {
    kb.coeff *= c;
  }
}

SkewOperator Multiply(const SkewOperator &o1, const SkewOperator &o2) {
  SkewOperator ret{};
  for (const auto &kb1 : o1.ketbras) {
    for (const auto &kb2 : o2.ketbras) {
      SkewOperator o3 = Simplify(Inner(kb1.bra, kb2.ket));
      SkewOperator o4 = Simplify(Inner(kb2.bra, kb1.ket));
      for (auto &kb3 : o3.ketbras) {
        for (auto &kb4 : o4.ketbras) {
          complex coeff = kb1.coeff * kb2.coeff * kb3.coeff * kb4.coeff;
          if (coeff != 0.0)
            ret.ketbras.push_back(
                {kb1.coeff * kb2.coeff * kb3.coeff * kb4.coeff,
                 Tensor(kb3.ket, kb4.ket), Tensor(kb3.bra, kb4.bra)});
        }
      }
    }
  }
  return Simplify(ret);
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
  assert(inner.ketbras.size() == 1);
  assert(GetSpace(inner.ketbras[0].bra) == QSpace{0});
  assert(GetSpace(inner.ketbras[0].ket) == QSpace{0});
  assert(inner.ketbras[0].coeff.imag() <= 1e-6);
  return inner.ketbras[0].coeff.real();
}

double ExpectedValue(const PauliOperator &pauli,
                     const std::shared_ptr<SumState> &state) {
  auto ket = Operate(pauli, state);
  SkewOperator inner = Inner(state, ket);
  assert(inner.ketbras.size() == 1);
  assert(GetSpace(inner.ketbras[0].bra) == QSpace{0});
  assert(GetSpace(inner.ketbras[0].ket) == QSpace{0});
  return inner.ketbras[0].coeff.real();
}
double ExpectedValue(const PauliHamiltonian &h,
                     const std::shared_ptr<SumState> &state) {
  auto ket = Operate(h, state);
  SkewOperator inner = Inner(state, ket);
  double self_inner = Norm(state);
  assert(inner.ketbras.size() == 1 || inner.ketbras.empty());
  if (inner.ketbras.empty())
    return 0.0;
  assert(GetSpace(inner.ketbras[0].bra) == QSpace{0});
  assert(GetSpace(inner.ketbras[0].ket) == QSpace{0});
  assert(inner.ketbras[0].coeff.imag() <= 1e-6);
  assert(self_inner != 0.0);
  return inner.ketbras[0].coeff.real() / self_inner;
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
