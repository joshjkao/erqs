#include "skewoperator.hpp"
#include "operations.hpp"
#include "optimization.hpp"
#include "quantumstate.hpp"
#include <cassert>
#include <iostream>
#include <ranges>

SkewOperator FromKet(std::shared_ptr<SumState> sum) {
  SkewOperator ret;
  std::shared_ptr<PureState> zero = MakePure(QSpace{0}, QSpace{0});
  for (const auto &[coeff, prod] : std::views::zip(sum->coeffs, sum->states)) {
    ret.AddKetBra({.coeff = coeff, .ket = prod, .bra = zero});
  }
  return ret;
}

SkewOperator FromBra(std::shared_ptr<SumState> sum) {
  SkewOperator ret;
  std::shared_ptr<PureState> zero = MakePure(QSpace{0}, QSpace{0});
  for (const auto &[coeff, prod] : std::views::zip(sum->coeffs, sum->states)) {
    ret.AddKetBra({.coeff = std::conj(coeff), .ket = zero, .bra = prod});
  }
  return ret;
}
SkewOperator FromKet(std::shared_ptr<ProductState> prod) {
  SkewOperator ret;
  ret.AddKetBra(
      {.coeff = 1, .ket = prod, .bra = MakePure(QSpace(0), QSpace(0))});
  return ret;
}

SkewOperator FromBra(std::shared_ptr<ProductState> prod) {
  SkewOperator ret;
  ret.AddKetBra(
      {.coeff = 1, .ket = MakePure(QSpace(0), QSpace(0)), .bra = prod});
  return ret;
}

void PrintToStream(std::ostream &os, const KetBra &kb) {
  os << kb.coeff << "|";
  std::visit([&](const auto e) { PrintToStream(os, e, 2); }, kb.ket);
  os << "><";
  std::visit([&](const auto e) { PrintToStream(os, e, 2); }, kb.bra);
  os << "|";
}

void PrintToStream(std::ostream &os, const SkewOperator &op) {
  for (const auto &kb : op) {
    PrintToStream(os, kb);
    os << "\n";
  }
}

void Print(const KetBra &kb) { PrintToStream(std::cout, kb); }

void Print(const SkewOperator &op) { PrintToStream(std::cout, op); }

std::ostream &operator<<(std::ostream &os, const KetBra &kb) {
  PrintToStream(os, kb);
  return os;
}

void SkewOperator::AddKetBra(const KetBra &kb) {
  ketbras.push_back(kb);
  n_kbs++;
}

void SkewOperator::AddKetBra(KetBra &&kb) {
  ketbras.emplace_back(kb);
  n_kbs++;
}

void Add(SkewOperator &o1, const SkewOperator &o2) {
  // assert(o1.GetBraSpace() == o1.GetBraSpace());
  // assert(o2.GetBraSpace() == o2.GetKetSpace());
  for (const auto &kb : o2) {
    o1.AddKetBra(kb);
  }
  o1 = Simplify(o1);
}

void Multiply(SkewOperator &o, const complex &c) {
  for (auto &kb : o) {
    kb.coeff *= c;
  }
}

SkewOperator Multiply(const SkewOperator &o1, const SkewOperator &o2) {
  SkewOperator ret{};
  for (const auto &kb1 : o1) {
    for (const auto &kb2 : o2) {
      SkewOperator o3 = Simplify(Inner(kb1.bra, kb2.ket));
      for (auto &kb3 : o3) {
        complex coeff = kb1.coeff * kb2.coeff * kb3.coeff;
        if (coeff != 0.0)
          ret.AddKetBra(
              {coeff, Tensor(kb3.ket, kb1.ket), Tensor(kb3.bra, kb2.bra)});
      }
    }
  }
  return Simplify(ret);
}
