#include "skewoperator.hpp"
#include "quantumstate.hpp"
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
