#pragma once
#include "common.hpp"
#include <vector>

struct KetBra {
  complex coeff;
  ptr_variant ket;
  ptr_variant bra;
};

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

class SkewOperator {
private:
  std::vector<KetBra> ketbras;

public:
  SkewOperator() = default;
  SkewOperator(std::vector<KetBra> kbs) : ketbras(kbs) {}

  void AddKetBra(const KetBra &kb) { ketbras.push_back(kb); }
  void AddKetBra(KetBra &&kb) { ketbras.emplace_back(kb); }

  auto begin() const { return ketbras.begin(); }
  auto begin() { return ketbras.begin(); }
  auto end() const { return ketbras.end(); }
  auto end() { return ketbras.end(); }

  auto size() const { return ketbras.size(); }
  auto empty() const { return ketbras.empty(); }

  auto operator[](size_t index) const { return ketbras[index]; }
  auto operator[](size_t index) -> KetBra & { return ketbras[index]; }
};

SkewOperator FromKet(std::shared_ptr<SumState> sum);
SkewOperator FromBra(std::shared_ptr<SumState> sum);

void PrintToStream(std::ostream &os, const KetBra &kb);
void PrintToStream(std::ostream &os, const SkewOperator &op);
void Print(const KetBra &kb);
void Print(const SkewOperator &op);
std::ostream &operator<<(std::ostream &os, const KetBra &kb);
