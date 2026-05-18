#pragma once
#include <bitset>
#include <complex>
#include <functional>
#include <memory>
#include <random>
#include <variant>
#include <vector>

constexpr size_t NQUBITS = 16;

// ---- BASIC TYPES ---- //
using complex = std::complex<double>;
using QSpace = std::bitset<NQUBITS>;
using BitString = std::bitset<NQUBITS>;

struct PureState;
struct ProductState;
struct SumState;

struct PureState {
  QSpace space;
  BitString bits;
};

struct ProductState {
  QSpace space;
  std::vector<std::shared_ptr<SumState>> states;
};

using ptr_variant =
    std::variant<std::shared_ptr<PureState>, std::shared_ptr<ProductState>>;
struct SumState {
  QSpace space;
  std::vector<complex> coeffs;
  std::vector<ptr_variant> states;
};

// ---- INNER PRODUCT TYPES ---- //
struct KetBra {
  complex coeff;
  ptr_variant ket;
  ptr_variant bra;
};

struct SkewOperator {
  std::vector<KetBra> ketbras;
};

struct PauliOperator {
  // assume these are always disjoint
  // a 1 means that gate acts on that bit
  BitString x;
  BitString y;
  BitString z;
};

struct PauliHamiltonian {
  std::vector<double> coeffs;
  std::vector<PauliOperator> operators;
};

// ---- FACTORIES ---- //
std::shared_ptr<PureState> MakePure(auto space, auto bits) {
  return std::make_shared<PureState>(QSpace(space), BitString(bits));
}
std::shared_ptr<SumState> MakeSum(const std::vector<complex> &coeffs,
                                  const std::vector<ptr_variant> &states);
std::shared_ptr<ProductState>
MakeProduct(const std::vector<std::shared_ptr<SumState>> &states);

std::shared_ptr<SumState> Clone(const std::shared_ptr<SumState> &root);

ptr_variant RandomProductState(size_t max_depth, size_t n_terms,
                               size_t n_factors, QSpace space,
                               std::mt19937 &gen);

std::shared_ptr<SumState> RandomSumState(size_t max_depth, size_t n_terms,
                                         size_t n_factors, QSpace space,
                                         std::mt19937 &gen);

// ---- GETTERS ---- //
inline const QSpace &GetSpace(const std::shared_ptr<PureState> &ptr) {
  return ptr->space;
}
inline const QSpace &GetSpace(const std::shared_ptr<ProductState> &ptr) {
  return ptr->space;
}
inline const QSpace &GetSpace(const std::shared_ptr<SumState> &ptr) {
  return ptr->space;
}
inline const QSpace &GetSpace(const ptr_variant &ptr) {
  return std::visit([](const auto &p) -> auto const & { return GetSpace(p); },
                    ptr);
}

// ---- DEBUGGING ---- //
void PrintToStream(std::ostream &os, const std::shared_ptr<PureState> &pure_ptr,
                   size_t indent = 0);
void PrintToStream(std::ostream &os,
                   const std::shared_ptr<ProductState> &product_ptr,
                   size_t indent = 0);
void PrintToStream(std::ostream &os, const std::shared_ptr<SumState> &sum_ptr,
                   size_t indent = 0);
void PrintToStream(std::ostream &os, const ptr_variant &ptr, size_t indent = 0);
void PrintToStream(std::ostream &os, const KetBra &kb);
void PrintToStream(std::ostream &os, const SkewOperator &op);
void Print(const std::shared_ptr<PureState> &pure_ptr, size_t indent = 0);
void Print(const std::shared_ptr<ProductState> &product_ptr, size_t indent = 0);
void Print(const std::shared_ptr<SumState> &sum_ptr, size_t indent = 0);
void Print(const ptr_variant &ptr, size_t indent = 0);
void Print(const KetBra &kb);
void Print(const SkewOperator &op);
std::ostream &operator<<(std::ostream &os, const KetBra &kb);

void Flatten(std::shared_ptr<SumState> &root);
complex Inner_slow(const std::shared_ptr<SumState> &p1,
                   const std::shared_ptr<SumState> &p2);
complex Inner_slow(const ptr_variant &p1, const ptr_variant &p2);

// ---- HELPER ALGORITHMS ---- //
struct Visitor {
  std::function<void(std::shared_ptr<SumState>)> sum_visitor;
  std::function<void(std::shared_ptr<ProductState>)> prod_visitor;
  std::function<void(std::shared_ptr<PureState>)> pure_visitor;
};
// invoke a function recursively to children first, then to self
void InvokeBottomUp(Visitor &visitor, std::shared_ptr<SumState> root);
// invoke a function on self first, then to children
void InvokeTopDown(Visitor &visitor, std::shared_ptr<SumState> root);
void InvokeOnChildren(Visitor &visitor, std::shared_ptr<SumState> state);
void InvokeOnChildren(Visitor &visitor, std::shared_ptr<ProductState> state);
std::vector<std::shared_ptr<SumState>>
CollectSums(std::shared_ptr<SumState> root);
std::vector<std::shared_ptr<ProductState>>
CollectProducts(std::shared_ptr<SumState> root);
std::vector<std::shared_ptr<PureState>>
CollectPures(std::shared_ptr<SumState> root);
std::shared_ptr<SumState> PickRandomSum(std::shared_ptr<SumState> &root,
                                        std::mt19937 &gen);
std::shared_ptr<ProductState> PickRandomProduct(std::shared_ptr<SumState> &root,
                                                std::mt19937 &gen);
std::shared_ptr<PureState> PickRandomPure(std::shared_ptr<SumState> &root,
                                          std::mt19937 &gen);
