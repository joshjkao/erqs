#include "quantumstate.h"
#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <nlopt.hpp>
#include <random>
#include <ranges>
#include <string>
#include <unordered_map>
#include <variant>

std::shared_ptr<SumState> MakeSum(const std::vector<complex> &coeffs,
                                  const std::vector<ptr_variant> &states) {
  QSpace this_space{0};
  if (coeffs.size() != states.size()) {
    throw std::runtime_error("coeffs and states must be the same size");
  }
  for (auto &state : states) {
    QSpace next_space =
        std::visit([](const auto &e) { return GetSpace(e); }, state);
    if (!this_space.any()) {
      this_space = next_space;
    } else {
      if (this_space != next_space) {
        throw std::runtime_error("sum states must be on same subspace");
      }
    }
  }
  return std::make_shared<SumState>(this_space, coeffs, states);
}

std::shared_ptr<ProductState>
MakeProduct(const std::vector<std::shared_ptr<SumState>> &states) {
  QSpace this_space{0};
  for (auto &state : states) {
    QSpace next_space = state->space;
    if ((this_space & next_space).any()) {
      throw std::runtime_error("product states must be on disjoint subspaces");
    } else {
      this_space |= next_space;
    }
  }
  return std::make_shared<ProductState>(this_space, states);
}

static ptr_variant rclone(const std::shared_ptr<PureState> &pure) {
  return std::make_shared<PureState>(*pure);
}
static std::shared_ptr<SumState> rclone(const std::shared_ptr<SumState> &sum);
static ptr_variant rclone(const std::shared_ptr<ProductState> &prod) {
  std::vector<std::shared_ptr<SumState>> states;
  for (const auto &state : prod->states) {
    states.push_back(rclone(state));
  }
  return std::make_shared<ProductState>(prod->space, states);
}
static std::shared_ptr<SumState> rclone(const std::shared_ptr<SumState> &sum) {
  std::vector<ptr_variant> states;
  for (const auto &state : sum->states) {
    states.push_back(
        std::visit([](const auto &e) { return rclone(e); }, state));
  }
  return std::make_shared<SumState>(sum->space, sum->coeffs, states);
}
std::shared_ptr<SumState> Clone(const std::shared_ptr<SumState> &root) {
  return rclone(root);
}

static std::string get_indentation(size_t indent) {
  return std::string(indent, ' ');
}
void PrintToStream(std::ostream &os, const std::shared_ptr<PureState> &pure_ptr,
                   size_t) {
  os << (*pure_ptr).space << " " << (*pure_ptr).bits;
}
void PrintToStream(std::ostream &os,
                   const std::shared_ptr<ProductState> &product_ptr,
                   size_t indent) {
  os << "\n";
  os << get_indentation(indent);
  os << "{\n";
  for (const auto &sum : (*product_ptr).states) {
    PrintToStream(os, sum, indent + 2);
    os << "\n";
  }
  os << get_indentation(indent);
  os << "}";
}
void PrintToStream(std::ostream &os, const std::shared_ptr<SumState> &sum_ptr,
                   size_t indent) {
  os << get_indentation(indent);
  os << "[\n";
  for (const auto &[c, var] :
       std::views::zip((*sum_ptr).coeffs, (*sum_ptr).states)) {
    os << get_indentation(indent + 2);
    os << c << " ";
    std::visit([&](const auto &s) { PrintToStream(os, s, indent + 2); }, var);
    os << "\n";
  }
  os << get_indentation(indent);
  os << "]";
}
void PrintToStream(std::ostream &os, const ptr_variant &ptr, size_t indent) {
  std::visit([&](const auto &p) { PrintToStream(os, p, indent); }, ptr);
}
void PrintToStream(std::ostream &os, const KetBra &kb) {
  os << kb.coeff << "|";
  std::visit([&](const auto e) { PrintToStream(os, e, 2); }, kb.ket);
  os << "><";
  std::visit([&](const auto e) { PrintToStream(os, e, 2); }, kb.bra);
  os << "|";
}
void PrintToStream(std::ostream &os, const SkewOperator &op) {
  for (const auto &kb : op.ketbras) {
    PrintToStream(os, kb);
    os << "\n";
  }
}

void Print(const std::shared_ptr<PureState> &pure_ptr, size_t) {
  PrintToStream(std::cout, pure_ptr);
}
void Print(const std::shared_ptr<ProductState> &product_ptr, size_t indent) {
  PrintToStream(std::cout, product_ptr, indent);
}
void Print(const std::shared_ptr<SumState> &sum_ptr, size_t indent) {
  PrintToStream(std::cout, sum_ptr, indent);
}
void Print(const ptr_variant &ptr, size_t indent) {
  PrintToStream(std::cout, ptr, indent);
}
void Print(const KetBra &kb) { PrintToStream(std::cout, kb); }
void Print(const SkewOperator &op) { PrintToStream(std::cout, op); }
std::ostream &operator<<(std::ostream &os, const KetBra &kb) {
  PrintToStream(os, kb);
  return os;
}

// needed to perform compression
struct PureStateHash {
  std::size_t operator()(const PureState &t) const {
    std::size_t seed = 0;
    auto hash_combine = [&seed](std::size_t hash_value) {
      seed ^= hash_value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    std::hash<BitString> bitset_hasher;
    hash_combine(bitset_hasher(t.space));
    hash_combine(bitset_hasher(t.bits));
    return seed;
  }
};
static std::vector<std::pair<complex, PureState>>
flatten_helper(std::shared_ptr<PureState> &pure_ptr) {
  return {{1.0, *pure_ptr}};
}
static std::vector<std::pair<complex, PureState>>
flatten_helper(std::shared_ptr<ProductState> &prod_ptr);
static std::vector<std::pair<complex, PureState>>
flatten_helper(std::shared_ptr<SumState> &sum_ptr) {
  std::unordered_map<PureState, complex, PureStateHash> map;
  std::vector<std::pair<complex, PureState>> ret;
  for (auto &&[c1, s1] :
       std::views::zip((*sum_ptr).coeffs, (*sum_ptr).states)) {
    auto term_flat = std::visit([](auto &e) { return flatten_helper(e); }, s1);
    for (auto &[c2, s2] : term_flat) {
      ret.push_back({c1 * c2, s2});
    }
  }
  return ret;
}
static std::vector<std::pair<complex, PureState>>
flatten_prod_helper(std::shared_ptr<ProductState> &prod_ptr) {
  if (prod_ptr->states.size() == 1) {
    auto ret = flatten_helper(prod_ptr->states.back());
    prod_ptr->states.pop_back();
    return ret;
  } else {
    std::vector<std::pair<complex, PureState>> ret;
    auto this_flat = flatten_helper(prod_ptr->states.back());
    prod_ptr->states.pop_back();
    auto next_flat = flatten_prod_helper(prod_ptr);
    for (auto &[c1, s1] : this_flat) {
      for (auto &[c2, s2] : next_flat) {
        ret.push_back({c1 * c2, {s1.space | s2.space, s1.bits | s2.bits}});
      }
    }
    return ret;
  }
}
static std::vector<std::pair<complex, PureState>>
flatten_helper(std::shared_ptr<ProductState> &prod_ptr) {
  return flatten_prod_helper(prod_ptr);
}
void Flatten(std::shared_ptr<SumState> &root) {
  auto flat = flatten_helper(root);
  std::unordered_map<PureState, complex, PureStateHash> map;
  for (auto &[coeff, state] : flat) {
    map[state] += coeff;
  }
  root->coeffs.clear();
  root->states.clear();
  for (auto &[state, coeff] : map) {
    root->coeffs.push_back(coeff);
    root->states.push_back(std::make_shared<PureState>(state));
  }
  // will error if empty
  assert(!flat.empty());
  root->space = flat[0].second.space;
}

complex Inner_slow(const ptr_variant &p1, const ptr_variant &p2) {
  auto clone1 = std::visit([](auto &p) { return rclone(p); }, p1);
  auto clone2 = std::visit([](auto &p) { return rclone(p); }, p2);
  auto flat1 = std::visit([](auto &p) { return flatten_helper(p); }, clone1);
  auto flat2 = std::visit([](auto &p) { return flatten_helper(p); }, clone2);
  complex ret = 0.0;
  for (const auto &[c1, s1] : flat1) {
    for (const auto &[c2, s2] : flat2) {
      if (s1.bits == s2.bits && s1.space == s2.space)
        ret += std::conj(c1) * c2;
    }
  }
  return ret;
}
complex Inner_slow(const std::shared_ptr<SumState> &p1,
                   const std::shared_ptr<SumState> &p2) {
  auto clone1 = Clone(p1);
  auto clone2 = Clone(p2);
  auto flat1 = flatten_helper(clone1);
  std::unordered_map<PureState, complex, PureStateHash> map1;
  for (const auto &[c, s] : flat1) {
    map1[s] += c;
  }
  auto flat2 = flatten_helper(clone2);
  std::unordered_map<PureState, complex, PureStateHash> map2;
  for (const auto &[c, s] : flat2) {
    map2[s] += c;
  }
  complex ret = 0.0;
  for (const auto &[s, c] : map2) {
    auto it = map1.find(s);
    if (it != map1.end()) {
      ret += std::conj(it->second) * c;
    }
  }
  return ret;
}

static std::bitset<NQUBITS> random_bitset(std::mt19937 &gen) {
  std::bernoulli_distribution dist{0.5};
  QSpace bits;
  for (size_t i = 0; i < NQUBITS; ++i) {
    bits[i] = dist(gen);
  }
  return bits;
}

static std::shared_ptr<SumState> random_sum(size_t max_depth, size_t n_terms,
                                            size_t n_factors, QSpace space,
                                            std::mt19937 &gen);
static ptr_variant random_product(size_t max_depth, size_t n_terms,
                                  size_t n_factors, QSpace space,
                                  std::mt19937 &gen) {
  if (max_depth <= 1) {
    QSpace bits = random_bitset(gen);
    bits &= space;
    return MakePure(space, bits);
  } else {
    std::vector<std::shared_ptr<SumState>> states;
    for (size_t i = 0; i < n_terms && space.any(); ++i) {
      QSpace new_space = random_bitset(gen);
      new_space &= space;
      if (new_space.none())
        continue;
      space ^= new_space;
      states.push_back(
          random_sum(max_depth - 1, n_terms, n_factors, new_space, gen));
    }
    if (space.any()) {
      states.push_back(
          random_sum(max_depth - 1, n_terms, n_factors, space, gen));
    }
    return MakeProduct(states);
  }
}
static std::shared_ptr<SumState> random_sum(size_t max_depth, size_t n_terms,
                                            size_t n_factors, QSpace space,
                                            std::mt19937 &gen) {
  std::uniform_real_distribution<double> dist{-1.0, 1.0};
  std::vector<complex> coeffs(n_terms);
  for (auto &e : coeffs) {
    e = complex{dist(gen), dist(gen)};
  }
  std::vector<ptr_variant> states{n_terms};
  for (auto &s : states) {
    s = random_product(max_depth - 1, n_terms, n_factors, space, gen);
  }
  return MakeSum(coeffs, states);
}
std::shared_ptr<SumState> RandomSumState(size_t max_depth, size_t n_terms,
                                         size_t n_factors, QSpace space,
                                         std::mt19937 &gen) {
  return random_sum(max_depth, n_terms, n_factors, space, gen);
}
ptr_variant RandomProductState(size_t max_depth, size_t n_terms,
                               size_t n_factors, QSpace space,
                               std::mt19937 &gen) {
  return random_product(max_depth, n_terms, n_factors, space, gen);
}

// EXPERIMENTAL: arbitrary visitors
static void recursive_invoke(Visitor &visitor,
                             std::shared_ptr<ProductState> prod);
static void recursive_invoke(Visitor &visitor, std::shared_ptr<PureState> pure);
static void recursive_invoke(Visitor &visitor, std::shared_ptr<SumState> sum) {
  for (auto var : sum->states) {
    std::visit([&](auto p) { recursive_invoke(visitor, p); }, var);
  }
  visitor.sum_visitor(sum);
}
static void recursive_invoke(Visitor &visitor,
                             std::shared_ptr<ProductState> prod) {
  for (auto sum : prod->states) {
    recursive_invoke(visitor, sum);
  }
  visitor.prod_visitor(prod);
}
static void recursive_invoke(Visitor &visitor,
                             std::shared_ptr<PureState> pure) {
  visitor.pure_visitor(pure);
}
void InvokeBottomUp(Visitor &visitor, std::shared_ptr<SumState> root) {
  recursive_invoke(visitor, root);
}
static void invoke(Visitor &visitor, std::shared_ptr<SumState> state) {
  visitor.sum_visitor(state);
}
static void invoke(Visitor &visitor, std::shared_ptr<ProductState> state) {
  visitor.prod_visitor(state);
}
static void invoke(Visitor &visitor, std::shared_ptr<PureState> state) {
  visitor.pure_visitor(state);
}
void InvokeOnChildren(Visitor &visitor, std::shared_ptr<SumState> state) {
  for (auto child : state->states) {
    std::visit([&](auto v) { invoke(visitor, v); }, child);
  }
}
void InvokeOnChildren(Visitor &visitor, std::shared_ptr<ProductState> state) {
  for (auto child : state->states) {
    invoke(visitor, child);
  }
}

std::vector<std::shared_ptr<SumState>>
CollectSums(std::shared_ptr<SumState> root) {
  std::vector<std::shared_ptr<SumState>> ret;
  Visitor visitor{[&](std::shared_ptr<SumState> p) { ret.push_back(p); },
                  [](std::shared_ptr<ProductState>) {},
                  [](std::shared_ptr<PureState>) {}};
  InvokeBottomUp(visitor, root);
  return ret;
}
std::vector<std::shared_ptr<ProductState>>
CollectProducts(std::shared_ptr<SumState> root) {
  std::vector<std::shared_ptr<ProductState>> ret;
  Visitor visitor{[](std::shared_ptr<SumState>) {},
                  [&](std::shared_ptr<ProductState> p) { ret.push_back(p); },
                  [](std::shared_ptr<PureState>) {}};
  InvokeBottomUp(visitor, root);
  return ret;
}
std::vector<std::shared_ptr<PureState>>
CollectPures(std::shared_ptr<SumState> root) {
  std::vector<std::shared_ptr<PureState>> ret;
  Visitor visitor{[](std::shared_ptr<SumState>) {},
                  [](std::shared_ptr<ProductState>) {},
                  [&](std::shared_ptr<PureState> p) { ret.push_back(p); }};
  InvokeBottomUp(visitor, root);
  return ret;
}

std::shared_ptr<SumState> PickRandomSum(std::shared_ptr<SumState> &root,
                                        std::mt19937 &gen) {
  auto sums = CollectSums(root);
  std::uniform_int_distribution<size_t> dist(0, sums.size() - 1);
  size_t rand_ind = dist(gen);
  return sums[rand_ind];
}
std::shared_ptr<ProductState> PickRandomProduct(std::shared_ptr<SumState> &root,
                                                std::mt19937 &gen) {
  auto prods = CollectProducts(root);
  if (prods.empty())
    return nullptr;
  std::uniform_int_distribution<size_t> dist(0, prods.size() - 1);
  size_t rand_ind = dist(gen);
  return prods[rand_ind];
}
std::shared_ptr<PureState> PickRandomPure(std::shared_ptr<SumState> &root,
                                          std::mt19937 &gen) {

  auto pures = CollectPures(root);
  std::uniform_int_distribution<size_t> dist(0, pures.size() - 1);
  size_t rand_ind = dist(gen);
  return pures[rand_ind];
}

bool Equals(const std::shared_ptr<PureState> &p1,
            const std::shared_ptr<PureState> &p2) {
  return p1->bits == p2->bits && p1->space == p2->space;
}

bool Equals_slow(const std::shared_ptr<SumState> &sum1,
                 const std::shared_ptr<SumState> &sum2) {
  auto self_inner_1 = Inner_slow(sum1, sum1);
  auto inner_12 = Inner_slow(sum1, sum2);
  return std::abs(self_inner_1.real() - inner_12.real()) <
         1e-4 * self_inner_1.real();
}
