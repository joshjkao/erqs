#include "quantumstate.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <nlopt.hpp>
#include <random>
#include <ranges>
#include <spanstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
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

static void skip_whitespace(std::string_view &sv) {
  auto first = sv.find_first_not_of(" \t\n\r");
  if (first != std::string_view::npos) {
    sv.remove_prefix(first);
  }
}
void sync_with_stream(std::string_view &view, std::ispanstream &iss) {
  auto pos = iss.tellg();
  if (pos == std::streampos(-1)) {
    throw std::runtime_error("Stream parsing failed: invalid position.");
  }
  view.remove_prefix(static_cast<std::size_t>(pos));
}

// 1. Pass by reference!
static std::pair<complex, ptr_variant>
parse_helper_var(std::string_view &input);

static std::shared_ptr<SumState> parse_helper_sum(std::string_view &input) {
  skip_whitespace(input);
  if (input.empty() || input.front() != '[') {
    throw std::runtime_error("parse_helper_sum: expected '['");
  }
  input.remove_prefix(1); // Consume '['

  std::vector<complex> coeffs;
  std::vector<ptr_variant> vars;

  // 4. Safely handle whitespace before the closing bracket
  while (true) {
    skip_whitespace(input);
    if (input.empty() || input.front() == ']') {
      break;
    }

    auto [coeff, var] = parse_helper_var(input);
    coeffs.push_back(coeff);
    vars.push_back(var);
  }

  if (input.empty() || input.front() != ']') {
    throw std::runtime_error("parse_helper_sum: expected ']'");
  }
  input.remove_prefix(1); // Consume ']'
  return MakeSum(coeffs, vars);
}

static std::pair<complex, ptr_variant>
parse_helper_var(std::string_view &input) {
  skip_whitespace(input);
  if (input.empty() || input.front() != '(') {
    throw std::runtime_error("parse_helper_var: expected '('");
  }

  std::ispanstream iss{input};
  complex c;
  iss >> c >> std::ws;

  if (iss.peek() == '{') {
    sync_with_stream(input, iss);

    // 3. Consume the '{' so the next parser doesn't choke on it
    input.remove_prefix(1);

    std::vector<std::shared_ptr<SumState>> states;
    while (true) {
      skip_whitespace(input);
      if (input.empty() || input.front() == '}') {
        break;
      }
      states.push_back(parse_helper_sum(input));
    }

    if (input.empty() || input.front() != '}') {
      throw std::runtime_error("parse_helper_var: expected '}'");
    }
    input.remove_prefix(1); // Consume '}'
    return {c, MakeProduct(states)};

  } else {
    QSpace space;
    BitString bits;
    iss >> space >> bits;

    // 2. Sync the stream for leaf nodes too!
    sync_with_stream(input, iss);

    return {c, MakePure(space, bits)};
  }
}
std::shared_ptr<SumState> FromString(std::string_view input) {
  return parse_helper_sum(input);
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
std::string Stringify(const std::shared_ptr<SumState> &state) {
  std::stringstream ss{};
  PrintToStream(ss, state);
  return ss.str();
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

static bool equals_literal_rec(const ptr_variant &var1,
                               const ptr_variant &var2) {
  using pureptr = std::shared_ptr<PureState>;
  using prodptr = std::shared_ptr<ProductState>;
  if (std::holds_alternative<pureptr>(var1) &&
      std::holds_alternative<pureptr>(var2)) {
    auto pure1 = std::get<pureptr>(var1);
    auto pure2 = std::get<pureptr>(var2);
    return *pure1 == *pure2;
  } else if (std::holds_alternative<prodptr>(var1) &&
             std::holds_alternative<prodptr>(var2)) {
    auto prod1 = std::get<prodptr>(var1);
    auto prod2 = std::get<prodptr>(var2);
    if (prod1->states.size() != prod2->states.size())
      return false;
    for (const auto &&[sum1, sum2] :
         std::views::zip(prod1->states, prod2->states)) {
      if (!Equals_literal(sum1, sum2))
        return false;
    }
    return true;
  } else {
    return false;
  }
}
bool Equals_literal(const std::shared_ptr<SumState> &sum1,
                    const std::shared_ptr<SumState> &sum2) {
  if (sum1->states.size() != sum2->states.size())
    return false;
  for (const auto &&[var1, var2] :
       std::views::zip(sum1->states, sum2->states)) {
    if (!equals_literal_rec(var1, var2))
      return false;
  }
  // this is needed to make it easier to compare stringified expressions to
  // string-parsed expressions
  constexpr double epsilon = 1e-5; // Adjust based on your ostream precision
  bool coeffs_match = std::ranges::equal(
      sum1->coeffs, sum2->coeffs, [](const complex &a, const complex &b) {
        return std::abs(a.real() - b.real()) < epsilon &&
               std::abs(a.imag() - b.imag()) < epsilon;
      });
  return coeffs_match;
}

bool Equals_literal_flat(const std::shared_ptr<SumState> &sum1,
                         const std::shared_ptr<SumState> &sum2) {
  std::unordered_map<PureState, complex, PureStateHash> pures1;
  for (const auto &[coeff, var1] :
       std::views::zip(sum1->coeffs, sum1->states)) {
    if (!std::holds_alternative<std::shared_ptr<PureState>>(var1))
      throw std::runtime_error(
          "Equals_literal_flat should only be used on flattened states");
    auto pure1 = std::get<std::shared_ptr<PureState>>(var1);
    pures1[*pure1] += coeff;
  }
  std::unordered_map<PureState, complex, PureStateHash> pures2;
  for (const auto &[coeff, var2] :
       std::views::zip(sum2->coeffs, sum2->states)) {
    if (!std::holds_alternative<std::shared_ptr<PureState>>(var2))
      throw std::runtime_error(
          "Equals_literal_flat should only be used on flattened states");
    auto pure2 = std::get<std::shared_ptr<PureState>>(var2);
    pures2[*pure2] += coeff;
  }
  if (pures1.size() != pures2.size())
    return false;

  constexpr double epsilon = 1e-5;
  return std::ranges::all_of(pures1, [&pures2, epsilon](const auto &pair) {
    const auto &[state, coeff1] = pair;
    auto it = pures2.find(state);
    if (it == pures2.end()) {
      return false; // Key missing in second map
    }
    const auto &coeff2 = it->second;
    return std::abs(coeff1.real() - coeff2.real()) < epsilon &&
           std::abs(coeff1.imag() - coeff2.imag()) < epsilon;
  });
}
