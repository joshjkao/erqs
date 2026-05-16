#include "optimization.h"
#include "operations.h"
#include "quantumstate.h"
// #include <iostream>
#include <cstddef>
#include <memory>
#include <nlopt.hpp>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <unordered_map>

static void set_coefficients(std::vector<complex> &,
                             std::shared_ptr<PureState> &) {
  return;
}
static void set_coefficients(std::vector<complex> &coeffs,
                             std::shared_ptr<SumState> &state);
static void set_coefficients(std::vector<complex> &coeffs,
                             std::shared_ptr<ProductState> &state) {
  for (auto &s : state->states) {
    set_coefficients(coeffs, s);
  }
}
static void set_coefficients(std::vector<complex> &coeffs,
                             std::shared_ptr<SumState> &state) {
  for (auto &c : state->coeffs) {
    c = coeffs.back();
    coeffs.pop_back();
  }
  for (auto &s : state->states) {
    std::visit([&](auto &e) { set_coefficients(coeffs, e); }, s);
  }
}
void SetCoefficients(const std::vector<complex> &coeffs,
                     std::shared_ptr<SumState> &state) {
  auto coeffscp = coeffs;
  set_coefficients(coeffscp, state);
}
void SetCoefficients(const std::vector<double> &coeffs,
                     std::shared_ptr<SumState> &state) {
  std::vector<complex> coeffs_complex;
  for (size_t i = 0; i < coeffs.size(); ++i) {
    double re = coeffs[i++];
    double im = coeffs[i];
    coeffs_complex.push_back({re, im});
  }
  SetCoefficients(coeffs_complex, state);
}

static size_t count_coefficients(const std::shared_ptr<PureState> &) {
  return 0;
}
static size_t count_coefficients(const std::shared_ptr<ProductState> &state) {
  size_t ret = 0;
  for (const auto &s : state->states) {
    ret += CountCoefficients(s);
  }
  return ret;
}
size_t CountCoefficients(const std::shared_ptr<SumState> &state) {
  size_t ret = state->coeffs.size();
  for (const auto &s : state->states) {
    ret += std::visit([](const auto &e) { return count_coefficients(e); }, s);
  }
  return ret;
}

static size_t count_nodes(const std::shared_ptr<PureState> &) { return 1; }
static size_t count_nodes(const std::shared_ptr<ProductState> &state) {
  size_t ret = 0;
  for (const auto &s : state->states) {
    ret += CountNodes(s);
  }
  return ret;
}
size_t CountNodes(const std::shared_ptr<SumState> &state) {
  size_t ret = 0;
  for (const auto &s : state->states) {
    ret += std::visit([](const auto &e) { return count_nodes(e); }, s);
  }
  return ret;
}

using Callback = std::function<double(const std::vector<double> &)>;
static double obj(const std::vector<double> &x, std::vector<double> &,
                  void *data) {
  Callback *f = static_cast<Callback *>(data);
  double ret = (*f)(x);
  return ret;
}
double OptimizeCoefficients(const PauliHamiltonian &H,
                            std::shared_ptr<SumState> &root) {
  size_t num_coeffs = CountCoefficients(root);
  Callback e = [&](const std::vector<double> &x) -> double {
    SetCoefficients(x, root);
    Normalize(root);
    double ret = ExpectedValue(H, root);
    return ret;
  };
  std::vector<double> grad;
  std::vector<complex> coeffs(num_coeffs);
  SetCoefficients(coeffs, root);
  nlopt::opt opt(nlopt::LN_COBYLA, static_cast<unsigned int>(num_coeffs) * 2);
  opt.set_min_objective(obj, static_cast<void *>(&e));
  opt.set_xtol_rel(1e-4);
  std::vector<double> x(num_coeffs * 2, 1.0);
  double minf;
  nlopt::result result = opt.optimize(x, minf);
  if (result > 0) {
    SetCoefficients(x, root);
    Normalize(root);
    return minf;
  } else {
    throw std::runtime_error("error during minimuzation step");
  }
}

static void prune(std::shared_ptr<PureState> &, double) { return; }
static void prune(std::shared_ptr<ProductState> &state, double threshold) {
  std::vector<std::shared_ptr<SumState>> new_states;
  for (auto &s : state->states) {
    Prune(s, threshold);
    if (s->states.empty())
      continue;
    new_states.push_back(s);
  }
  state->states = new_states;
}
static bool prune_should_keep(std::shared_ptr<PureState>) { return true; }
static bool prune_should_keep(std::shared_ptr<ProductState> p) {
  return !p->states.empty();
}
static bool prune_should_keep(ptr_variant v) {
  return std::visit([](auto p) { return prune_should_keep(p); }, v);
}
void Prune(std::shared_ptr<SumState> &root, double threshold) {
  for (auto &p : root->states) {
    std::visit([=](auto &s) { prune(s, threshold); }, p);
  }
  std::vector<complex> coeffs_new;
  std::vector<ptr_variant> states_new;
  for (size_t i = 0; i < root->coeffs.size(); ++i) {
    double re = root->coeffs[i].real();
    double im = root->coeffs[i].imag();
    double mag = (re * re + im * std::conj(im)).real();
    if (mag >= threshold && prune_should_keep(root->states[i])) {
      coeffs_new.push_back(root->coeffs[i]);
      states_new.push_back(root->states[i]);
    }
  }
  root->coeffs = coeffs_new;
  root->states = states_new;
}

void AddRandomTerm(std::shared_ptr<SumState> &root, std::mt19937 &gen) {
  std::shared_ptr<SumState> state = PickRandomSum(root, gen);
  ptr_variant rand_state = RandomProductState(2, 3, 3, state->space, gen);
  state->coeffs.push_back({1.0, 0.0});
  state->states.push_back(rand_state);
}

void UnifyRandom(std::shared_ptr<SumState> &root, std::mt19937 &gen) {
  std::shared_ptr<ProductState> state = PickRandomProduct(root, gen);
  if (state == nullptr)
    return;
  if (state->states.size() < 2)
    return;
  std::uniform_int_distribution<size_t> dist1(0uz, state->states.size() - 1uz);
  size_t i1 = dist1(gen);
  std::shared_ptr<SumState> s1 = state->states[i1];
  state->states.erase(state->states.begin() + static_cast<std::ptrdiff_t>(i1));
  std::uniform_int_distribution<size_t> dist2(0uz, state->states.size() - 1uz);
  size_t i2 = dist2(gen);
  std::shared_ptr<SumState> s2 = state->states[i2];
  state->states.erase(state->states.begin() + static_cast<std::ptrdiff_t>(i2));
  std::shared_ptr<SumState> new_state = Tensor(s1, s2);
  state->states.push_back(new_state);
}

// Really don't like this solution
std::vector<std::shared_ptr<SumState>>
get_states_if_product(const std::shared_ptr<PureState> &) {
  return {};
}
std::vector<std::shared_ptr<SumState>>
get_states_if_product(const std::shared_ptr<ProductState> &p) {
  return p->states;
}
std::vector<std::shared_ptr<SumState>>
get_states_if_product(const ptr_variant &var) {
  return std::visit([](const auto &v) { return get_states_if_product(v); },
                    var);
}
static void remove_singles_rec(std::shared_ptr<PureState> &);
static void remove_singles_rec(std::shared_ptr<ProductState> &state);
static void remove_singles_rec(std::shared_ptr<SumState> &state);
static void remove_singles_rec(ptr_variant &var) {
  std::visit([](auto &p) { remove_singles_rec(p); }, var);
}
static void remove_singles_rec(std::shared_ptr<PureState> &) { return; }
static void remove_singles_rec(std::shared_ptr<SumState> &state) {
  for (auto var : state->states) {
    remove_singles_rec(var);
  }
  std::vector<ptr_variant> states_new;
  // dropping a coefficient, technically
  for (auto var : state->states) {
    auto var_states = get_states_if_product(var);
    // the variant is a pure state
    if (var_states.empty())
      states_new.push_back(var);
    // the variant has more than one factor state
    else if (var_states.size() > 1)
      states_new.push_back(var);
    // the variant is a product state with a single term
    else {
      for (auto child_var : var_states[0]->states) {
        states_new.push_back(child_var);
      }
    }
  }
  state->states = states_new;
}
static void remove_singles_rec(std::shared_ptr<ProductState> &state) {
  // handle recursive call first
  for (auto &sum : state->states) {
    remove_singles_rec(sum);
  }
  // we're going to replace this state's vector
  std::vector<std::shared_ptr<SumState>> states_new;
  for (auto &sum : state->states) {
    if (sum->states.size() == 1) {
      auto sum_states = get_states_if_product(sum->states[0]);
      if (sum_states.size() == 0) {
        // sum state containing a single pure state
        // we can't remove this because pure states must
        // be contained in
        // a sum state
        states_new.push_back(sum);
      } else {
        // sum state containing a single product state
        complex coeff = sum->coeffs[0];
        // the coeff should get multiplied only onto one
        // of the sum
        // states (the first one)
        for (auto &sum1 : sum_states | std::views::take(1)) {
          for (auto &c : sum1->coeffs) {
            c *= coeff;
          }
          states_new.push_back(sum1);
        }
        for (auto &sum1 : sum_states | std::views::drop(1)) {
          states_new.push_back(sum1);
        }
      }
    } else {
      // this state should be copied over without removing
      states_new.push_back(sum);
    }
  }
  // the old vector's destructor will handle deallocating states that
  // should be removed
  state->states = states_new;
}

void RemoveSingles(std::shared_ptr<SumState> &root) {
  for (auto &var : root->states) {
    remove_singles_rec(var);
  }
}

static complex normalize_rec(std::shared_ptr<PureState>) { return 1.0; }
static complex normalize_rec(std::shared_ptr<ProductState> state);
static complex normalize_rec(ptr_variant var) {
  return std::visit([](auto p) { return normalize_rec(p); }, var);
}
static complex normalize_rec(std::shared_ptr<SumState> state) {
  for (auto &&[c, v] : std::views::zip(state->coeffs, state->states))
    c *= normalize_rec(v);
  complex norm2 = 0.0;
  for (auto c : state->coeffs)
    norm2 += c * std::conj(c);
  complex norm = sqrt(norm2.real());
  for (auto &c : state->coeffs)
    c /= norm;
  return norm;
}
static complex normalize_rec(std::shared_ptr<ProductState> state) {
  complex prod = 1.0;
  for (auto s : state->states) {
    prod *= normalize_rec(s);
  }
  return prod;
}
void Normalize(std::shared_ptr<SumState> &root) { normalize_rec(root); }

// needed to perform compression
using PureTuple = std::tuple<QSpace, QSpace>;
using CoeffPureTuple = std::tuple<complex, QSpace, QSpace>;
struct PureTupleHash {
  std::size_t operator()(const PureTuple &t) const {
    std::size_t seed = 0;
    auto hash_combine = [&seed](std::size_t hash_value) {
      seed ^= hash_value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    std::hash<BitString> bitset_hasher;
    hash_combine(bitset_hasher(std::get<0>(t)));
    hash_combine(bitset_hasher(std::get<1>(t)));
    return seed;
  }
};
static std::optional<CoeffPureTuple>
get_tuple_if_pure(std::shared_ptr<SumState> &ptr) {
  // assume we're already as simple as possible
  if (ptr->states.size() != 1)
    return std::nullopt;
  auto pure_ptr_ptr = std::get_if<std::shared_ptr<PureState>>(&ptr->states[0]);
  if (!pure_ptr_ptr)
    return std::nullopt;
  return CoeffPureTuple{ptr->coeffs[0], (*pure_ptr_ptr)->space,
                        (*pure_ptr_ptr)->bits};
}
static std::optional<CoeffPureTuple>
get_tuple_if_pure(std::shared_ptr<PureState> &ptr) {
  return CoeffPureTuple{complex{1.0}, ptr->space, ptr->bits};
}
static std::optional<CoeffPureTuple>
get_tuple_if_pure(std::shared_ptr<ProductState> &ptr) {
  if (ptr->states.size() != 1)
    return std::nullopt;
  return get_tuple_if_pure(ptr->states[0]);
}
std::shared_ptr<SumState> Simplify(std::shared_ptr<SumState> &ptr);
ptr_variant Simplify(std::shared_ptr<PureState> &ptr) { return ptr; }
ptr_variant Simplify(std::shared_ptr<ProductState> &ptr) {
  complex coeff{1.0};
  QSpace space{0};
  QSpace bits{0};
  std::vector<std::shared_ptr<SumState>> sums;
  for (auto &sum : ptr->states) {
    sum = Simplify(sum);
    auto tuple_if_pure = get_tuple_if_pure(sum);
    if (tuple_if_pure) {
      auto &[tup_coeff, tup_space, tup_bits] = tuple_if_pure.value();
      coeff *= coeff;
      space |= tup_space;
      bits |= tup_bits;
    } else {
      sums.push_back(sum);
    }
  }

  auto temp_pure = MakePure(space, bits);
  if (sums.empty()) {
    // sums empty means we should return the pure state directly instead of
    // wrapping it in a product
    return temp_pure;
  } else {
    // otherwise, we have nontrivial sum states
    if (space.any() || coeff != complex{1.0}) {
      // if our pure state is just a representation of the number 1, keeping it
      // is redundant
      auto temp_sum = MakeSum({coeff}, {temp_pure});
      sums.push_back(temp_sum);
    }
    return MakeProduct(sums);
  }
}
ptr_variant Simplify(ptr_variant &ptr) {
  return std::visit([](auto &v) { return Simplify(v); }, ptr);
}
std::shared_ptr<SumState> Simplify(std::shared_ptr<SumState> &ptr) {
  std::unordered_map<PureTuple, complex, PureTupleHash> pures;
  std::vector<complex> coeffs;
  std::vector<ptr_variant> vars;
  for (auto &&[coeff, var] : std::views::zip(ptr->coeffs, ptr->states)) {
    var = Simplify(var);
    auto tuple_if_pure =
        std::visit([](auto &v) { return get_tuple_if_pure(v); }, var);
    if (tuple_if_pure) {
      auto &[tup_coeff, tup_space, tup_bits] = tuple_if_pure.value();
      PureTuple tuple{tup_space, tup_bits};
      pures[tuple] += coeff * tup_coeff;
    } else {
      coeffs.push_back(coeff);
      vars.push_back(var);
    }
  }
  for (const auto &[tuple, coeff] : pures) {
    coeffs.push_back(coeff);
    auto &[space, bits] = tuple;
    vars.push_back(MakePure(space, bits));
  }
  return MakeSum(coeffs, vars);
}

static std::optional<QSpace> bits_if_pure(std::shared_ptr<PureState> p) {
  return p->bits;
}
static std::optional<QSpace> bits_if_pure(std::shared_ptr<ProductState> p);
static std::optional<QSpace> bits_if_pure(ptr_variant var) {
  return std::visit([](auto v) { return bits_if_pure(v); }, var);
}
static std::optional<QSpace> bits_if_pure(std::shared_ptr<ProductState> p) {
  if (p->states.size() != 1)
    return std::nullopt;
  else if (p->states[0]->states.size() != 1)
    return std::nullopt;
  else
    return bits_if_pure(p->states[0]->states[0]);
}
void SquashPures(std::shared_ptr<SumState> &root) {
  auto sums = CollectSums(root);
  for (auto sum : sums) {
    std::unordered_map<QSpace, complex> pures;
    std::vector<ptr_variant> states_new;
    std::vector<complex> coeffs_new;
    for (auto &&[c, var] : std::views::zip(sum->coeffs, sum->states)) {
      auto bits = std::visit([](auto v) { return bits_if_pure(v); }, var);
      if (bits) {
        if (pures.contains(bits.value()))
          pures[bits.value()] += c;
        else
          pures[bits.value()] = c;
      } else {
        coeffs_new.push_back(c);
        states_new.push_back(var);
      }
    }
    for (auto [bits, c] : pures) {
      coeffs_new.push_back(c);
      states_new.push_back(MakePure(sum->space, bits));
    }
    sum->coeffs = coeffs_new;
    sum->states = states_new;
  }
  // auto prods = CollectProducts(root);
  // for (auto prod : prods) {
  //   std::vector<std::shared_ptr<SumState>> states_new;
  //   QSpace space_new{0};
  //   QSpace bits_new{0};
  //   complex coeff_new = 1.0;
  //   for (auto s : prod->states) {
  //     if (s->states.size() != 1) {
  //       states_new.push_back(s);
  //       continue;
  //     }
  //     auto bits = bits_if_pure(s->states[0]);
  //     if (bits) {
  //       bits_new |= bits.value();
  //       space_new |= GetSpace(s->states[0]);
  //       coeff_new *= s->coeffs[0];
  //     } else {
  //       states_new.push_back(s);
  //     }
  //   }
  //   if (space_new.any()) {
  //     auto pure_new = MakePure(space_new, bits_new);
  //     auto sum_new = MakeSum({coeff_new}, {pure_new});
  //     states_new.push_back(sum_new);
  //   }
  //   prod->states = states_new;
  // }
}

// void Factorize(std::shared_ptr<SumState> &root) {
// 	// unclear how to do this, it's not very well described in the notes
// 	// "wherever separability emerges" is what he says, maybe i'll have to
// ask
// 	// him
// 	// maybe look at sum states with more than one pure state, then look for
// the
// 	// bits they share
// 	// for example, say we have a 11 and a 10: this could be rewritten as 1
// 	// tensor 0+1
// 	// it goes from two coefficients to two coefficients, though instead of
// two
// 	// pure states in a sum state it becomes a sumstate with a product state
// 	// between two sums
// }
