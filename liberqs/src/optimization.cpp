#include "optimization.hpp"
#include "operations.hpp"
#include "quantumstate.hpp"
#include <cmath>
#include <cstddef>
#include <memory>
#include <nlopt.hpp>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <variant>

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
    if (coeffs.empty())
      throw std::runtime_error("not enough coefficients to set");
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
  if (!coeffscp.empty())
    throw std::runtime_error("too many coefficients passed to set");
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

void SetCoefficients_shallow(const std::vector<complex> &coeffs,
                             std::shared_ptr<SumState> &state) {
  if (state->coeffs.size() != coeffs.size())
    throw std::runtime_error("list size mismatch");
  state->coeffs = coeffs;
}
void SetCoefficients_shallow(const std::vector<double> &coeffs,
                             std::shared_ptr<SumState> &state) {
  std::vector<complex> coeffs_complex;
  for (size_t i = 0; i < coeffs.size(); ++i) {
    double re = coeffs[i++];
    double im = coeffs[i];
    coeffs_complex.push_back({re, im});
  }
  SetCoefficients_shallow(coeffs_complex, state);
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
};

size_t CountCoefficients_shallow(const std::shared_ptr<SumState> &state) {
  return state->coeffs.size();
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
    double ret = ExpectedValue(H, root);
    return ret;
  };
  std::vector<double> grad;
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

double OptimizeCoefficients_local(const PauliHamiltonian &H,
                                  std::shared_ptr<SumState> &root,
                                  std::shared_ptr<SumState> &state) {
  size_t num_coeffs = state->coeffs.size();
  Callback e = [&](const std::vector<double> &x) -> double {
    SetCoefficients_shallow(x, state);
    double ret = ExpectedValue(H, root);
    return ret;
  };
  std::vector<double> grad;
  nlopt::opt opt(nlopt::LN_COBYLA, static_cast<unsigned int>(num_coeffs) * 2);
  opt.set_min_objective(obj, static_cast<void *>(&e));
  opt.set_xtol_rel(1e-4);
  std::vector<double> x(num_coeffs * 2, 1.0);
  double minf;
  nlopt::result result = opt.optimize(x, minf);
  if (result > 0) {
    SetCoefficients_shallow(x, state);
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

    if (std::sqrt(mag) >= threshold && prune_should_keep(root->states[i])) {
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
  std::vector<complex> coeffs_new;
  for (auto &&[coeff, var] : std::views::zip(state->coeffs, state->states)) {
    auto var_states = get_states_if_product(var);
    // the variant is a pure state
    if (var_states.empty()) {
      coeffs_new.push_back(coeff);
      states_new.push_back(var);
    }
    // the variant has more than one factor state
    else if (var_states.size() > 1) {
      coeffs_new.push_back(coeff);
      states_new.push_back(var);
    }
    // the variant is a product state with a single term
    else {
      for (auto &&[child_coeff, child_var] :
           std::views::zip(var_states[0]->coeffs, var_states[0]->states)) {
        coeffs_new.push_back(coeff * child_coeff);
        states_new.push_back(child_var);
      }
    }
  }
  state->coeffs = coeffs_new;
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

void Normalize_slow(std::shared_ptr<SumState> &root) {
  auto self_inner = Inner_slow(root, root);
  complex norm = std::sqrt(self_inner.real());
  for (auto &coeff : root->coeffs) {
    coeff /= norm;
  }
}

void Normalize(std::shared_ptr<SumState> &root) {
  auto self_inner = Inner(root, root);
  complex norm = std::sqrt(self_inner[0].coeff.real());
  for (auto &coeff : root->coeffs) {
    coeff /= norm;
  }
}

using CoeffPureTuple = std::tuple<complex, QSpace, QSpace>;
template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

// Forward declare to allow mutual recursion
std::optional<CoeffPureTuple> ExtractPure(const ptr_variant &var);
std::optional<CoeffPureTuple> ExtractPure(const std::shared_ptr<SumState> &s);

std::optional<CoeffPureTuple> ExtractPure(const std::shared_ptr<SumState> &s) {
  if (s->states.size() == 1) {
    if (auto *p = std::get_if<std::shared_ptr<PureState>>(&s->states[0])) {
      return CoeffPureTuple{s->coeffs[0], (*p)->space, (*p)->bits};
    }
  }
  return std::nullopt;
}

std::optional<CoeffPureTuple> ExtractPure(const ptr_variant &var) {
  return std::visit(overloaded{[](const std::shared_ptr<PureState> &p)
                                   -> std::optional<CoeffPureTuple> {
                                 return CoeffPureTuple{complex{1.0}, p->space,
                                                       p->bits};
                               },
                               [](const std::shared_ptr<ProductState> &prod)
                                   -> std::optional<CoeffPureTuple> {
                                 if (prod->states.size() == 1) {
                                   return ExtractPure(prod->states[0]);
                                 }
                                 return std::nullopt;
                               }},
                    var);
}

static std::shared_ptr<SumState> simplify(std::shared_ptr<SumState> &ptr);
static ptr_variant simplify(std::shared_ptr<PureState> &ptr) { return ptr; }

static ptr_variant simplify(std::shared_ptr<ProductState> &ptr) {
  complex coeff{1.0};
  QSpace space{0}, bits{0};
  std::vector<std::shared_ptr<SumState>> sums;

  for (auto &sum : ptr->states) {
    sum = simplify(sum); // Simplify child sums

    // Calls ExtractPure(const std::shared_ptr<SumState>&)
    if (auto pure_opt = ExtractPure(sum)) {
      auto &[tup_coeff, tup_space, tup_bits] = *pure_opt;
      coeff *= tup_coeff;
      space |= tup_space;
      bits |= tup_bits;
    } else {
      sums.push_back(sum);
    }
  }

  auto temp_pure = MakePure(space, bits);
  if (sums.empty()) {
    if (coeff == complex{1.0})
      return temp_pure;
    return MakeProduct({MakeSum({coeff}, {temp_pure})});
  } else {
    if (space.any() || coeff != complex{1.0}) {
      sums.push_back(MakeSum({coeff}, {temp_pure}));
    }
    return MakeProduct(sums);
  }
}

static ptr_variant simplify(ptr_variant &ptr) {
  return std::visit([](auto &v) { return simplify(v); }, ptr);
}

static std::shared_ptr<SumState> simplify(std::shared_ptr<SumState> &ptr) {
  std::unordered_map<PureState, complex, PureStateHash> pures;
  std::vector<complex> coeffs;
  std::vector<ptr_variant> vars;

  auto add_term = [&](complex c, const ptr_variant &v) {
    // Calls ExtractPure(const ptr_variant&)
    if (auto pure_opt = ExtractPure(v)) {
      auto &[tup_coeff, tup_space, tup_bits] = *pure_opt;
      pures[PureState{tup_space, tup_bits}] += c * tup_coeff;
    } else {
      coeffs.push_back(c);
      vars.push_back(v);
    }
  };

  for (auto &&[coeff, var] : std::views::zip(ptr->coeffs, ptr->states)) {
    var = simplify(var); // Simplify the variant

    // Strictly match what is actually inside ptr_variant
    std::visit(overloaded{[&](const std::shared_ptr<ProductState> &prod) {
                            if (prod->states.size() == 1) {
                              // prod->states[0] is a SumState, so zip its
                              // coeffs and variants
                              for (auto &&[child_coeff, child_var] :
                                   std::views::zip(prod->states[0]->coeffs,
                                                   prod->states[0]->states)) {
                                add_term(coeff * child_coeff, child_var);
                              }
                            } else {
                              add_term(coeff, var);
                            }
                          },
                          [&](const std::shared_ptr<PureState> &) {
                            add_term(coeff, var);
                          }},
               var);
  }

  for (const auto &[tuple, c] : pures) {
    coeffs.push_back(c);
    vars.push_back(MakePure(tuple.space, tuple.bits));
  }

  return MakeSum(coeffs, vars);
}

void Simplify(ptr_variant &var) { var = simplify(var); }

void Simplify(std::shared_ptr<SumState> &ptr) { ptr = simplify(ptr); }

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
}
