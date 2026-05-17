#include "validation.h"
#include "quantumstate.h"
#include <iostream>
#include <memory>
#include <type_traits>
#include <variant>

static bool validate_sum_helper(const std::shared_ptr<SumState> &sum,
                                const ValidationArgs &args, bool is_root);

bool Validate(const std::shared_ptr<PureState> &pure,
              const ValidationArgs &args) {
  bool log = args.log_to_stdout;
  bool ret = true;
  if (args.check_no_single_sum_states) {
    // unused
  }
  if (args.check_no_single_prod_states) {
    // unused
  }
  if (args.check_product_subspaces) {
    // unused
  }
  if (args.check_sum_subspaces) {
    // unused
  }
  if (args.check_pure_support) {
    if ((pure->bits & pure->space) != pure->bits) {
      if (log) {
        std::cout << "pure state bits must be a subset of its subspace";
        Print(pure);
      }
      ret = false;
    }
  }
  if (args.check_sum_list_sizes) {
    // unused
  }
  return ret;
}

bool Validate(const std::shared_ptr<ProductState> &prod,
              const ValidationArgs &args) {
  bool log = args.log_to_stdout;
  bool ret = true;
  for (const auto &sum : prod->states) {
    if (!validate_sum_helper(sum, args, false))
      ret = false;
  }

  if (args.check_no_single_sum_states) {
    if (prod->states.empty()) {
      if (log) {
        std::cout << "product states must not be empty";
        Print(prod);
      }
      ret = false;
    } else if (prod->states.size() == 1) {
      if (log) {
        std::cout << "product states holding a single sum state are redundant";
        Print(prod);
      }
      ret = false;
    }
  }
  if (args.check_no_single_prod_states) {
    // unused
  }
  if (args.check_product_subspaces) {
    // unused
  }
  if (args.check_sum_subspaces) {
    QSpace space{0};
    for (const auto &child : prod->states) {
      if ((space & child->space).any()) {
        if (log) {
          std::cout << "product states must be over disjoint subspaces";
          ret = false;
        }
      }
      space |= child->space;
    }
    if (space != prod->space) {
      if (log) {
        std::cout
            << "declared subspace of product state doesn't match the union "
               "of its children";
      }
      ret = false;
    }
  }
  if (args.check_pure_support) {
  }
  if (args.check_sum_list_sizes) {
    // unused
  }
  return ret;
}

static bool validate_sum_helper(const std::shared_ptr<SumState> &sum,
                                const ValidationArgs &args, bool is_root) {
  bool log = args.log_to_stdout;
  bool ret = true;
  for (const auto &var : sum->states) {
    bool child_ret =
        std::visit([&](const auto &v) { return Validate(v, args); }, var);
    if (!child_ret)
      ret = false;
  }

  if (args.check_no_single_sum_states) {
    // unused
  }
  if (args.check_no_single_prod_states) {
    if (sum->states.empty()) {
      if (log) {
        std::cout << "sum states must not be empty";
      }
      ret = false;
    } else if (sum->states.size() == 1) {
      if (std::holds_alternative<std::shared_ptr<ProductState>>(
              sum->states[0]) &&
          !is_root) {
        if (log) {
          std::cout
              << "sum states holding a single product state are redundant";
        }
        ret = false;
      }
    }
  }
  if (args.check_product_subspaces) {
    // unused
  }
  if (args.check_sum_subspaces) {
    QSpace decl_space = sum->space;
    for (const auto &var : sum->states) {
      if (GetSpace(var) != decl_space) {
        if (log) {
          std::cout << "child space of a sum doesn't match its declared space";
        }
        ret = false;
      }
    }
  }
  if (args.check_pure_support) {
    // unused
  }
  if (args.check_sum_list_sizes) {
    if (sum->states.size() != sum->coeffs.size()) {
      if (log) {
        std::cout << "sum states must hold the same number of coefficients and "
                     "states";
      }
      ret = false;
    }
  }
  return ret;
}

bool Validate(const std::shared_ptr<SumState> &sum,
              const ValidationArgs &args) {
  return validate_sum_helper(sum, args, true);
}

bool Validate(const SkewOperator &op, const ValidationArgs &args);

bool CheckEqual_slow(const std::shared_ptr<SumState> &sum1,
                     const std::shared_ptr<SumState> &sum2) {
  auto self_inner_1 = Inner_slow(sum1, sum1);
  auto inner_12 = Inner_slow(sum1, sum2);
  return (self_inner_1.real() - inner_12.real()) < 1e-4;
}
