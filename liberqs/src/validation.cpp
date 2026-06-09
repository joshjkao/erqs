#include "validation.hpp"
#include "quantumstate.hpp"
#include <iostream>
#include <memory>
#include <unordered_set>
#include <variant>

bool Validate(const std::shared_ptr<PureState> &pure,
              const ValidationArgs &args) {
  bool ret = true;
  auto handle_error = [&](const std::string &msg) {
    if (args.log_to_stdout) {
      std::cout << msg << "\n";
      if (args.log_state_on_error) {
        Print(pure);
        std::cout << "\n";
      }
    }
    ret = false;
  };

  if (args.check_pure_support) {
    if ((pure->bits & pure->space) != pure->bits) {
      handle_error("pure state bits must be a subset of its subspace");
    }
  }
  return ret;
}

bool Validate(const std::shared_ptr<ProductState> &prod,
              const ValidationArgs &args) {
  bool ret = true;
  auto handle_error = [&](const std::string &msg) {
    if (args.log_to_stdout) {
      std::cout << msg << "\n";
      if (args.log_state_on_error) {
        Print(prod);
        std::cout << "\n";
      }
    }
    ret = false;
  };

  for (const auto &sum : prod->states) {
    if (!Validate(sum, args))
      ret = false;
  }

  if (args.check_no_single_sum_states) {
    if (prod->states.empty()) {
      handle_error("product states must not be empty");
    } else if (prod->states.size() == 1) {
      handle_error("product states holding a single sum state are redundant");
    }
  }
  if (args.check_sum_subspaces) {
    QSpace space{0};
    for (const auto &child : prod->states) {
      if ((space & child->space).any()) {
        handle_error("product states must be over disjoint subspaces");
      }
      space |= child->space;
    }
    if (space != prod->space) {
      handle_error("declared subspace of product state doesn't match the union "
                   "of its children");
    }
  }
  return ret;
}

bool Validate(const std::shared_ptr<SumState> &sum,
              const ValidationArgs &args) {
  bool ret = true;
  auto handle_error = [&](const std::string &msg) {
    if (args.log_to_stdout) {
      std::cout << msg << "\n";
      if (args.log_state_on_error) {
        Print(sum);
        std::cout << "\n";
      }
    }
    ret = false;
  };

  for (const auto &var : sum->states) {
    bool child_ret =
        std::visit([&](const auto &v) { return Validate(v, args); }, var);
    if (!child_ret)
      ret = false;
  }

  if (args.check_no_single_prod_states) {
    if (sum->states.empty()) {
      handle_error("sum states must not be empty");
    } else if (sum->states.size() == 1) {
      if (std::holds_alternative<std::shared_ptr<ProductState>>(
              sum->states[0])) {
        handle_error("sum states holding a single product state are redundant");
      }
    }
  }
  if (args.check_sum_subspaces) {
    QSpace decl_space = sum->space;
    for (const auto &var : sum->states) {
      if (GetSpace(var) != decl_space) {
        handle_error("child space of a sum doesn't match its declared space");
      }
    }
  }
  if (args.check_sum_list_sizes) {
    if (sum->states.size() != sum->coeffs.size()) {
      handle_error(
          "sum states must hold the same number of coefficients and states");
    }
  }
  if (args.check_redundant_pure_states) {
    std::unordered_set<PureState, PureStateHash> pures;
    for (const auto &var : sum->states) {
      if (std::holds_alternative<std::shared_ptr<PureState>>(var)) {
        std::shared_ptr<PureState> pure_ptr =
            std::get<std::shared_ptr<PureState>>(var);
        if (pures.contains(*pure_ptr)) {
          handle_error("sums shouldn't contain duplicate pure states");
        } else {
          pures.insert(*pure_ptr);
        }
      }
    }
  }
  return ret;
}

bool Validate(const KetBra &kb, const ValidationArgs &args) {
  bool log = args.log_to_stdout;
  bool ret = true;
  bool check_ket =
      std::visit([&](const auto &var) { return Validate(var, args); }, kb.ket);
  bool check_bra =
      std::visit([&](const auto &var) { return Validate(var, args); }, kb.bra);
  if (!check_ket | !check_bra)
    ret = false;

  if (args.check_ketbra_subspaces) {
    if ((GetSpace(kb.bra) & GetSpace(kb.ket)).any()) {
      if (log) {
        std::cout << "ketbra ket and bra subspaces must be disjoint";
      }
      ret = false;
    }
  }

  return ret;
}

bool Validate(const SkewOperator &op, const ValidationArgs &args) {
  bool ret = true;
  if (op.empty())
    return ret;
  auto handle_error = [&](const std::string &msg) {
    if (args.log_to_stdout)
      std::cout << msg << "\n";
    ret = false;
  };

  for (const auto &kb : op) {
    if (!Validate(kb, args))
      ret = false;
  }

  if (args.check_skewop_redundant_constants) {
    auto kb_is_constant = [](const KetBra &kb) {
      return (GetSpace(kb.ket) | GetSpace(kb.bra)).none();
    };
    bool has_a_const = false;
    for (const auto &kb : op) {
      if (kb_is_constant(kb) && !has_a_const)
        has_a_const = true;
      else if (kb_is_constant(kb)) {
        handle_error("skewoperator as redundant constant terms");
      }
    }
  }
  if (args.check_skewop_redundant_pures) {
    // todo
  }
  if (args.check_skewop_subspaces) {
    QSpace bspace = op.GetBraSpace();
    QSpace kspace = op.GetKetSpace();
    for (const auto &kb : op) {
      if (GetSpace(kb.bra) != bspace)
        handle_error("skewoperator bra space doesn't match");
      if (GetSpace(kb.ket) != kspace)
        handle_error("skewoperator ket space doesn't match");
    }
  }

  return ret;
}

bool Validate(const PauliOperator &op, const ValidationArgs &args) {
  bool ret = true;
  auto handle_error = [&](const std::string &msg) {
    if (args.log_to_stdout)
      std::cout << msg << "\n";
    ret = false;
  };

  if (args.check_pauliop_disjoint) {
    if ((op.x & op.y).any() || (op.x & op.z).any() || (op.y & op.z).any()) {
      handle_error("pauli operators must be disjoint");
    }
  }

  return ret;
}

bool Validate(const PauliHamiltonian &h, const ValidationArgs &args) {
  bool ret = true;
  auto handle_error = [&](const std::string &msg) {
    if (args.log_to_stdout)
      std::cout << msg << "\n";
    ret = false;
  };

  for (const auto &op : h.operators) {
    if (!Validate(op, args))
      ret = false;
  }

  if (args.check_hamiltonian_list_sizes) {
    if (h.coeffs.size() != h.operators.size()) {
      handle_error("hamiltonians must have an equal number of coefficients and "
                   "operators");
    }
  }

  return ret;
}
