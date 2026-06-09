#pragma once
#include "hamiltonian.hpp"
#include "quantumstate.hpp"
#include "skewoperator.hpp"

struct ValidationArgs {
  bool log_to_stdout = true; // write to std::cout if an error is found
  bool log_state_on_error = true;

  bool check_pure_support = true; // bits are a subset of space

  bool check_no_single_prod_states = true; // single states
  bool check_product_subspaces = true;     // spaces must be disjoint

  bool check_no_single_sum_states = true; // single child states
  bool check_sum_subspaces = true;        // spaces must be identical
  bool check_sum_list_sizes = true;       // same number of coeffients as states
  bool check_redundant_pure_states = true;

  bool check_ketbra_subspaces = true; // spaces must be disjoint

  bool check_skewop_redundant_constants = true;
  bool check_skewop_redundant_pures = true;
  bool check_skewop_subspaces = true;

  bool check_pauliop_disjoint = true;
  bool check_hamiltonian_list_sizes = true;
};

inline constinit ValidationArgs CHECK_ALL{};

inline constinit ValidationArgs CHECK_ALL_QUIET{.log_to_stdout = false,
                                                .log_state_on_error = false};

inline constinit ValidationArgs CHECK_CORRECTNESS_ONLY{
    .check_no_single_prod_states = false,
    .check_no_single_sum_states = false,
    .check_redundant_pure_states = false,
    .check_skewop_redundant_constants = false,
    .check_skewop_redundant_pures = false};

bool Validate(const std::shared_ptr<PureState> &pure,
              const ValidationArgs &args);
bool Validate(const std::shared_ptr<ProductState> &prod,
              const ValidationArgs &args);
bool Validate(const std::shared_ptr<SumState> &sum, const ValidationArgs &args);
bool Validate(const KetBra &kb, const ValidationArgs &args);
bool Validate(const SkewOperator &op, const ValidationArgs &args);
bool Validate(const PauliOperator &p, const ValidationArgs &args);
bool Validate(const PauliHamiltonian &h, const ValidationArgs &args);
