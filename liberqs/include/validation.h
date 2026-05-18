#pragma once
#include "quantumstate.h"

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
};

bool Validate(const std::shared_ptr<PureState> &pure,
              const ValidationArgs &args);
bool Validate(const std::shared_ptr<ProductState> &prod,
              const ValidationArgs &args);
bool Validate(const std::shared_ptr<SumState> &sum, const ValidationArgs &args);
bool Validate(const KetBra &kb, const ValidationArgs &args);
bool Validate(const SkewOperator &op, const ValidationArgs &args);
