#pragma once
#include "quantumstate.h"

struct ValidationArgs {
  bool log_to_stdout = true;
  bool check_no_single_sum_states = true;
  bool check_no_single_prod_states = true;
  bool check_product_subspaces = true;
  bool check_sum_subspaces = true;
  bool check_pure_support = true;
  bool check_sum_list_sizes = true;
};

bool Validate(const std::shared_ptr<PureState> &pure,
              const ValidationArgs &args);
bool Validate(const std::shared_ptr<ProductState> &prod,
              const ValidationArgs &args);
bool Validate(const std::shared_ptr<SumState> &sum, const ValidationArgs &args);
bool Validate(const SkewOperator &op, const ValidationArgs &args);

bool CheckEqual_slow(const std::shared_ptr<SumState> &sum1,
                     const std::shared_ptr<SumState> &sum2);
