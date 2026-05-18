#include "optimization.h"
#include "quantumstate.h"
#include "validation.h"
#include <iostream>
#include <random>

int main() {
  size_t max_depth = 5, n_terms = 2, n_factors = 1;
  bool test_equality = true;
  QSpace space = ~QSpace{0};

  std::random_device rd{};
  std::mt19937 gen{rd()};

  auto random_state = RandomSumState(max_depth, n_terms, n_factors, space, gen);
  auto original = Clone(random_state);

  Print(original);

  auto simplified = Simplify(random_state);
  bool simplified_is_valid = Validate(simplified, ValidationArgs{});

  bool simplified_equals_original = false;
  if (test_equality)
    simplified_equals_original = Equals_slow(simplified, original);

  std::cout << "{\n";
  std::cout << "\t\"is_valid\": " << simplified_is_valid << ",\n";
  std::cout << "\t\"matches_original\": " << simplified_equals_original << "\n";
  std::cout << "}";
}
