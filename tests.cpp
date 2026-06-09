#include "optimization.h"
#include "quantumstate.h"
#include "validation.h"
#include <gtest/gtest.h>
#include <random>
#include <stdexcept>

TEST(Factories, MakeSum) {
  auto pure1 = MakePure("0111", "0000");
  auto pure2 = MakePure("1110", "1110");
  auto pure3 = MakePure("1100", "1100");
  auto pure4 = MakePure("1100", "0000");
  EXPECT_THROW(MakeSum({0.2, 0.2}, {pure1, pure2}), std::runtime_error);
  EXPECT_NO_THROW(MakeSum({2, 5}, {pure3, pure4}));
  EXPECT_THROW(MakeSum({}, {}), std::runtime_error);
}

TEST(Factories, MakeProduct) {
  auto pure1 = MakePure("0011", "0000");
  auto pure2 = MakePure("0011", "0011");
  auto pure3 = MakePure("1100", "1100");
  auto pure4 = MakePure("1100", "0000");
  auto pure5 = MakePure("1010", "1010");
  auto pure6 = MakePure("1010", "0010");
  auto sum1 = MakeSum({5, 6}, {pure1, pure2});
  auto sum2 = MakeSum({7, 8}, {pure3, pure4});
  auto sum3 = MakeSum({9, 1}, {pure5, pure6});
  EXPECT_THROW(MakeProduct({sum1, sum3}), std::runtime_error);
  EXPECT_NO_THROW(MakeProduct({sum1, sum2}));
  EXPECT_THROW(MakeProduct({}), std::runtime_error);
}

TEST(Equality, Pure) {
  auto pure1 = MakePure("1111", "0000");
  auto pure2 = MakePure("1111", "0000");
  auto pure3 = MakePure("1110", "0000");
  auto pure4 = MakePure("1111", "1000");
  EXPECT_TRUE(Equals(pure1, pure1));
  EXPECT_TRUE(Equals(pure1, pure2));
  EXPECT_FALSE(Equals(pure1, pure3));
  EXPECT_FALSE(Equals(pure1, pure4));
}

TEST(Equality, Sum) {
  auto pure1 = MakePure("1111", "0000");
  auto pure2 = MakePure("1111", "1010");
  auto sum1 = MakeSum({complex{2.3, -0.75}, complex{1.2, 4.5}}, {pure1, pure2});
  auto sum2 = Clone(sum1);
  EXPECT_TRUE(Equals_slow(sum1, sum2));
  EXPECT_TRUE(Equals_slow(sum1, sum1));
  sum2 = MakeSum({complex{1.3, -1.75}, complex{1.0, 5.5}}, {pure1, pure2});
  EXPECT_FALSE(Equals_slow(sum1, sum2));
}

TEST(Equality, Literal_flat) {
  auto pure1 = MakePure("1111", "0000");
  auto pure2 = MakePure("1111", "0000");
  auto sum1 = MakeSum({complex{2.3, -0.75}, complex{1.2, 4.5}}, {pure1, pure2});
  EXPECT_TRUE(Equals_literal_flat(sum1, sum1));
}

TEST(Equality, Literal) {
  auto pure1 = MakePure("1111", "0000");
  auto pure2 = MakePure("1111", "1010");
  auto sum1 = MakeSum({1}, {pure1});
  auto sum2 = MakeSum({1}, {pure1});
  EXPECT_TRUE(Equals_literal(sum1, sum2));
  auto sum3 = MakeSum({2}, {pure1});
  EXPECT_FALSE(Equals_literal(sum1, sum3));
  auto sum4 = MakeSum({1, 2, 3}, {pure1, pure1, pure2});
  EXPECT_FALSE(Equals_literal(sum1, sum4));
  auto sum5 = MakeSum({1}, {pure2});
  EXPECT_FALSE(Equals_literal(sum1, sum5));

  std::random_device rd{};
  std::mt19937 gen{rd()};
  auto rand1 = RandomSumState(7, 7, 7, ~QSpace{0}, gen);
  auto rand2 = Clone(rand1);
  EXPECT_TRUE(Equals_literal(rand1, rand2));
}

TEST(FromString, Small) {
  auto pure1 = MakePure("1111", "0101");
  auto pure2 = MakePure("1111", "1010");
  auto sum1 = MakeSum({1.0}, {pure1});
  auto from_str = FromString("[(1.0) 1111 0101]");
  EXPECT_TRUE(Equals_literal(sum1, from_str));
  auto sum2 = MakeSum({2.4, 5.7}, {pure1, pure2});
  auto from_str2 = FromString("[(2.4) 1111 0101 (5.7) 1111 1010]");
  EXPECT_TRUE(Equals_literal(sum2, from_str2));
  auto pure3 = MakePure("110000", "010000");
  auto pure4 = MakePure("110000", "100000");
  auto sum3 = MakeSum({1.1, -1.0}, {pure3, pure4});
  auto prod1 = MakeProduct({sum2, sum3});
  auto pure5 = MakePure("111111", "010101");
  auto sum4 = MakeSum({-10, 3}, {pure5, prod1});
  auto from_str3 =
      FromString("[(-10) 111111 010101 (3) {[(2.4) 1111 0101 (5.7) 1111 1010] "
                 "[(1.1) 110000 010000 (1) 110000 100000]}]");
  auto from_str4 =
      FromString("[(-10) 111111 010101 (3) {[(2.4) 1111 0101 (5.7) 1111 1010] "
                 "[(1.1) 110000 010000 (-1) 110000 100000]}]");
  EXPECT_FALSE(Equals_literal(sum4, from_str3));
  EXPECT_TRUE(Equals_literal(sum4, from_str4));
}

TEST(Flatten, Small) {
  auto pure1 = MakePure("1111", "0101");
  auto pure2 = MakePure("1111", "1010");
  auto sum1 = MakeSum({1.0}, {pure1});
  auto flat1 = Clone(sum1);
  Flatten(flat1);
  auto flat1_from_str = FromString("[(1.0) 1111 0101]");
  EXPECT_TRUE(Equals_literal_flat(flat1, flat1_from_str));

  auto sum2 = MakeSum({2.4, 5.7}, {pure1, pure2});
  auto flat2 = Clone(sum2);
  Flatten(flat2);
  auto flat_from_str = FromString("[(2.4) 1111 0101 (5.7) 1111 1010]");
  EXPECT_TRUE(Equals_literal_flat(flat2, flat_from_str));

  auto pure3 = MakePure("110000", "010000");
  auto pure4 = MakePure("110000", "100000");
  auto sum3 = MakeSum({1.1, -1.0}, {pure3, pure4});
  auto prod1 = MakeProduct({sum2, sum3});
  auto pure5 = MakePure("111111", "010101");
  auto sum4 = MakeSum({-10, 3}, {pure5, prod1});

  auto flat4 = Clone(sum4);
  Flatten(flat4);

  auto flat4_from_str =
      FromString("[(-2.08) 111111 010101 (-7.2) 111111 100101 "
                 "(18.81) 111111 011010 (-17.1) 111111 101010]");

  EXPECT_TRUE(Equals_literal_flat(flat4, flat4_from_str));
}

TEST(Flatten, ComplexCoefficients) {
  // Subspace 1: the first two qubits
  auto pure1 = MakePure("1100", "0100");

  // Subspace 2: the last two qubits (disjoint from pure1)
  auto pure2 = MakePure("0011", "0010");

  // sum1 = (1 + 2i) * pure1
  auto sum1 = MakeSum({{1.0, 2.0}}, {pure1});
  // sum2 = (0 - 1i) * pure2
  auto sum2 = MakeSum({{0.0, -1.0}}, {pure2});

  // Now valid! The spaces "1100" and "0011" are disjoint.
  auto prod = MakeProduct({sum1, sum2});

  auto root_sum = MakeSum({1.0}, {prod});
  auto flat = Clone(root_sum);
  Flatten(flat);

  // Math: 1.0 * (1 + 2i) * (0 - 1i) = 2 - i
  // Tensor product of spaces "1100" and "0011" -> "1111"
  // Tensor product of bits "0100" and "0010" -> "0110"
  auto expected = FromString("[(2.0,-1.0) 1111 0110]");

  EXPECT_TRUE(Equals_literal_flat(flat, expected));
}

TEST(Flatten, TermCancellation) {
  // Both pure states share the exact same subspace "111"
  auto pure1 = MakePure("111", "010");
  auto pure2 = MakePure("111", "101");

  // Valid: All states in the sum belong to subspace "111"
  auto root_sum = MakeSum({5.0, 2.0, -5.0, 3.0}, {pure1, pure2, pure1, pure2});

  auto flat = Clone(root_sum);
  Flatten(flat);

  // pure1 terms perfectly cancel, pure2 terms accumulate to 5.0
  auto expected = FromString("[ (0.0) 111 010 (5.0) 111 101 ]");

  EXPECT_TRUE(Equals_literal_flat(flat, expected));
}

TEST(Flatten, DeepNesting) {
  // Three disjoint subspaces of a 6-bit system
  auto pure1 = MakePure("110000", "010000"); // Qubits 0, 1
  auto pure2 = MakePure("001100", "001000"); // Qubits 2, 3
  auto pure3 = MakePure("000011", "000001"); // Qubits 4, 5

  // Level 1
  auto inner_sum = MakeSum({4.0}, {pure1});

  // Level 2: pure2 (001100) x pure1 (110000) -> 111100
  auto middle_prod = MakeProduct({MakeSum({3.0}, {pure2}), inner_sum});
  auto middle_sum = MakeSum({1.0}, {middle_prod});

  // Level 3: pure3 (000011) x middle (111100) -> 111111
  auto outer_prod = MakeProduct({MakeSum({2.0}, {pure3}), middle_sum});
  auto root_sum = MakeSum({1.0}, {outer_prod});

  auto flat = Clone(root_sum);
  Flatten(flat);

  // Math: 1.0 * (2.0 * (3.0 * 4.0)) = 24.0
  // Bits combined: 010000 | 001000 | 000001 = 011001
  auto expected = FromString("[ (24.0) 111111 011001 ]");

  EXPECT_TRUE(Equals_literal_flat(flat, expected));
}

TEST(FromString, Random) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  auto rand1 = RandomSumState(7, 7, 7, ~QSpace{0}, gen);
  std::string str = Stringify(rand1);
  auto clone = FromString(str);
  std::string strcpy = Stringify(clone);
  EXPECT_EQ(str, strcpy);
  EXPECT_TRUE(Equals_literal(rand1, clone));
}

TEST(Validate, Pure_Support) {
  auto pure1 = MakePure("1111", "0000");
  auto pure2 = MakePure("0000", "0000");
  auto pure3 = MakePure("0011", "0110");
  EXPECT_TRUE(Validate(pure1, ValidationArgs{}));
  EXPECT_TRUE(Validate(pure2, ValidationArgs{}));
  EXPECT_FALSE(Validate(pure3, ValidationArgs{.log_to_stdout = false}));
}

TEST(Validate, Single_Sum_States) {
  auto pure1 = MakePure("001111", "000000");
  auto sum1 = MakeSum({complex{0.5}}, {pure1});
  auto prod1 = MakeProduct({sum1});
  EXPECT_FALSE(Validate(prod1, ValidationArgs{.log_to_stdout = false}));
  auto pure2 = MakePure("110000", "110000");
  auto sum2 = MakeSum({complex{1.5}}, {pure2});
  auto prod2 = MakeProduct({sum1, sum2});
  EXPECT_TRUE(Validate(prod2, ValidationArgs{}));
}

TEST(Validate, Product_Subspaces) {
  auto pure1 = MakePure("001111", "000000");
  auto pure2 = MakePure("110000", "110000");
  auto sum1 = MakeSum({complex{0.5}}, {pure1});
  auto sum2 = MakeSum({complex{1.5}}, {pure2});
  auto prod2 = MakeProduct({sum1, sum2});
  EXPECT_TRUE(Validate(prod2, ValidationArgs{}));
  prod2->states[0]->space = QSpace{"111111"};
  EXPECT_FALSE(Validate(prod2, ValidationArgs{.log_to_stdout = false}));
}

TEST(Validate, Single_Product_States) {
  // single pure state, okay
  auto pure1 = MakePure("001111", "000000");
  auto sum1 = MakeSum({complex{0.5}}, {pure1});
  EXPECT_TRUE(Validate(sum1, ValidationArgs{}));
  // single pure state, okay
  auto pure2 = MakePure("110000", "110000");
  auto sum2 = MakeSum({complex{2}}, {pure2});
  EXPECT_TRUE(Validate(sum2, ValidationArgs{}));
  // single product, bad
  auto prod1 = MakeProduct({sum1, sum2});
  auto sum3 = MakeSum({complex{1}}, {prod1});
  EXPECT_FALSE(Validate(sum3, ValidationArgs{.log_to_stdout = false}));
  // two products, okay
  auto sum4 = MakeSum({complex{1}, complex{2}}, {prod1, prod1});
  EXPECT_TRUE(Validate(sum4, ValidationArgs{}));
}

TEST(Validate, Sum_Subspaces) {
  auto pure1 = MakePure("1111", "1111");
  auto pure2 = MakePure("1111", "0000");
  auto sum1 = MakeSum({complex{1}, complex{2}}, {pure1, pure2});
  EXPECT_TRUE(Validate(sum1, ValidationArgs{}));
  pure1->space = QSpace{"11111"};
  EXPECT_FALSE(Validate(sum1, ValidationArgs{.log_to_stdout = false}));
}

TEST(Validate, Sum_List_Sizes) {
  auto pure1 = MakePure("1111", "1111");
  auto pure2 = MakePure("1111", "0000");
  auto sum1 = MakeSum({complex{1}, complex{2}}, {pure1, pure2});
  EXPECT_TRUE(Validate(sum1, ValidationArgs{}));
  sum1->coeffs.pop_back();
  EXPECT_FALSE(Validate(sum1, ValidationArgs{.log_to_stdout = false}));
  sum1->states.pop_back();
  EXPECT_TRUE(Validate(sum1, ValidationArgs{}));
}

TEST(Validate, Ketbra_Subspaces) {
  auto pure1 = MakePure("01111", "01111");
  auto pure2 = MakePure("10000", "10000");
  KetBra kb{.coeff = complex{1.0}, .ket = pure1, .bra = pure2};
  EXPECT_TRUE(Validate(kb, ValidationArgs{}));
  auto pure3 = MakePure("01100", "01100");
  kb.bra = pure3;
  EXPECT_FALSE(Validate(kb, ValidationArgs{.log_to_stdout = false}));
}

TEST(Validate, Skewop_Redundant_Constants) {
  // todo
}

TEST(Validate, Skewop_Redundant_Pures) {
  // todo
}

TEST(Validate, PauliOp) {
  PauliOperator op1{
      .x = BitString{"1000"}, .y = BitString{"0110"}, .z = BitString{"0001"}};
  EXPECT_TRUE(Validate(op1, CHECK_ALL));
  PauliOperator op2{
      .x = BitString{"1100"}, .y = BitString{"0110"}, .z = BitString{"0001"}};
  EXPECT_FALSE(Validate(op2, CHECK_ALL_QUIET));
}

TEST(Validate, Hamiltonian) {
  std::vector<double> coeffs{1, 2};
  std::vector<PauliOperator> ops{PauliOperator{
      .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0000"}}};
  PauliHamiltonian h1{coeffs, ops};
  EXPECT_FALSE(Validate(h1, CHECK_ALL_QUIET));
  ops.push_back(PauliOperator{
      .x = BitString{"1111"}, .y = BitString{"1111"}, .z = BitString{"1111"}});
  PauliHamiltonian h2{coeffs, ops};
  EXPECT_FALSE(Validate(h2, CHECK_ALL_QUIET));
  ops.pop_back();
  ops.push_back(PauliOperator{
      .x = BitString{"0000"}, .y = BitString{"0001"}, .z = BitString{"1110"}});
  PauliHamiltonian h3{coeffs, ops};
  EXPECT_TRUE(Validate(h3, CHECK_ALL));
}

// We know this is broken, but we should just use Simplify instead
// TEST(RemoveSingles, BasicSum) {
//   auto pure1 = MakePure("1111", "1111");
//   auto sum1 = MakeSum({1}, {pure1});
//   auto prod1 = MakeProduct({sum1});
//   auto sum2 = MakeSum({1}, {prod1});
//   Print(sum2);
//   EXPECT_FALSE(Validate(sum2, ValidationArgs{.log_to_stdout = false}));
//   RemoveSingles(sum2);
//   Print(sum2);
//   EXPECT_TRUE(Validate(sum2, ValidationArgs{}));
// }

TEST(Simplify, Redundant_Pures) {
  auto pure1 = MakePure("00001111", "00001111");
  auto pure2 = MakePure("00001111", "00001111");
  auto sum1 = MakeSum({2, 3}, {pure1, pure2});
  EXPECT_FALSE(Validate(sum1, ValidationArgs{.log_to_stdout = false}));
  auto sum1_simple = Clone(sum1);
  Simplify(sum1_simple);
  EXPECT_TRUE(Validate(sum1_simple, CHECK_ALL));
  EXPECT_TRUE(sum1_simple->coeffs.size() == 1);
  EXPECT_TRUE(Equals_slow(sum1_simple, sum1));
}

TEST(Simplify, SingleStates) {
  auto pure1 = MakePure("1111", "1111");
  auto sum1 = MakeSum({2}, {pure1});
  auto prod1 = MakeProduct({sum1});
  auto sum2 = MakeSum({3}, {prod1});
  auto prod2 = MakeProduct({sum2});
  auto sum3 = MakeSum({4}, {prod2});
  EXPECT_FALSE(Validate(sum3, ValidationArgs{.log_to_stdout = false}));
  Simplify(sum3);
  EXPECT_TRUE(Validate(sum3, ValidationArgs{}));
}

TEST(Simplify, Sum) {
  auto pure1 = MakePure("00001111", "00001111");
  auto pure2 = MakePure("00001111", "00001111");
  auto pure3 = MakePure("11110000", "11110000");
  auto pure4 = MakePure("11110000", "11110000");
  auto pure5 = MakePure("00000000", "00000000");
  auto sum1 = MakeSum({complex{2.1}, complex{3.2}}, {pure1, pure2});
  auto sum2 = MakeSum({complex{0.1}, complex{-3.2, .1}}, {pure3, pure4});
  auto sum3 = MakeSum({complex{2, 3}}, {pure5});
  auto prod1 = MakeProduct({sum1, sum2, sum3});
  auto root = MakeSum({complex{1.0, 0.0}}, {prod1});
  EXPECT_FALSE(Validate(root, ValidationArgs{.log_to_stdout = false}));
  auto root_simple = Clone(root);
  Simplify(root_simple);
  EXPECT_TRUE(Validate(root_simple, ValidationArgs{}));
  EXPECT_TRUE(Equals_slow(root, root_simple));
}

TEST(Simplify, BigTree) {
  std::random_device rd{};
  std::mt19937 gen(rd());
  auto rand1 = RandomSumState(7, 7, 7, ~QSpace{0}, gen);
  auto simp1 = Clone(rand1);
  Simplify(simp1);
  EXPECT_TRUE(Validate(simp1, CHECK_ALL));
}

TEST(Inner, Slow) {
  auto pure1 = MakePure("00001111", "00001111");
  auto pure2 = MakePure("00001111", "00000111");
  auto sum1 = MakeSum({0.4, 10.1}, {pure1, pure2});
  auto sum2 = MakeSum({0.4, 10.1}, {pure1, pure2});
  EXPECT_EQ(Inner_slow(sum1, sum2), 0.4 * 0.4 + 10.1 * 10.1);
  auto sum3 = MakeSum({0.4, 10.1}, {pure1, pure1});
  auto sum4 = MakeSum({0.4, 10.1}, {pure1, pure1});
  EXPECT_EQ(Inner_slow(sum3, sum4).real(), (0.4 + 10.1) * (0.4 + 10.1));
}
