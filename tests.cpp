#include "optimization.h"
#include "quantumstate.h"
#include "validation.h"
#include <gtest/gtest.h>

TEST(Factories, MakeSum) {
  auto pure1 = MakePure("0111", "0000");
  auto pure2 = MakePure("1110", "1110");
  auto pure3 = MakePure("1100", "1100");
  auto pure4 = MakePure("1100", "0000");
  EXPECT_THROW(MakeSum({0.2, 0.2}, {pure1, pure2}), std::runtime_error);
  EXPECT_NO_THROW(MakeSum({2, 5}, {pure3, pure4}));
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
  auto sum1_simple = Simplify(sum1);
  EXPECT_TRUE(Validate(sum1_simple, ValidationArgs{}));
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
  sum3 = Simplify(sum3);
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
  auto root_simple = Simplify(root);
  EXPECT_TRUE(Validate(root_simple, ValidationArgs{}));
  EXPECT_TRUE(Equals_slow(root, root_simple));
}

TEST(Inner, Slow) {
  auto pure1 = MakePure("00001111", "00001111");
  auto pure2 = MakePure("00001111", "00000111");
  auto sum1 = MakeSum({0.4, 10.1}, {pure1, pure2});
  auto sum2 = MakeSum({0.4, 10.1}, {pure1, pure2});
  EXPECT_EQ(Inner_slow(sum1, sum2), 0.4 * 0.4 + 10.1 * 10.1);
  auto sum3 = MakeSum({0.4, 10.1}, {pure1, pure1});
  auto sum4 = MakeSum({0.4, 10.1}, {pure1, pure1});
  EXPECT_EQ(Inner_slow(sum3, sum4), (0.4 + 10.1) * (0.4 + 10.1));
}
