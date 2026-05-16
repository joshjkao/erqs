#include "quantumstate.h"
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
	EXPECT_EQ(*pure1, *pure1);
	EXPECT_EQ(*pure1, *pure2);
	EXPECT_NE(*pure1, *pure3);
	EXPECT_NE(*pure1, *pure4);
}

TEST(Equality, Product) {
	auto pure1 = MakePure("0011", "0000");
	auto pure2 = MakePure("0011", "0011");
	auto pure3 = MakePure("1100", "1100");
	auto pure4 = MakePure("1100", "0000");
	auto pure5 = MakePure("1010", "1010");
	auto pure6 = MakePure("1010", "0010");
	auto pure7 = MakePure("0101", "0001");
	auto pure8 = MakePure("0101", "0100");
	auto sum1 = MakeSum({5, 6}, {pure1, pure2});
	auto sum2 = MakeSum({7, 8}, {pure3, pure4});
	auto sum3 = MakeSum({9, 1}, {pure5, pure6});
	auto sum4 = MakeSum({11, 13}, {pure7, pure8});
	auto prod1 = MakeProduct({sum1, sum2});
	auto prod2 = MakeProduct({sum3, sum4});
	EXPECT_FALSE(Equals(*prod1, *pure1));
	EXPECT_NE(*prod1, *prod2);
	EXPECT_EQ(*prod1, *prod1);
}

TEST(Equality, Sum) {
	auto pure1 = MakePure("0011", "0000");
	auto pure2 = MakePure("0011", "0011");
	auto pure3 = MakePure("1100", "1100");
	auto pure4 = MakePure("1100", "0000");
	auto pure5 = MakePure("1010", "1010");
	auto pure6 = MakePure("1010", "0010");
	auto pure7 = MakePure("0101", "0001");
	auto pure8 = MakePure("0101", "0100");
	auto sum1 = MakeSum({5, 6}, {pure1, pure2});
	auto sum2 = MakeSum({7, 8}, {pure3, pure4});
	auto sum3 = MakeSum({9, 1}, {pure5, pure6});
	auto sum4 = MakeSum({11, 13}, {pure7, pure8});
	EXPECT_NE(sum1, sum2);
	EXPECT_EQ(sum1, sum1);
	EXPECT_EQ(sum3, sum3);
}

TEST(Equality, KetBra) {
	auto pure1 = MakePure("0011", "0000");
	auto pure2 = MakePure("0011", "0011");
	auto pure3 = MakePure("1100", "1100");
	auto pure4 = MakePure("1100", "0000");
	auto pure5 = MakePure("1010", "1010");
	auto pure6 = MakePure("1010", "0010");
	auto pure7 = MakePure("0101", "0001");
	auto pure8 = MakePure("0101", "0100");
	auto sum1 = MakeSum({1, 2}, {pure1, pure2});
	auto sum2 = MakeSum({0.5, 1}, {pure1, pure2});
	auto prod1 = MakeProduct({sum1});
	auto prod2 = MakeProduct({sum2});
	auto kb1 = KetBra{.coeff = 1, .ket = pure1, .bra = pure3};
	auto kb2 = KetBra{.coeff = 1, .ket = pure1, .bra = pure3};
	auto kb3 = KetBra{.coeff = 1, .ket = prod1, .bra = pure3};
	auto kb4 = KetBra{.coeff = 2, .ket = prod2, .bra = pure3};
	EXPECT_EQ(kb1, kb2);
	EXPECT_NE(kb1, kb3);
	// EXPECT_EQ(kb3, kb4);
}

// TEST(Equality, SkewOperator) {
// 	// this test is less trivial because ketbra terms can be
// 	// in any order and searching for them might be pretty expensive
// 	// (exponential operator per equality test)
// }

TEST(Inner, PurePure) {
	auto pure1 = MakePure("1100", "0000");
	auto pure2 = MakePure("0110", "0010");
	auto pure3 = MakePure("1100", "0000");
	auto pure4 = MakePure("1100", "1000");
	auto pure5 = MakePure("1110", "1000");
	KetBra kb{.coeff = 1.0,
	          .ket = MakePure("0010", "0010"),
	          .bra = MakePure("1000", "0000")};
	KetBra one{.coeff = 1.0,
	           .ket = MakePure("0000", "0000"),
	           .bra = MakePure("0000", "0000")};
	KetBra zero{.coeff = 0.0,
	            .ket = MakePure("0000", "0000"),
	            .bra = MakePure("0000", "0000")};
	SkewOperator inner12 = Inner(pure1, pure2);
	SkewOperator inner13 = Inner(pure1, pure3);
	SkewOperator inner14 = Inner(pure1, pure4);
	SkewOperator inner15 = Inner(pure1, pure5);
	EXPECT_TRUE(inner12.ketbras.size() == 1);
	EXPECT_TRUE(inner13.ketbras.size() == 1);
	EXPECT_TRUE(inner14.ketbras.size() == 1);
	EXPECT_TRUE(inner15.ketbras.size() == 1);
	EXPECT_EQ(inner12.ketbras[0], kb);
	EXPECT_EQ(inner13.ketbras[0], one);
	EXPECT_EQ(inner14.ketbras[0], zero);
	EXPECT_EQ(inner15.ketbras[0], zero);
}

TEST(Inner, SumSum) {
	// completely joint
	auto pure1 = MakePure("1100", "0000");
	auto pure2 = MakePure("1100", "1100");
	auto pure3 = MakePure("1100", "0100");
	auto pure4 = MakePure("1100", "1000");
	auto sum1 = MakeSum({2, 3, 5, 7}, {pure1, pure2, pure3, pure4});
	SkewOperator in1 = Inner(sum1, sum1);
	in1 = Simplify(in1);
	// EXPECT_TRUE(Simplify(in1));
	EXPECT_EQ(in1.ketbras.size(), 1uz);
	EXPECT_EQ(std::real(in1.ketbras[0].coeff), 87);

	// completely disjoint
	auto pure5 = MakePure("0011", "0000");
	auto pure6 = MakePure("0011", "0011");
	auto pure7 = MakePure("0011", "0001");
	auto pure8 = MakePure("0011", "0010");
	auto sum2 = MakeSum({2, 3, 5, 7}, {pure5, pure6, pure7, pure8});
	SkewOperator in12 = Inner(sum1, sum2);
	in12 = Simplify(in12);
	// EXPECT_FALSE(Simplify(in12));
	EXPECT_EQ(in12.ketbras.size(), 17uz);

	// partially joint
	auto pure9 = MakePure("0110", "0000");
	auto pure10 = MakePure("0110", "0110");
	auto pure11 = MakePure("0110", "0010");
	auto pure12 = MakePure("0110", "0100");
	auto sum3 = MakeSum({2, 3, 5, 7}, {pure9, pure10, pure11, pure12});
	SkewOperator in13 = Inner(sum1, sum3);
	in13 = Simplify(in13);
	// EXPECT_FALSE(Simplify(in13));
	EXPECT_EQ(in13.ketbras.size(), 9uz);
}

TEST(Inner, ProductProduct) {
	// completely joint
	auto pure1 = MakePure("1100", "0000");
	auto pure2 = MakePure("1100", "1100");
	auto pure3 = MakePure("0011", "0001");
	auto pure4 = MakePure("0011", "0010");
	auto sum1 = MakeSum({2, 3}, {pure1, pure2});
	auto sum2 = MakeSum({5, 7}, {pure3, pure3});
	auto prod1 = MakeProduct({sum1, sum2});
	SkewOperator in1 = Inner(prod1, prod1);
	in1 = Simplify(in1);
	// EXPECT_TRUE(Simplify(in1));
	EXPECT_EQ(in1.ketbras.size(), 1uz);

	// completely disjoint
	auto pure5 = MakePure("1000", "0000");
	auto pure6 = MakePure("1000", "1000");
	auto pure7 = MakePure("0100", "0100");
	auto pure8 = MakePure("0100", "0000");
	auto pure9 = MakePure("0010", "0000");
	auto pure10 = MakePure("0010", "0010");
	auto pure11 = MakePure("0001", "0001");
	auto pure12 = MakePure("0001", "0000");
	auto sum3 = MakeSum({2, 3}, {pure5, pure6});
	auto sum4 = MakeSum({5, 7}, {pure7, pure8});
	auto prod2 = MakeProduct({sum3, sum4});
	auto sum5 = MakeSum({11, 13}, {pure9, pure10});
	auto sum6 = MakeSum({17, 19}, {pure11, pure12});
	auto prod3 = MakeProduct({sum5, sum6});
	SkewOperator in23 = Inner(prod2, prod3);
	in23 = Simplify(in23);
	// EXPECT_FALSE(Simplify(in23));
	EXPECT_EQ(in23.ketbras.size(), 2uz);
	// Print(in23);

	auto sum7 = MakeSum({23, 29}, {prod2, prod3});
	auto in77 = Inner(sum7, sum7);
	in77 = Simplify(in77);
	// EXPECT_FALSE(Simplify(in77));
	Print(in77);
}

TEST(Inner, Random) {
	std::random_device rd{};
	std::mt19937 gen{rd()};
	auto root = RandomState(3, gen);
	Print(root);
	SkewOperator in = Inner(root, root);
	in = Simplify(in);
	// EXPECT_TRUE(Simplify(in));
}
