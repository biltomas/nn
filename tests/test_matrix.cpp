#include "catch2.hpp"
#include "../Matrix.hpp"


TEST_CASE("Matrix - basics") {
	std::vector<float> data1 {1, 2,
							  3, 4};
	Matrix<float> matrix1(data1, 2);

	SECTION("Matrix get, set, size") {
		REQUIRE(matrix1.rows() == 2);
		REQUIRE(matrix1.cols() == 2);
		REQUIRE(matrix1[{0, 1}] == 2);
		REQUIRE(matrix1[{1, 1}] == 4);
		matrix1[{1, 1}] = 0;	
		REQUIRE(matrix1[{1, 1}] == 0);
	}

	SECTION("Out of bounds access") {
		REQUIRE_THROWS_AS((matrix1[{2, 1}]), std::out_of_range);
		REQUIRE_THROWS_AS((matrix1[{-1, 0}]), std::out_of_range);
	}
}
