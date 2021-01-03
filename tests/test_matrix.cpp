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


TEST_CASE("Matrix - unary operators") {
    std::vector<float> data {1, 2, 3,
                             4, 5, 6};
    Matrix<float> matrix(data, 2);
    SECTION("Transposition") {
        std::vector<float> expected_data {1, 4,
                                          2, 5,
                                          3, 6};
        Matrix<float> expected_matrix (expected_data, 3);
        auto result = matrix.transpose();
        REQUIRE(result.rows() == 3);
        REQUIRE(result.cols() == 2);
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 2; col++) {
                CHECK(result[{row, col}] == expected_matrix[{row, col}]);
            }
        }
    }
}


TEST_CASE("Matrix - binary operators") {
    std::vector<float> data1 {1, 2,
                              3, 4,
                              5, 6};
    std::vector<float> data2 {1, 2, 3,
                              4, 5, 6};
    std::vector<float> data3 {1, 2,
                              3, 4};
    Matrix<float> mat1(data1, 3);
    Matrix<float> mat2(data2, 2);
    Matrix<float> mat3(data3, 2);

    SECTION("Matrix multiplication - basic") {
        const std::vector<float> expected_data {9, 12, 15,
                                                19, 26, 33,
                                                29, 40, 51};
        const Matrix<float> expected_matrix(expected_data, 3);
        Matrix<float> result(3, 3);
        matmul(mat1, mat2, result);
        REQUIRE(result.rows() == 3);
        REQUIRE(result.cols() == 3);
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 2; col++) {
                CHECK(result[{row, col}] == expected_matrix[{row, col}]);
            }
        }
    }
    SECTION("Matrix multiplication - undefined multiplication") {
        Matrix<float> result(3, 3);
        REQUIRE_THROWS_AS(matmul(mat2, mat3, result), std::invalid_argument);
    }
}
