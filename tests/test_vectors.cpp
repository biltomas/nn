#define CATCH_CONFIG_MAIN
#include "catch2.hpp"
#include "../RowVector.hpp"

#define TEST_ROWVEC_SIZE 10

TEST_CASE("RowVector - get and set") {
    RowVector<float> test_vector(TEST_ROWVEC_SIZE);

    SECTION("Basics - get, length, set") {
        REQUIRE(test_vector.length() == 10);
        for (int i = 0; i < TEST_ROWVEC_SIZE; i++) {
            test_vector.setValue(i, i);
        }
        for (int i = 0; i < TEST_ROWVEC_SIZE; i++) {
            REQUIRE(test_vector.coeffRef(i) == i);
        }
    }

    SECTION("Out of range get") {
        REQUIRE_THROWS_AS(test_vector.coeffRef(TEST_ROWVEC_SIZE), std::out_of_range);
        REQUIRE_THROWS_AS(test_vector.coeffRef(-1), std::out_of_range);
    }
}


TEST_CASE("RowVector - operators") {
    std::vector<float> vec1 {1, 2, 3, 4};
    std::vector<float> vec2 {4, 3, 2, 1};
    std::vector<float> vec3 {1, 2};
    RowVector<float> row_vec_1(vec1);
    RowVector<float> row_vec_2(vec2);
    RowVector<float> row_vec_3(vec3);
    SECTION("Dot product") {
        REQUIRE(row_vec_1.dot(row_vec_2) == 20);
        REQUIRE_THROWS_AS(row_vec_1.dot(row_vec_3), std::out_of_range);
    }

    SECTION("Matrix Multiplication") {
        std::vector<float> mat {1, 1,
                                2, 3};
        Matrix<float> matrix(mat, 2);
        std::vector<float> vector {5, 1};
        RowVector<float> row_vec(vector);
        RowVector<float> result(2);
        matmul(row_vec, matrix, result);
        REQUIRE(result.length() == 2);
        std::vector<float> expected_result {7, 8};
        for (int i = 0; i < 2; i++) {
            CHECK(expected_result[i] == result.coeffRef(i));
        }
    }

    SECTION("Substract") {
        auto result = row_vec_1 - row_vec_2;
        REQUIRE(result.length() == 4);
        const std::vector<float> expected {-3, -1, 1, 3};
        for (int i = 0; i < 4; i++) {
            CHECK(expected[i] == result.coeffRef(i));
        }
        REQUIRE_THROWS_AS(row_vec_1 - row_vec_3, std::invalid_argument);
    }
}
