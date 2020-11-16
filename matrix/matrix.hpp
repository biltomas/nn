#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <iostream>
#include <utility>
#include <vector>

class matrix {
public:
	using value_type = float;

private:
	std::vector<value_type> matrix_;
	std::size_t rows_; // number of rows
	std::size_t cols_; // number of columns

public:
	matrix() {

	}
	matrix(std::size_t rows, std::size_t cols)
	    : matrix_(rows * cols)
	    , rows_(rows)
	    , cols_(cols) {}

	matrix(const std::vector<value_type> &data, std::size_t rows)
	    : matrix_(data.size() % rows == 0
	              ? data
	              : throw std::invalid_argument("Input data has incompatible size"))
	    , rows_(rows)
	    , cols_(data.size() / rows) {}

	// Returns number of rows, number of columns
	std::pair<std::size_t, std::size_t> size() const {
		return { rows_, cols_ };
	}

	// Returns number of rows
	std::size_t rows() const { return rows_; }

	// Returns number of columns
	std::size_t cols() const { return cols_; }
	std::vector<value_type> to_vector() const { return matrix_; }

	/******* TODO - Operators *******/
	// operator[] takes a pair { row index, column index } where
	// - row index is in range [0, rows)
	// - column index is in range [0, cols)

	bool operator==(const matrix& m1) const;


	bool operator!=(const matrix& m1) const;

	
	matrix& operator+=(const matrix& m2);

	matrix& operator-=(const matrix& m2);

	matrix& operator*=(int k);

	value_type& operator[](std::pair<size_t, size_t> index);
	
	const value_type& operator[](std::pair<size_t, size_t> index) const;
};

std::ostream &operator<<(std::ostream &, const matrix &);

/****** TODO - Operators *******/

matrix operator+(matrix m1, const matrix& m2);

matrix operator-(matrix m1, const matrix& m2);

matrix operator*(matrix m1, const int k);

matrix operator*(const int k, matrix m1);
matrix operator*(matrix m1, const matrix& m2);

#endif
