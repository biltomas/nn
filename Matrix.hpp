#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <iostream>
#include <utility>
#include <vector>

template <class T>
class Matrix {
    static_assert(std::is_fundamental<T>::value,
                    "Non-primitive type in Matrix");

private:
	std::vector<T> matrix_;
	std::size_t rows_; // number of rows
	std::size_t cols_; // number of columns
    bool transposed;

public:
	Matrix() {

	}

	Matrix(size_t rows, std::size_t cols)
	    : matrix_(rows * cols)
	    , rows_(rows)
	    , cols_(cols)
        , transposed(false) {}

	Matrix(const std::vector<T> &data, std::size_t rows)
	    : matrix_(data.size() % rows == 0
	              ? data
	              : throw std::invalid_argument("Input data has incompatible size"))
	    , rows_(rows)
	    , cols_(data.size() / rows)
        , transposed(false) {}

	// Returns number of rows, number of columns
	std::pair<size_t, size_t> size() const {
		return { rows_, cols_ };
	}

	// Returns number of rows
	size_t rows() const { return rows_; }

	// Returns number of columns
	size_t cols() const { return cols_; }
	std::vector<T>& to_vector() { return matrix_; }

	/******* TODO - Operators *******/
	// operator[] takes a pair { row index, column index } where
	// - row index is in range [0, rows)
	// - column index is in range [0, cols)

	bool operator==(const Matrix<T>& m1) const;


	bool operator!=(const Matrix<T>& m1) const;
    

    //Random initialization of the entire matrix
    void setRandom();

    //zero-initialization of the entire matrix
	void setZero();

    //randomly initializes the matrix with offset of "number"
	void setNumber(float number);

	
	Matrix<T>& operator+=(const Matrix<T>& m2);

	Matrix<T>& operator-=(const Matrix<T>& m2);

	Matrix<T>& operator*=(T k);


    //transposes the matrix
    //the memory itself remains unchanged to avoid unnecessary
    //operations
	Matrix<T>& transpose();

	T& operator[](std::pair<size_t, size_t> index);
	
	const T& operator[](std::pair<size_t, size_t> index) const;

};

template <typename T>
std::ostream &operator<<(std::ostream &, const Matrix<T> &);

template <typename T>
Matrix<T> operator+(Matrix<T> m1, const Matrix<T>& m2);

template <typename T, typename U>
Matrix<T> operator-(Matrix<T> m1, const Matrix<T>& m2);

//multiplication of matrix with a constant
template <typename T>
Matrix<T> operator*(Matrix<T> m1, const T k);

template <typename T>
Matrix<T> operator*(const T k, Matrix<T> m1);

//multiplication of two matrices
template <typename T>
Matrix<T> operator*(Matrix<T>& m1, Matrix<T>& m2);

#include "Matrix.tpp"

#endif
